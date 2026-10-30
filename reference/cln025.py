"""CLN025 folding campaign on GROMACS: build, equilibrate, chunked production with resume, QC, analysis.

Used by the notebooks in <PROJECT>/notebooks/. Settings follow the review doc "Modern MD Practice" (Sections 1, 2, 8).
Every GROMACS call goes through `sh`, which raises with the tail of the output when a command fails.
"""
from __future__ import annotations

import glob
import json
import os
import re
import shutil
import subprocess

import numpy as np

FF_REPO = "https://github.com/paulrobustelli/Force-Fields.git"
FF_COMMIT = "25e729d64065aa4f5f61dc5e9036f8e6dd070a36"   # pinned: the commit this module was tested with
SEQUENCE = "YYDPETGTWY"

# Section 8.1 of the review. Stage dicts below override these keys.
BASE_MDP = {
    "integrator": "md", "dt": 0.002, "nsteps": 0,
    "cutoff-scheme": "Verlet", "verlet-buffer-tolerance": 0.005,
    "coulombtype": "PME", "rcoulomb": 1.0, "rvdw": 1.0, "DispCorr": "EnerPres", "fourierspacing": 0.12,
    "constraints": "h-bonds", "constraint-algorithm": "LINCS", "lincs-order": 4, "lincs-iter": 1,
    "tcoupl": "v-rescale", "tc-grps": "System", "tau-t": 0.5, "ref-t": 340,
    "pcoupl": "C-rescale", "tau-p": 2.0, "ref-p": 1.0, "compressibility": 4.5e-5,
    "comm-mode": "Linear", "nstcomm": 100,
    "nstxout-compressed": 5000, "nstenergy": 5000, "nstlog": 5000, "nstcalcenergy": 100,
    "continuation": "yes", "gen-vel": "no",
}
STAGES = {
    "em": {"integrator": "steep", "emtol": 500, "nsteps": 20000, "tcoupl": "no", "pcoupl": "no",
           "continuation": "no", "nstxout-compressed": 0, "nstenergy": 500},
    "nvt": {"nsteps": 50000, "pcoupl": "no", "continuation": "no", "gen-vel": "yes", "define": "-DPOSRES",
            "nstxout-compressed": 0},
    "npt_posres": {"nsteps": 100000, "define": "-DPOSRES", "refcoord-scaling": "com", "nstxout-compressed": 0},
    "npt_free": {"nsteps": 100000, "nstxout-compressed": 0},
}


def sh(cmd, cwd, env=None, stdin=None, timeout=None) -> str:
    """Run a command; return stdout+stderr; raise with the output tail if it fails."""
    r = subprocess.run(cmd, cwd=cwd, env=env, input=stdin, capture_output=True, text=True, timeout=timeout,
                       shell=isinstance(cmd, str))
    out = r.stdout + r.stderr
    if r.returncode != 0:
        raise RuntimeError(f"Command failed ({r.returncode}) in {cwd}: {cmd}\n{out[-3000:]}")
    return out


def write_mdp(path, stage=None, **overrides):
    params = {**BASE_MDP, **(STAGES[stage] if stage else {}), **overrides}
    with open(path, "w") as f:
        f.writelines(f"{k:<24} = {v}\n" for k, v in params.items())
    return params


# ------------------------------------------------------------------ inputs

def fetch_forcefield(dest_dir: str) -> str:
    """Clone the pinned Robustelli force-field repo once; return the path to a99SBdisp.ff."""
    repo = os.path.join(dest_dir, "Force-Fields")
    if not os.path.isdir(repo):
        sh(["git", "clone", "-q", FF_REPO, repo], cwd=dest_dir)
        sh(["git", "checkout", "-q", FF_COMMIT], cwd=repo)
    return os.path.join(repo, "Gromacs_FFs", "a99SBdisp.ff")


def clean_crystal(pdb_text: str, chain: str = "A") -> str:
    """Keep protein ATOM records of one chain, first altloc only (' ' or 'A'); drop waters and HETATM."""
    keep = [l for l in pdb_text.splitlines()
            if l.startswith("ATOM") and l[21] == chain and l[16] in (" ", "A")]
    seq = []
    for l in keep:
        key = (l[22:27], l[17:20])
        if not seq or seq[-1] != key:
            seq.append(key)
    three = {"TYR": "Y", "ASP": "D", "PRO": "P", "GLU": "E", "THR": "T", "GLY": "G", "TRP": "W"}
    found = "".join(three.get(r, "?") for _, r in seq)
    assert found == SEQUENCE, f"Chain {chain} has sequence {found}, expected {SEQUENCE}"
    return "\n".join(l[:16] + " " + l[17:] for l in keep) + "\nEND\n"


def build_extended(path: str) -> str:
    import PeptideBuilder
    from Bio.PDB import PDBIO
    io = PDBIO()
    io.set_structure(PeptideBuilder.make_extended_structure(SEQUENCE))
    io.save(path)
    return path


# ------------------------------------------------------------------ build and equilibrate

def build_system(start_pdb: str, wd: str, ff_dir: str, gmx: str, env=None, box_nm: float = 6.0,
                 salt_M: float = 0.15) -> dict:
    """pdb2gmx -> fixed dodecahedron box -> solvate -> ions. Same box for every start structure."""
    os.makedirs(wd, exist_ok=True)
    if not os.path.isdir(os.path.join(wd, "a99SBdisp.ff")):
        shutil.copytree(ff_dir, os.path.join(wd, "a99SBdisp.ff"))
    shutil.copy(start_pdb, os.path.join(wd, "start.pdb"))
    run = lambda args, stdin=None: sh([gmx] + args, wd, env, stdin)   # noqa: E731
    run(["pdb2gmx", "-f", "start.pdb", "-o", "prot.gro", "-p", "topol.top", "-ff", "a99SBdisp",
         "-water", "select", "-ignh"], stdin="1\n")
    run(["editconf", "-f", "prot.gro", "-o", "box.gro", "-c", "-bt", "dodecahedron", "-box", str(box_nm)])
    run(["solvate", "-cp", "box.gro", "-cs", "a99SBdisp.ff/a99SBdisp_water.gro", "-o", "solv.gro", "-p", "topol.top"])
    write_mdp(os.path.join(wd, "ions.mdp"), "em")
    run(["grompp", "-f", "ions.mdp", "-c", "solv.gro", "-p", "topol.top", "-o", "ions.tpr", "-maxwarn", "1"])
    run(["genion", "-s", "ions.tpr", "-o", "ions.gro", "-p", "topol.top", "-pname", "NA", "-nname", "CL",
         "-neutral", "-conc", str(salt_M)], stdin="SOL\n")
    natoms = int(open(os.path.join(wd, "ions.gro")).readlines()[1])
    return {"gro": os.path.join(wd, "ions.gro"), "top": os.path.join(wd, "topol.top"), "natoms": natoms}


def equilibrate(wd: str, gmx: str, env=None, temp: float = 340, seed: int = 1, mdrun_extra=()) -> str:
    """EM -> NVT (posres) -> NPT (posres) -> NPT (free). Returns the path of the final .gro; .cpt next to it."""
    prev, cpt = "ions.gro", None
    for stage in ("em", "nvt", "npt_posres", "npt_free"):
        extra = {"ref-t": temp, "gen-temp": temp, "gen-seed": seed} if stage == "nvt" else {"ref-t": temp}
        write_mdp(os.path.join(wd, f"{stage}.mdp"), stage, **extra)
        args = ["grompp", "-f", f"{stage}.mdp", "-c", prev, "-r", prev, "-p", "topol.top", "-o", f"{stage}.tpr"]
        if cpt:
            args += ["-t", cpt]
        sh([gmx] + args, wd, env)
        sh([gmx, "mdrun", "-deffnm", stage] + list(mdrun_extra if stage != "em" else ()), wd, env)
        prev, cpt = f"{stage}.gro", (f"{stage}.cpt" if stage != "em" else None)
    return os.path.join(wd, prev)


def make_production_tpr(src_wd: str, rep_dir: str, gmx: str, env=None, ns: float = 1000, temp: float = 340,
                        seed: int = 1, hmr: bool = False, mdp_extra=None) -> str:
    """Production tpr in rep_dir from the equilibrated system in src_wd. New velocities per seed (independent replica).

    HMR (factor 3, 4 fs) only when asked: Section 1.3 says not for runs whose output is a rate.
    """
    src_wd, rep_dir = os.path.abspath(src_wd), os.path.abspath(rep_dir)
    os.makedirs(rep_dir, exist_ok=True)
    dt = 0.004 if hmr else 0.002
    extra = {"nsteps": int(round(ns * 1000 / dt)), "dt": dt, "ref-t": temp, "continuation": "no",
             "gen-vel": "yes", "gen-temp": temp, "gen-seed": seed, **(mdp_extra or {})}
    if hmr:
        extra["mass-repartition-factor"] = 3
    write_mdp(os.path.join(rep_dir, "md.mdp"), **extra)
    for name in ("topol.top", "posre.itp"):
        if os.path.exists(os.path.join(src_wd, name)):
            shutil.copy(os.path.join(src_wd, name), rep_dir)
    link = os.path.join(rep_dir, "a99SBdisp.ff")
    if not os.path.exists(link):
        shutil.copytree(os.path.join(src_wd, "a99SBdisp.ff"), link)
    sh([gmx, "grompp", "-f", "md.mdp", "-c", os.path.join(src_wd, "npt_free.gro"), "-p", "topol.top",
        "-o", "md.tpr"], rep_dir, env)
    return os.path.join(rep_dir, "md.tpr")


# ------------------------------------------------------------------ production in chunks

def mdrun_command(gmx: str, rep_dirs: list, maxh: float, ntomp: int, gpu: bool, plumed: str = None) -> list:
    """One replica: plain mdrun. Several: one MPI job with -multidir (all replicas share the GPU)."""
    cmd = [gmx, "mdrun", "-deffnm", "md", "-cpi", "md.cpt", "-noappend", "-maxh", str(maxh), "-ntomp", str(ntomp)]
    if gpu:
        cmd += ["-nb", "gpu", "-pme", "gpu", "-bonded", "gpu"] + ([] if plumed else ["-update", "gpu"])  # PLUMED needs CPU update
    if plumed:
        cmd += ["-plumed", plumed]
    if len(rep_dirs) == 1:
        return cmd
    return ["mpirun", "-np", str(len(rep_dirs)), "--oversubscribe"] + cmd + ["-multidir"] + list(rep_dirs)


def done_ns(rep_dir: str) -> float:
    """Simulated time (ns) at the latest checkpoint, 0 if none: last 'Writing checkpoint, step N' x dt."""
    dt = float(re.search(r"^dt\s*=\s*(\S+)", open(os.path.join(rep_dir, "md.mdp")).read(), re.M).group(1))
    steps = [int(x) for log in glob.glob(os.path.join(rep_dir, "md*.log"))
             for x in re.findall(r"Writing checkpoint, step (\d+)", open(log, errors="ignore").read())]
    return max(steps, default=0) * dt / 1000.0


def sync(src: str, dst: str, patterns=("md*.xtc", "md*.edr", "md*.log", "md.cpt", "md_prev.cpt", "md.tpr",
                                       "md.mdp", "topol.top", "*.dat", "COLVAR*", "KERNELS*")):
    """Copy files whose size changed. Parts are separate files (-noappend), so each big file is copied once."""
    os.makedirs(dst, exist_ok=True)
    copied = []
    for pat in patterns:
        for f in glob.glob(os.path.join(src, pat)):
            target = os.path.join(dst, os.path.basename(f))
            if not os.path.exists(target) or os.path.getsize(target) != os.path.getsize(f):
                shutil.copy2(f, target + ".tmp")
                os.replace(target + ".tmp", target)          # never leave a half-written file on Drive
                copied.append(os.path.basename(f))
    return copied


def restore_replica(drive_dir: str, local_dir: str, skip=(".xtc",)) -> list:
    """Copy a replica back from Drive to resume it. Trajectory parts stay on Drive: -noappend does not need them."""
    os.makedirs(local_dir, exist_ok=True)
    got = []
    for f in sorted(os.listdir(drive_dir)):
        src = os.path.join(drive_dir, f)
        if os.path.isfile(src) and not f.endswith(skip) and not os.path.exists(os.path.join(local_dir, f)):
            shutil.copy2(src, local_dir)
            got.append(f)
    return got


# ------------------------------------------------------------------ arm B: OPES multithermal (review Section 6)

def write_opes_multithermal(rep_dir: str, temp_max: float, pace: int = 500, stride: int = 500,
                            restart: bool = False) -> str:
    """PLUMED input: OPES_EXPANDED over a temperature range from the thermostat temperature up to temp_max.

    Needs the PLUMED 'opes' module (default-off in PLUMED 2.10). STATE lets a resumed chunk continue the bias.
    """
    lines = (["RESTART"] if restart else []) + [
        "ene: ENERGY",
        f"ecv: ECV_MULTITHERMAL ARG=ene TEMP_MAX={temp_max}",
        f"opes: OPES_EXPANDED ARG=ecv.* PACE={pace} FILE=DELTAFS STATE_WFILE=STATE STATE_WSTRIDE={100 * pace}"
        + (" STATE_RFILE=STATE" if restart else ""),
        f"PRINT STRIDE={stride} ARG=ene,opes.bias FILE=COLVAR FMT=%.6f",
    ]
    path = os.path.join(rep_dir, "plumed.dat")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def read_colvar(path: str):
    """PLUMED COLVAR as a DataFrame. A restart repeats the header; those lines are skipped."""
    import pandas as pd
    names = open(path).readline().split()[2:]
    rows = [l.split() for l in open(path) if not l.startswith("#") and l.strip()]
    df = pd.DataFrame(np.array(rows, dtype=float), columns=names)
    return df.drop_duplicates("time", keep="last").reset_index(drop=True)


def multithermal_logweights(ene, bias, temp0: float, temp: float):
    """log-weights that reweight an OPES multithermal run (thermostat temp0) to temperature temp.

    w = exp(beta0 * bias) * exp(-(beta - beta0) * U), energies in kJ/mol.
    """
    kb = 0.0083144626
    b0, b = 1 / (kb * temp0), 1 / (kb * temp)
    lw = b0 * np.asarray(bias) - (b - b0) * np.asarray(ene)
    return lw - lw.max()


def weighted_fraction(mask, logw) -> float:
    w = np.exp(logw)
    return float((w * mask).sum() / w.sum())


# ------------------------------------------------------------------ QC (review Section 2)

def log_qc(rep_dir: str) -> dict:
    """LINCS warnings from all md*.log parts; performance of the last part."""
    warn, perf = 0, None
    for log in sorted(glob.glob(os.path.join(rep_dir, "md*.log"))):
        text = open(log, errors="ignore").read()
        warn += len(re.findall(r"LINCS WARNING", text))
        m = re.search(r"Performance:\s+(\S+)", text)
        perf = float(m.group(1)) if m else perf
    return {"lincs_warnings": warn, "ns_per_day": perf}


def read_xvg(path: str) -> np.ndarray:
    return np.loadtxt([l for l in open(path) if not l.startswith(("#", "@"))], ndmin=2)


def parts(rep_dir: str, ext: str) -> list:
    """Trajectory/energy parts in time order. Empty parts (a chunk shorter than the output interval) are skipped."""
    return sorted(os.path.basename(f) for f in glob.glob(os.path.join(rep_dir, f"md*.{ext}"))
                  if not os.path.basename(f).startswith("md_whole") and os.path.getsize(f) > 0)


def read_xvg_named(path: str) -> dict:
    """{legend: column} from a gmx energy .xvg; the first column is 'time'."""
    names = ["time"] + re.findall(r'^@ s\d+ legend "(.*)"', open(path).read(), re.M)
    return dict(zip(names, read_xvg(path).T))


def energy_qc(rep_dir: str, gmx: str, env=None) -> dict:
    """Conserved-energy drift, kinetic-energy width test, and temperature, pressure, density statistics.

    Drift: one line fitted to Conserved-En. over all parts (per-part log values cover a few ps and are noise).
    KE width: in the canonical ensemble sd(T)/mean(T) = sqrt(2/ndof). A Berendsen-like thermostat gives a
    narrower distribution (ratio < 1). The test is a quick stand-in for physical_validation, whose GROMACS
    parser does not support the triclinic (dodecahedron) box.
    """
    terms = ("Kinetic-En.", "Conserved-En.", "Temperature", "Pressure", "Density")
    rows = []
    for edr in parts(rep_dir, "edr"):
        sh([gmx, "energy", "-f", edr, "-o", "qc_energy.xvg"], rep_dir, env, stdin="\n".join(terms) + "\n\n")
        rows.append(read_xvg_named(os.path.join(rep_dir, "qc_energy.xvg")))
    d = {k: np.concatenate([r[k] for r in rows]) for k in rows[0]}
    keep = np.unique(d["time"], return_index=True)[1]          # parts share their boundary frame
    d = {k: v[keep] for k, v in d.items()}
    natoms = int(re.search(r"natoms\s*=\s*(\d+)", sh([gmx, "dump", "-s", "md.tpr"], rep_dir, env)).group(1))
    t = d["time"]
    fit = lambda y: float(np.polyfit(t, y, 1)[0]) if len(t) > 2 else float("nan")   # noqa: E731
    ndof = float(np.mean(2 * d["Kinetic En."] / (0.0083144626 * d["Temperature"])))
    out = {"span_ns": float(t[-1] - t[0]) / 1000, "natoms": natoms, "n_frames": int(len(t)),
           "drift_kJ_mol_ps_atom": fit(d["Conserved En."]) / natoms,
           "ke_width_ratio": float(d["Temperature"].std() / d["Temperature"].mean() / np.sqrt(2 / ndof)),
           "ndof": ndof}
    for key, name in (("Temperature", "temperature_K"), ("Pressure", "pressure_bar"), ("Density", "density_kg_m3")):
        out[name] = {"mean": float(d[key].mean()), "sd": float(d[key].std()), "slope_per_ns": fit(d[key]) * 1000}
    return out


def min_image_distance(rep_dir: str, gmx: str, env=None) -> float:
    """Smallest distance between the protein and its own periodic image (nm), over all parts."""
    best = np.inf
    for xtc in parts(rep_dir, "xtc"):
        sh([gmx, "mindist", "-s", "md.tpr", "-f", xtc, "-pi", "-od", "qc_mindist.xvg"], rep_dir, env, stdin="Protein\n")
        best = min(best, float(read_xvg(os.path.join(rep_dir, "qc_mindist.xvg"))[:, 1].min()))
    return best


def qc_report(rep_dir: str, gmx: str, env=None, cutoff_nm: float = 1.0) -> dict:
    """Section 2.1 checks that need no extra runs. 'flags' lists anything that fails."""
    r = {**log_qc(rep_dir), "energy": energy_qc(rep_dir, gmx, env)}
    r["min_image_nm"] = min_image_distance(rep_dir, gmx, env)
    flags = []
    if r["lincs_warnings"]:
        flags.append(f"{r['lincs_warnings']} LINCS warnings")
    e = r["energy"]
    if e["span_ns"] >= 1.0 and abs(e["drift_kJ_mol_ps_atom"]) > 1e-3:     # below 1 ns the fit is noise
        flags.append(f"conserved-energy drift {e['drift_kJ_mol_ps_atom']:.1e} kJ/mol/ps/atom")
    if e["n_frames"] >= 200 and not 0.85 < e["ke_width_ratio"] < 1.15:  # 15 %: rough band for correlated frames
        flags.append(f"kinetic-energy width is {e['ke_width_ratio']:.2f} x the canonical value")
    if r["min_image_nm"] < cutoff_nm:
        flags.append(f"protein sees its own image at {r['min_image_nm']:.2f} nm")
    r["flags"] = flags
    with open(os.path.join(rep_dir, "qc.json"), "w") as f:
        json.dump(r, f, indent=2)
    return r


# ------------------------------------------------------------------ analysis (review Sections 3-4)

def whole_protein_xtc(rep_dir: str, gmx: str, env=None) -> str:
    """Concatenate parts, make the protein whole and centred, keep protein only."""
    xtcs = parts(rep_dir, "xtc")
    sh([gmx, "trjcat", "-f"] + xtcs + ["-o", "cat.xtc"], rep_dir, env)   # sorts by time, drops boundary duplicates
    sh([gmx, "trjconv", "-s", "md.tpr", "-f", "cat.xtc", "-o", "md_whole.xtc", "-pbc", "mol", "-center"],
       rep_dir, env, stdin="Protein\nProtein\n")
    sh([gmx, "trjconv", "-s", "md.tpr", "-f", "cat.xtc", "-o", "protein.pdb", "-dump", "0", "-pbc", "mol"],
       rep_dir, env, stdin="Protein\n")
    os.remove(os.path.join(rep_dir, "cat.xtc"))
    return os.path.join(rep_dir, "md_whole.xtc")


def best_hummer_q(traj, native, beta: float = 50.0, lam: float = 1.8, cutoff: float = 0.45):
    """Fraction of native contacts, Best, Hummer & Eaton PNAS 2013 (nm units: beta = 5 /A = 50 /nm)."""
    import mdtraj as md
    heavy = native.topology.select("element != H")
    pairs = np.array([(i, j) for k, i in enumerate(heavy) for j in heavy[k + 1:]
                      if abs(native.topology.atom(i).residue.index - native.topology.atom(j).residue.index) > 3])
    r0 = md.compute_distances(native[0], pairs)[0]
    native_pairs, r0 = pairs[r0 < cutoff], r0[r0 < cutoff]
    assert len(native_pairs), "native structure has no contacts: pass the folded structure as native"
    r = md.compute_distances(traj, native_pairs)
    return (1.0 / (1 + np.exp(beta * (r - lam * r0)))).mean(axis=1)


def count_transitions(q, lo: float = 0.3, hi: float = 0.8) -> dict:
    """Two-threshold transition count (review Section 3.3). Frames between lo and hi keep the last state."""
    state, last, folds, unfolds = np.full(len(q), -1), -1, 0, 0
    for i, x in enumerate(q):
        new = 1 if x > hi else 0 if x < lo else last
        if last == 0 and new == 1:
            folds += 1
        if last == 1 and new == 0:
            unfolds += 1
        state[i] = last = new
    assigned = state >= 0
    return {"folds": folds, "unfolds": unfolds, "state": state,
            "folded_fraction": float((state[assigned] == 1).mean()) if assigned.any() else float("nan")}


def stat_inefficiency(x) -> dict:
    """Equilibration cut, statistical inefficiency g and effective samples (pymbar)."""
    from pymbar import timeseries
    t0, g, neff = timeseries.detect_equilibration(np.asarray(x, dtype=float))
    return {"t0_frame": int(t0), "g": float(g), "n_eff": float(neff)}


def replica_timeseries(rep_dir: str, native_pdb: str, gmx: str, env=None, stride: int = 1):
    """Per-frame Q, Calpha RMSD to native (nm), Rg (nm) and strand fraction (DSSP 'E'). Writes timeseries.csv."""
    import mdtraj as md
    import pandas as pd
    xtc = whole_protein_xtc(rep_dir, gmx, env)
    traj = md.load(xtc, top=os.path.join(rep_dir, "protein.pdb"), stride=stride)
    native = md.load(native_pdb)
    names = lambda t: [(a.residue.index, a.name) for a in t.topology.atoms]    # noqa: E731
    assert names(native) == names(traj), "native.pdb and the trajectory have different atoms or atom order"
    ca = traj.topology.select("name CA")
    df = pd.DataFrame({
        "time_ns": traj.time / 1000.0,
        "Q": best_hummer_q(traj, native),
        "rmsd_ca_nm": md.rmsd(traj, native, atom_indices=ca),
        "rg_nm": md.compute_rg(traj),
        "strand_fraction": (md.compute_dssp(traj, simplified=True) == "E").mean(axis=1),
    })
    df.to_csv(os.path.join(rep_dir, "timeseries.csv"), index=False)
    return df


def campaign_summary(tables: dict, lo: float = 0.3, hi: float = 0.8, skip_ns: float = 0.0) -> dict:
    """Folded fraction per replica and across replicas (mean +- standard error), transition counts, g of Q."""
    per = {}
    for name, df in tables.items():
        q = df.loc[df.time_ns >= skip_ns, "Q"].to_numpy()
        t = count_transitions(q, lo, hi)
        per[name] = {"ns": float(df.time_ns.iloc[-1] - df.time_ns.iloc[0]), "folds": t["folds"],
                     "unfolds": t["unfolds"], "folded_fraction": t["folded_fraction"], **stat_inefficiency(q)}
    ff = np.array([v["folded_fraction"] for v in per.values()])
    ff = ff[np.isfinite(ff)]
    return {"replicas": per, "transitions_total": int(sum(v["folds"] + v["unfolds"] for v in per.values())),
            "folded_fraction_mean": float(ff.mean()) if len(ff) else float("nan"),
            "folded_fraction_se": float(ff.std(ddof=1) / np.sqrt(len(ff))) if len(ff) > 1 else float("nan")}


def tica_msm(trajs: list, lag: int, n_clusters: int = 100, dim: int = 2, its_lags=(1, 2, 5, 10, 20, 50), seed=0):
    """Calpha-distance features -> TICA -> k-means -> MSM (deeptime). Returns the models and implied timescales.

    Lags are in frames. Read the flat region of the implied timescales before trusting the MSM (Section 3.4).
    """
    import mdtraj as md
    from deeptime.clustering import KMeans
    from deeptime.decomposition import TICA, VAMP
    from deeptime.markov.msm import MaximumLikelihoodMSM
    from deeptime.markov import TransitionCountEstimator
    ca = trajs[0].topology.select("name CA")
    pairs = [(i, j) for k, i in enumerate(ca) for j in ca[k + 2:]]
    feats = [md.compute_distances(t, pairs).astype(np.float64) for t in trajs]
    vamp2 = VAMP(lagtime=lag, dim=dim).fit(feats).fetch_model().score(r=2)
    tica = TICA(lagtime=lag, dim=dim).fit(feats).fetch_model()
    proj = [tica.transform(f) for f in feats]
    km = KMeans(n_clusters=n_clusters, fixed_seed=seed, max_iter=200).fit(np.concatenate(proj)).fetch_model()
    dtrajs = [km.transform(p) for p in proj]
    its = {}
    for tau in its_lags:
        try:
            counts = TransitionCountEstimator(tau, "sliding").fit(dtrajs).fetch_model().submodel_largest()
            msm = MaximumLikelihoodMSM().fit(counts).fetch_model()
            its[tau] = msm.timescales(k=min(3, msm.n_states - 1)).tolist()
        except ValueError:          # too little data at this lag (early in a campaign): no estimate
            its[tau] = []
    counts = TransitionCountEstimator(lag, "sliding").fit(dtrajs).fetch_model().submodel_largest()
    msm = MaximumLikelihoodMSM().fit(counts).fetch_model()
    return {"vamp2": float(vamp2), "tica": tica, "proj": proj, "dtrajs": dtrajs, "msm": msm,
            "implied_timescales_frames": its}


def demo():
    """Self-check of the pure-Python parts (no GROMACS needed)."""
    q = np.array([0.1, 0.5, 0.9, 0.85, 0.5, 0.2, 0.5, 0.9])
    t = count_transitions(q)
    assert (t["folds"], t["unfolds"]) == (2, 1), t
    assert list(t["state"]) == [0, 0, 1, 1, 1, 0, 0, 1]
    assert abs(t["folded_fraction"] - 0.5) < 1e-9
    noisy = np.array([0.5, 0.79, 0.81, 0.79, 0.81])        # never leaves the middle band: no state, no transition
    assert count_transitions(noisy)["folds"] == 0
    cmd = mdrun_command("gmx_mpi", ["r0", "r1"], 1.0, 2, True)
    assert cmd[:3] == ["mpirun", "-np", "2"] and cmd[-3:] == ["-multidir", "r0", "r1"]
    assert "-multidir" not in mdrun_command("gmx", ["r0"], 1.0, 2, False)
    pdb = ("ATOM      1  N   TYR A   1      0.000   0.000   0.000  1.00  0.00           N\n"
           "ATOM      2  CA ATYR A   1      1.000   0.000   0.000  0.50  0.00           C\n"
           "ATOM      3  CA BTYR A   1      1.100   0.000   0.000  0.50  0.00           C\n")
    try:
        clean_crystal(pdb)
    except AssertionError as e:
        assert "expected YYDPETGTWY" in str(e)
    # multithermal reweighting: at temp == temp0 the weights reduce to exp(beta0 * bias)
    ene, bias = np.array([-1000.0, -1010.0]), np.array([0.0, 2.0])
    lw = multithermal_logweights(ene, bias, 340, 340)
    assert np.allclose(lw - lw[0], [0, 2 / (0.0083144626 * 340)])
    hot = multithermal_logweights(ene, np.zeros(2), 340, 400)          # a hotter target favours the higher energy
    assert hot[0] > hot[1]
    assert abs(weighted_fraction(np.array([1, 0]), np.log([0.25, 0.75])) - 0.25) < 1e-12
    assert "-update" not in mdrun_command("gmx", ["r0"], 1.0, 2, True, plumed="plumed.dat")
    print("cln025.demo: all checks passed")


if __name__ == "__main__":
    demo()
