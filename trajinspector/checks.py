"""Gates 1-4: metrics per replica and for the campaign -> a results table.

Results hold the *metric*, not the verdict. `evaluate_results` applies the current rules, so editing a threshold
re-renders the report without recomputing anything.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable

import MDAnalysis as mda
import numpy as np
from MDAnalysis.transformations import unwrap

from . import analysis as A
from . import rules as R
from .readers import Replica, open_universe

CAMPAIGN = "campaign"
KB_KJ_PER_MOL_K = 0.0083144626
HMR_VARIANT = "hmr"
BOND_CONSTRAINTS = ("h-bonds", "h-angles", "all-bonds", "all-angles")
ANGSTROM_PER_NM = 10.0
HMR_DETECT_FACTOR = 1.5   # mean solute-H mass above 1.5 x 1.008 means repartitioning (typical HMR: 3-4 amu); rule of thumb


@dataclass(frozen=True)
class Result:
    rule_id: str
    replica: str          # replica name or CAMPAIGN
    metric: object        # number, str, bool or None
    detail: str = ""      # human-readable value with units
    reason: str = ""      # why the metric is missing (status N/A) or a caveat
    variant: str | None = None   # selects an alternative threshold, e.g. value_hmr


@dataclass
class Series:
    """Per-replica arrays kept for plots and the viewer. Time in ps, lengths in nm."""
    time_ps: np.ndarray
    q: np.ndarray
    rg_nm: np.ndarray
    ca_features: np.ndarray
    dt_ps: float
    n_atoms: int


@dataclass
class Campaign:
    results: list[Result] = field(default_factory=list)
    series: dict[str, Series] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)   # caveats shown at the top of the report
    its_ns: dict[int, float | None] = field(default_factory=dict)
    replicas: list[Replica] = field(default_factory=list)
    msm_lag_ns: float | None = None

    def add(self, *results: Result) -> None:
        self.results.extend(results)


@dataclass(frozen=True)
class Row:
    result: Result
    rule: R.Rule
    status: str


def evaluate_results(ruleset: R.RuleSet, results: list[Result]) -> list[Row]:
    rows = []
    for res in results:
        rule = ruleset.get(res.rule_id)
        if res.variant == HMR_VARIANT:
            rule = replace(rule, value=rule.extra["value_hmr"])
        status = R.NA if res.metric is None and rule.kind not in ("manual", "info") else R.evaluate(rule, res.metric)
        rows.append(Row(res, rule, status))
    return rows


# ------------------------------------------------------------------ trajectory preparation

@dataclass
class Prepared:
    replica: Replica
    universe: mda.Universe
    group: mda.AtomGroup
    frames: np.ndarray


def prepare(replica: Replica, selection: str, max_frames: int) -> Prepared:
    """Open the replica, subsample to at most max_frames and unwrap the selection (guessing bonds if absent)."""
    universe, frames = open_universe(replica)
    group = universe.select_atoms(selection)
    if len(group) == 0:
        raise ValueError(f"{replica.name}: selection '{selection}' matches no atoms")
    if not hasattr(universe, "bonds"):
        group.guess_bonds()   # a PDB topology has no bonds; unwrap needs them
    universe.trajectory.add_transformations(unwrap(group))
    stride = math.ceil(len(frames) / max_frames)
    return Prepared(replica, universe, group, frames[::stride])


@dataclass(frozen=True)
class Reference:
    keys: list
    xyz_angstrom: np.ndarray
    source: str


def _heavy(group: mda.AtomGroup) -> mda.AtomGroup:
    return group.select_atoms("not name H*")


def reference_from_file(path: Path, selection: str) -> Reference:
    heavy = _heavy(mda.Universe(str(path)).select_atoms(selection))
    if len(heavy) == 0:
        raise ValueError(f"{path}: selection '{selection}' matches no heavy atoms")
    return Reference(A.atom_keys(heavy.resindices, heavy.names), heavy.positions.copy(), str(path))


def reference_from_frame0(prep: Prepared) -> Reference:
    prep.universe.trajectory[prep.frames[0]]
    heavy = _heavy(prep.group)
    return Reference(A.atom_keys(heavy.resindices, heavy.names), heavy.positions.copy(),
                     f"frame 0 of replica {prep.replica.name}")


# ------------------------------------------------------------------ per replica

def _trajectory_pass(prep: Prepared, contacts: A.NativeContacts, params: dict, ca_pairs: np.ndarray) -> dict:
    """One sweep over the analysed frames: Q, Rg, C-alpha features, frame-to-frame RMSD jump, min image distance."""
    from MDAnalysis.analysis.rms import rmsd
    ca = prep.group.select_atoms("name CA")
    q, rg, feats, jumps, images, times = [], [], [], [], [], []
    previous = None
    for n, frame in enumerate(prep.frames):
        ts = prep.universe.trajectory[frame]
        times.append(ts.time)
        q.append(A.best_hummer_q(prep.group.positions, contacts))
        rg.append(prep.group.radius_of_gyration() / ANGSTROM_PER_NM)
        pos = ca.positions
        feats.append(np.linalg.norm(pos[ca_pairs[:, 0]] - pos[ca_pairs[:, 1]], axis=1))
        if previous is not None:
            jumps.append(rmsd(pos, previous, superposition=True) / ANGSTROM_PER_NM)
        previous = pos.copy()
        if n % params["image_stride"] == 0 or n == 0:
            images.append(A.min_image_distance_angstrom(prep.group.positions, ts.dimensions) / ANGSTROM_PER_NM)
    return {"time_ps": np.array(times), "q": np.array(q), "rg": np.array(rg), "feats": np.array(feats),
            "jumps": np.array(jumps), "images": np.array(images)}


def _setup_results(rep: Replica) -> list[Result]:
    s, name = rep.setup, rep.name
    dt_fs = s.dt_ps * 1000.0
    is_hmr = None if s.hmr_factor is None else s.hmr_factor > HMR_DETECT_FACTOR
    constraints_ok = s.constraints in BOND_CONSTRAINTS or dt_fs <= 1.0
    return [
        Result("thermostat", name, s.thermostat, s.thermostat),
        Result("barostat", name, s.barostat, s.barostat),
        Result("forcefield_water", name, None, "", "water model is not recorded in the engine output"),
        Result("dt_max_fs", name, dt_fs, f"{dt_fs:g} fs", variant=HMR_VARIANT if is_hmr else None),
        Result("hmr_rates", name, is_hmr, "n/a" if is_hmr is None else f"mean solute-H mass x{s.hmr_factor:.2f}",
               "no hydrogens in the selection" if is_hmr is None else ""),
        Result("constraints", name, constraints_ok, f"constraints={s.constraints}, dt={dt_fs:g} fs"),
        Result("cutoff_min_nm", name, s.cutoff_nm, f"{s.cutoff_nm:.2f} nm"),
    ]



def _energy_results(rep: Replica, n_atoms: int, params: dict) -> list[Result]:
    e, name, out = rep.energies, rep.name, []
    span_ns = float(e.time_ps[-1] - e.time_ps[0]) / 1000.0
    # drift
    if e.conserved is None:
        out.append(Result("drift", name, None, reason=f"no conserved energy with {rep.setup.thermostat} thermostat"))
    elif span_ns < params["drift_min_ns"]:
        out.append(Result("drift", name, None, reason=f"run is {span_ns:.3g} ns, fit needs >= {params['drift_min_ns']} ns"))
    else:
        d = abs(A.slope(e.time_ps, e.conserved)) / n_atoms
        out.append(Result("drift", name, d, f"{d:.2e} kJ/mol/ps/atom over {span_ns:.2f} ns"))
    # kinetic-energy width
    ok = np.isfinite(e.temperature) & np.isfinite(e.kinetic)
    if ok.sum() < params["ke_min_frames"]:
        out.append(Result("ke_width", name, None, reason=f"{int(ok.sum())} energy frames, test needs >= {params['ke_min_frames']}"))
    else:
        ndof = float(np.mean(2 * e.kinetic[ok] / (KB_KJ_PER_MOL_K * e.temperature[ok])))
        ratio = float(e.temperature[ok].std() / e.temperature[ok].mean() / np.sqrt(2 / ndof))
        out.append(Result("ke_width", name, ratio, f"{ratio:.3f} (ndof {ndof:.0f})"))
    # warnings
    w = rep.setup.warnings_count
    out.append(Result("warnings", name, w, "" if w is None else f"{w} warnings", "" if w is not None else f"{rep.engine} does not log them"))
    # density
    start = int(len(e.time_ps) * params["density_equil_fraction"])
    s = A.slope(e.time_ps[start:], e.density_kg_m3[start:]) * 1000.0
    if math.isnan(s):
        out.append(Result("density_drift", name, None, reason="no density series or fewer than 3 points after equilibration"))
    else:
        out.append(Result("density_drift", name, abs(s), f"{abs(s):.2f} kg/m3/ns (last {1 - params['density_equil_fraction']:.0%} of run)"))
    return out


def analyze_replica(prep: Prepared, reference: Reference, params: dict, campaign: Campaign) -> None:
    rep = prep.replica
    keys = A.atom_keys(prep.group.resindices, prep.group.names)
    contacts, dropped = A.native_contacts(reference.keys, reference.xyz_angstrom, keys)
    if dropped:
        campaign.notes.append(f"{rep.name}: {dropped} reference heavy atoms have no match in the trajectory (matched by residue index + atom name); their contacts are dropped.")
    ca_pairs, subsampled = A.ca_feature_pairs(len(prep.group.select_atoms("name CA")))
    if subsampled:
        campaign.notes.append(f"shortcut: C-alpha features subsampled to {A.MAX_CA_FEATURE_ATOMS} atoms for TICA/MSM.")
    out = _trajectory_pass(prep, contacts, params, ca_pairs)
    dt_ps = float(np.median(np.diff(out["time_ps"]))) if len(out["time_ps"]) > 1 else float("nan")
    campaign.series[rep.name] = Series(out["time_ps"], out["q"], out["rg"], out["feats"], dt_ps, len(prep.group))
    cutoff = rep.setup.cutoff_nm
    box_pad = out["images"][0] / 2
    min_image = float(out["images"].min())
    jump = float(out["jumps"].max()) if len(out["jumps"]) else None
    campaign.add(
        *_setup_results(rep), *_energy_results(rep, len(prep.group), params),
        Result("box_padding_nm", rep.name, box_pad, f"{box_pad:.2f} nm (frame-0 min image distance {out['images'][0]:.2f} nm)"),
        Result("min_image", rep.name, min_image - cutoff, f"min image distance {min_image:.2f} nm vs cutoff {cutoff:.2f} nm, {len(out['images'])} frames sampled"),
        Result("rmsd_jump", rep.name, jump, "" if jump is None else f"{jump:.3f} nm (stride {int(prep.frames[1] - prep.frames[0]) if len(prep.frames) > 1 else 1})",
               "" if jump is not None else "single frame"),
    )


# ------------------------------------------------------------------ gate 3 / 4 from Q

def _q_results(name: str, s: Series, params: dict) -> tuple[list[Result], dict]:
    """Per-replica Q statistics. Returns results and a dict with n_eff, folded fraction and transition counts.

    Everything except the INFO row uses all frames: pymbar's t0 is not a safe cut for a two-state series
    (see analysis.stat_inefficiency). Drop an initial segment by choosing a later reference/start instead.
    """
    lo, hi = params["q_unfolded"], params["q_folded"]
    stats = A.stat_inefficiency(s.q)
    t0 = stats["t0_frame"]
    trans = A.count_transitions(s.q, lo, hi)
    ratio = A.block_se_ratio(s.q, params["min_blocks"])
    between = float(((s.q > lo) & (s.q < hi)).mean())
    no_g = math.isnan(stats["n_eff"])
    out = [
        Result("equilibration", name, t0, f"pymbar t0 = {s.time_ps[t0] / 1000:.4g} ns (frame {t0} of {len(s.q)}); not applied"),
        Result("n_eff_min", name, None if no_g else stats["n_eff"],
               "" if no_g else f"N_eff = {stats['n_eff']:.1f} (g = {stats['g']:.1f}, {len(s.q)} frames)",
               "Q is constant: statistical inefficiency undefined" if no_g else ""),
        Result("block_plateau", name, None if math.isnan(ratio) else ratio, "" if math.isnan(ratio) else f"{ratio:.2f}",
               f"needs >= {2 * params['min_blocks']} frames and some variance" if math.isnan(ratio) else ""),
        Result("q_bimodal", name, between, f"{between:.1%} of frames between {lo} and {hi}"),
    ]
    return out, {"n_eff": stats["n_eff"], "p": trans["folded_fraction"], "transitions": trans["folds"] + trans["unfolds"]}


def _agreement(per_rep: dict[str, dict]) -> Result:
    names = list(per_rep)
    if len(names) < 2:
        return Result("replica_agreement", CAMPAIGN, None, reason="needs at least 2 replicas")
    worst, pair = 0.0, ("", "")
    for a, b in ((names[i], names[j]) for i in range(len(names)) for j in range(i + 1, len(names))):
        pa, pb = per_rep[a]["p"], per_rep[b]["p"]
        if math.isnan(pa) or math.isnan(pb) or math.isnan(per_rep[a]['n_eff']) or math.isnan(per_rep[b]['n_eff']):
            return Result("replica_agreement", CAMPAIGN, None, reason=f"{a} or {b} has no folded fraction or N_eff (Q never left the middle zone, or is constant)")
        se2 = sum(p * (1 - p) / max(per_rep[n]["n_eff"], 1.0) for p, n in ((pa, a), (pb, b)))
        z = abs(pa - pb) / math.sqrt(se2) if se2 > 0 else (0.0 if pa == pb else math.inf)
        if z >= worst:
            worst, pair = z, (a, b)
    return Result("replica_agreement", CAMPAIGN, worst, f"max |ff_i - ff_j| / sqrt(se_i^2+se_j^2) = {worst:.2f} ({pair[0]} vs {pair[1]})")


def _msm_results(campaign: Campaign, params: dict) -> None:
    series = list(campaign.series.values())
    dt_ps = float(np.nanmedian([s.dt_ps for s in series]))
    lag = max(1, round(params["msm_lag_ns"] * 1000.0 / dt_ps))
    campaign.msm_lag_ns = lag * dt_ps / 1000.0
    try:
        summary = A.msm_summary([s.ca_features for s in series], lag, params["seed"])
    except ValueError as exc:
        reason = f"MSM not estimable at lag {lag} frames ({lag * dt_ps / 1000:.3g} ns): {exc}"
        campaign.add(*(Result(i, CAMPAIGN, None, reason=reason) for i in ("vamp2", "its_flat", "ck_test")))
        return
    its = summary["its_frames"]
    campaign.its_ns = {int(k): None if v is None else v * dt_ps / 1000.0 for k, v in its.items()}
    campaign.add(Result("vamp2", CAMPAIGN, summary["vamp2"], f"VAMP-2 = {summary['vamp2']:.3f} at lag {lag} frames"))
    a, b = its.get(lag), its.get(2 * lag)
    if a is None or b is None:
        campaign.add(Result("its_flat", CAMPAIGN, None, reason=f"no timescale at lag {lag} or {2 * lag} frames (too little data)"))
    else:
        change = abs(b - a) / a
        campaign.add(Result("its_flat", CAMPAIGN, change, f"slowest ITS {a * dt_ps / 1000:.3g} -> {b * dt_ps / 1000:.3g} ns ({change:.1%})"))
    ck = summary["ck_error"]
    campaign.add(Result("ck_test", CAMPAIGN, ck, "" if ck is None else f"max CK error {ck:.3f}",
                        "" if ck is not None else "too little data for the Chapman-Kolmogorov test"))


def run_checks(replicas: list[Replica], ruleset: R.RuleSet, selection: str = "protein", reference: Path | None = None,
               progress: Callable[[str], None] = lambda message: None) -> Campaign:
    """Run gates 1-4. `reference` is a structure file for Q; without it frame 0 of the first replica is used."""
    params = ruleset.params
    campaign = Campaign(replicas=list(replicas))
    preps = []
    for rep in replicas:
        progress(f"opening {rep.name}")
        preps.append(prepare(rep, selection, params["max_frames"]))
    if reference is not None:
        ref = reference_from_file(Path(reference), selection)
    else:
        ref = reference_from_frame0(preps[0])
        campaign.notes.append("No reference structure given: native contacts for Q come from frame 0 of the first replica (WARN).")
    per_rep = {}
    for prep in preps:
        progress(f"analysing {prep.replica.name}")
        analyze_replica(prep, ref, params, campaign)
        results, per_rep[prep.replica.name] = _q_results(prep.replica.name, campaign.series[prep.replica.name], params)
        campaign.add(*results)
    total = sum(v["transitions"] for v in per_rep.values())
    campaign.add(
        Result("replicas_min", CAMPAIGN, len(replicas), f"{len(replicas)} replicas"),
        _agreement(per_rep),
        *(Result(i, CAMPAIGN, total, f"{total} transitions (all replicas)") for i in ("transitions", "transitions_tm", "transitions_rates")),
        Result("group_temps", CAMPAIGN, None, reason="not in v1: needs velocities"),
        Result("fes", CAMPAIGN, None, "see the plot", "plot only: no barriers from projections"),
    )
    progress("building MSM")
    _msm_results(campaign, params)
    return campaign
