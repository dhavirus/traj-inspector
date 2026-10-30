"""Engine detection and normalized Setup / Energies per replica (GROMACS, Amber, OpenMM)."""
from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import MDAnalysis as mda
import numpy as np

KCAL_TO_KJ = 4.184
HYDROGEN_MASS_AMU = 1.008
HEAVY_MASS_AMU = 5.0   # shortcut: OpenMM constraint class from masses; a particle above this is not H even with HMR (<= 4 amu)
ENERGY_FIELDS = ("time_ps", "temperature", "kinetic", "potential", "conserved", "density_kg_m3", "volume")


class ReaderError(RuntimeError):
    """The folder cannot be read; the message says which file or format is the problem."""


@dataclass(frozen=True)
class Setup:
    engine: str
    thermostat: str            # berendsen, v-rescale, nose-hoover, langevin, andersen, none, ...
    barostat: str              # berendsen, c-rescale, parrinello-rahman, monte-carlo, none, ...
    dt_ps: float
    constraints: str           # none, h-bonds, all-bonds, ...
    cutoff_nm: float
    hmr_factor: float | None   # mean solute-H mass / 1.008; None when the selection has no hydrogens
    warnings_count: int | None  # None = the engine does not log them


@dataclass(frozen=True)
class Energies:
    """Arrays share one length; absent quantities are NaN. Energies in kJ/mol, volume in nm^3."""
    time_ps: np.ndarray
    temperature: np.ndarray
    kinetic: np.ndarray
    potential: np.ndarray
    conserved: np.ndarray | None   # None when the integrator has no conserved quantity (Langevin, Berendsen, ...)
    density_kg_m3: np.ndarray
    volume: np.ndarray


@dataclass(frozen=True)
class Replica:
    name: str
    engine: str
    topology: Path
    trajectories: tuple[Path, ...]
    setup: Setup
    energies: Energies


def _natural(path: Path) -> list:
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", path.name)]


def _files(folder: Path, *patterns: str) -> list[Path]:
    """Non-empty files matching the patterns, natural order (md.part0002 before md.part0010)."""
    found = {p for pat in patterns for p in folder.glob(pat) if p.is_file() and p.stat().st_size > 0}
    return sorted(found, key=_natural)


def detect_engine(folder: Path) -> str:
    if _files(folder, "*.tpr"):
        return "gromacs"
    if _files(folder, "system.xml", "integrator.xml"):
        return "openmm"
    if _files(folder, "mdout*", "*.mdout", "*.out"):
        return "amber"
    raise ReaderError(f"{folder}: no .tpr (GROMACS), system.xml/integrator.xml (OpenMM) or mdout (Amber) found")


def _one(candidates: list[Path], what: str, folder: Path, prefer: str = "md") -> Path:
    if not candidates:
        raise ReaderError(f"{folder}: no {what} file")
    if len(candidates) > 1:
        preferred = [c for c in candidates if c.name.startswith(prefer)]
        if len(preferred) != 1:
            raise ReaderError(f"{folder}: several {what} files {[c.name for c in candidates]}; keep one")
        return preferred[0]
    return candidates[0]


def _kv_pairs(text: str) -> dict[str, str]:
    """key = value pairs, lower-cased; '=' with comma/whitespace separators (Amber namelist, GROMACS mdp)."""
    return {k.lower(): v for k, v in re.findall(r"(\w+)\s*=\s*([^,\s]+)", text)}


# ------------------------------------------------------------------ hydrogen mass

def hmr_factor(universe: mda.Universe, selection: str) -> float | None:
    """Mean mass of solute hydrogens / 1.008. Water H stays 1.008 (rigid), so only the selection is looked at."""
    hydrogens = universe.select_atoms(f"({selection}) and name H*")
    if len(hydrogens) == 0:
        return None
    return float(hydrogens.masses.mean() / HYDROGEN_MASS_AMU)


# ------------------------------------------------------------------ GROMACS

def _read_mdp(path: Path) -> dict[str, str]:
    pairs = {}
    for line in path.read_text().splitlines():
        line = line.split(";")[0]
        if "=" in line:
            key, value = line.split("=", 1)
            pairs[key.strip().lower().replace("_", "-")] = value.strip().lower()
    return pairs


_GMX_THERMOSTAT = {"no": "none", "v-rescale": "v-rescale", "berendsen": "berendsen", "nose-hoover": "nose-hoover",
                   "andersen": "andersen", "andersen-massive": "andersen"}
_GMX_BAROSTAT = {"no": "none", "c-rescale": "c-rescale", "berendsen": "berendsen",
                 "parrinello-rahman": "parrinello-rahman", "mttk": "mttk"}


def _gromacs_setup(folder: Path, universe: mda.Universe, selection: str) -> Setup:
    mdp_path = _one(_files(folder, "*.mdp"), "mdp", folder, prefer="mdout")
    mdp = _read_mdp(mdp_path)
    integrator = mdp.get("integrator", "md")
    thermostat = "langevin" if integrator in ("sd", "bd") else _GMX_THERMOSTAT.get(mdp.get("tcoupl", "no"))
    barostat = _GMX_BAROSTAT.get(mdp.get("pcoupl", "no"))
    if thermostat is None or barostat is None:
        raise ReaderError(f"{mdp_path}: unknown tcoupl/pcoupl value ({mdp.get('tcoupl')}, {mdp.get('pcoupl')})")
    logs = _files(folder, "*.log")
    warnings = sum(len(re.findall(r"LINCS WARNING", p.read_text(errors="ignore"))) for p in logs)
    return Setup("gromacs", thermostat, barostat, float(mdp["dt"]), mdp.get("constraints", "none"),
                 max(float(mdp.get("rvdw", 1.0)), float(mdp.get("rcoulomb", 1.0))),
                 hmr_factor(universe, selection), warnings)


def _gromacs_energies(folder: Path) -> Energies:
    import pyedr
    parts = [pyedr.edr_to_dict(str(p)) for p in _files(folder, "*.edr")]
    if not parts:
        raise ReaderError(f"{folder}: no non-empty .edr file")
    keys = ("Time", "Temperature", "Kinetic En.", "Potential", "Conserved En.", "Density", "Volume")
    table = {k: np.concatenate([p.get(k, np.full(len(p["Time"]), np.nan)) for p in parts]) for k in keys}
    time, keep = np.unique(table["Time"], return_index=True)   # parts share their boundary frame
    conserved = table["Conserved En."][keep]
    return Energies(time, table["Temperature"][keep], table["Kinetic En."][keep], table["Potential"][keep],
                    None if np.isnan(conserved).all() else conserved, table["Density"][keep], table["Volume"][keep])


# ------------------------------------------------------------------ Amber

_AMBER_THERMOSTAT = {0: "none", 1: "berendsen", 2: "andersen", 3: "langevin"}
_AMBER_CONSTRAINTS = {1: "none", 2: "h-bonds", 3: "all-bonds"}
_AMBER_BLOCK = re.compile(r"NSTEP\s*=\s*\d+\s+TIME\(PS\)\s*=\s*([-\d.]+)\s+TEMP\(K\)\s*=\s*([-\d.]+)")


def _amber_outputs(folder: Path) -> list[tuple[Path, str]]:
    """mdout files of production runs (imin=0), time order by name; minimization outputs are skipped."""
    outs = []
    for path in _files(folder, "mdout*", "*.mdout", "*.out"):
        text = path.read_text(errors="ignore")
        echo = text.split("Here is the input file:")[-1].split("Nature and format")[0] if "Here is the input file:" in text else ""
        if echo and _kv_pairs(echo).get("imin") == "0":
            outs.append((path, text))
    if not outs:
        raise ReaderError(f"{folder}: no Amber mdout of an MD run (imin=0) with the echoed input file")
    return outs


def _amber_setup(folder: Path, universe: mda.Universe, selection: str) -> Setup:
    outs = _amber_outputs(folder)
    text = outs[0][1]
    echo = text.split("Here is the input file:")[1].split("Nature and format")[0]
    control = text.split("ATOMIC COORDINATES")[0].split("2.  CONTROL")[-1]
    pairs = {**_kv_pairs(control), **_kv_pairs(echo)}   # control data supplies defaults (cut, dt); the echo wins
    ntt, ntp = int(pairs.get("ntt", 0)), int(pairs.get("ntp", 0))
    if ntt not in _AMBER_THERMOSTAT:
        raise ReaderError(f"{outs[0][0]}: unsupported ntt={ntt}")
    barostat = "none" if ntp == 0 else ("monte-carlo" if "Monte-Carlo Barostat" in text else "berendsen")
    warnings = sum(t.count("vlimit exceeded") for _, t in outs)
    return Setup("amber", _AMBER_THERMOSTAT[ntt], barostat, float(pairs["dt"]),
                 _AMBER_CONSTRAINTS[int(pairs.get("ntc", 1))], float(pairs["cut"]) / 10.0,
                 hmr_factor(universe, selection), warnings)


def _amber_value(block: str, name: str) -> float:
    m = re.search(rf"\b{name}\s*=\s*([-\d.]+)", block)
    return float(m.group(1)) if m else np.nan


def _amber_energies(folder: Path, setup_thermostat: str, setup_barostat: str) -> Energies:
    rows = []
    for _, text in _amber_outputs(folder):
        body = text.split("A V E R A G E S")[0]
        for block in re.split(r"(?=\sNSTEP\s*=)", body)[1:]:
            m = _AMBER_BLOCK.search(block)
            if not m:
                continue
            rows.append([float(m.group(1)), float(m.group(2)), _amber_value(block, "EKtot") * KCAL_TO_KJ, _amber_value(block, "EPtot") * KCAL_TO_KJ,
                         _amber_value(block, "Etot") * KCAL_TO_KJ, _amber_value(block, "Density") * 1000.0, _amber_value(block, "VOLUME") * 1e-3])
    if not rows:
        raise ReaderError(f"{folder}: no 'NSTEP = ... TIME(PS) = ... TEMP(K) = ...' blocks in the mdout")
    t, temp, ke, pe, etot, dens, vol = np.array(rows).T
    time, keep = np.unique(t, return_index=True)
    nve = setup_thermostat == "none" and setup_barostat == "none"   # Etot is conserved only without any coupling
    return Energies(time, temp[keep], ke[keep], pe[keep], etot[keep] if nve else None, dens[keep], vol[keep])


# ------------------------------------------------------------------ OpenMM

_OMM_THERMOSTAT = {"LangevinMiddleIntegrator": "langevin", "LangevinIntegrator": "langevin",
                   "VerletIntegrator": "none", "NoseHooverIntegrator": "nose-hoover", "BrownianIntegrator": "brownian"}


def _openmm_masses(system_root: ET.Element) -> np.ndarray:
    return np.array([float(p.get("mass")) for p in system_root.find("Particles")])


def _openmm_constraints(root: ET.Element, masses: np.ndarray) -> str:
    cons = root.find("Constraints")
    if cons is None or len(cons) == 0:
        return "none"
    heavy_pairs = any(masses[int(c.get("p1"))] > HEAVY_MASS_AMU and masses[int(c.get("p2"))] > HEAVY_MASS_AMU for c in cons)
    return "all-bonds" if heavy_pairs else "h-bonds"


def _openmm_setup(folder: Path, universe: mda.Universe, selection: str) -> tuple[Setup, np.ndarray]:
    system_path = _one(_files(folder, "system.xml"), "system.xml", folder)
    integrator_path = _one(_files(folder, "integrator.xml"), "integrator.xml", folder)
    system, integrator = ET.parse(system_path).getroot(), ET.parse(integrator_path).getroot()
    forces = {f.get("type"): f for f in system.find("Forces")}
    masses = _openmm_masses(system)
    thermostat = "andersen" if "AndersenThermostat" in forces else _OMM_THERMOSTAT.get(integrator.get("type"))
    if thermostat is None:
        raise ReaderError(f"{integrator_path}: unsupported integrator type '{integrator.get('type')}'")
    has_baro = any(t.startswith("MonteCarlo") and "Barostat" in t for t in forces)
    if "NonbondedForce" not in forces:
        raise ReaderError(f"{system_path}: no NonbondedForce, cannot read the cutoff")
    universe.atoms.masses = masses   # a PDB topology has no masses; MDAnalysis would guess H = 1.008 and hide HMR
    return Setup("openmm", thermostat, "monte-carlo" if has_baro else "none", float(integrator.get("stepSize")),
                 _openmm_constraints(system, masses), float(forces["NonbondedForce"].get("cutoff")),
                 hmr_factor(universe, selection), None), masses


_OMM_COLUMNS = {"Time (ps)": "time_ps", "Temperature (K)": "temperature", "Kinetic Energy (kJ/mole)": "kinetic",
                "Potential Energy (kJ/mole)": "potential", "Density (g/mL)": "density_kg_m3", "Box Volume (nm^3)": "volume"}


def _openmm_energies(folder: Path, thermostat: str, barostat: str) -> Energies:
    path = _one(_files(folder, "*.csv"), "StateDataReporter csv", folder, prefer="state")
    lines = path.read_text().splitlines()
    header = [h.strip('"') for h in lines[0].lstrip("#").split(",")]
    if "Time (ps)" not in header:
        raise ReaderError(f"{path}: header has no 'Time (ps)' column: {header}")
    data = np.loadtxt(lines[1:], delimiter=",", ndmin=2)
    col = {name: data[:, header.index(name)] if name in header else np.full(len(data), np.nan) for name in _OMM_COLUMNS}
    out = {_OMM_COLUMNS[name]: values for name, values in col.items()}
    out["density_kg_m3"] = out["density_kg_m3"] * 1000.0
    nve = thermostat == "none" and barostat == "none" and "Total Energy (kJ/mole)" in header
    time, keep = np.unique(out["time_ps"], return_index=True)
    out = {k: v[keep] for k, v in out.items()}
    out["conserved"] = data[:, header.index("Total Energy (kJ/mole)")][keep] if nve else None
    return Energies(**out)


# ------------------------------------------------------------------ public

def _topology_and_trajectories(folder: Path, engine: str) -> tuple[Path, list[Path]]:
    if engine == "gromacs":
        return _one(_files(folder, "*.tpr"), "tpr", folder), _files(folder, "*.xtc", "*.trr")
    if engine == "amber":
        return _one(_files(folder, "*.prmtop", "*.parm7"), "prmtop/parm7", folder), _files(folder, "*.nc", "*.mdcrd", "*.dcd")
    topology = _files(folder, "*.pdb", "*.cif", "*.prmtop", "*.parm7")
    return _one(topology, "topology (pdb/cif/prmtop)", folder, prefer="topology"), _files(folder, "*.dcd", "*.xtc")


def read_replica(folder: str | Path, selection: str = "protein") -> Replica:
    """Detect the engine of `folder`, read Setup and Energies. Raises ReaderError naming the missing/odd file."""
    folder = Path(folder)
    engine = detect_engine(folder)
    topology, trajectories = _topology_and_trajectories(folder, engine)
    if not trajectories:
        raise ReaderError(f"{folder}: no non-empty trajectory file for {engine}")
    universe = mda.Universe(str(topology))
    if engine == "gromacs":
        setup, energies = _gromacs_setup(folder, universe, selection), _gromacs_energies(folder)
    elif engine == "amber":
        setup = _amber_setup(folder, universe, selection)
        energies = _amber_energies(folder, setup.thermostat, setup.barostat)
    else:
        setup, _ = _openmm_setup(folder, universe, selection)
        energies = _openmm_energies(folder, setup.thermostat, setup.barostat)
    return Replica(folder.name, engine, topology, tuple(trajectories), setup, energies)


def open_universe(replica: Replica) -> tuple[mda.Universe, np.ndarray]:
    """Universe over all trajectory parts plus the indices of frames to use.

    Parts are chained plainly (continuous=True fails on parts with a single frame); frames whose time is not
    after the previous kept frame (boundary duplicates between -noappend parts) are dropped.
    """
    universe = mda.Universe(str(replica.topology), [str(p) for p in replica.trajectories])
    if replica.engine == "openmm":
        system = ET.parse(_one(_files(replica.topology.parent, "system.xml"), "system.xml", replica.topology.parent)).getroot()
        universe.atoms.masses = _openmm_masses(system)
    times = np.array([ts.time for ts in universe.trajectory])
    running_max = np.maximum.accumulate(np.concatenate([[-np.inf], times[:-1]]))
    return universe, np.flatnonzero(times > running_max)
