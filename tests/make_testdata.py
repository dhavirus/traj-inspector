"""Regenerate the CLN025 (YYDPETGTWY) test systems under tests/data/{amber,gromacs,openmm}.

Not run by pytest: the outputs are committed. Needs `gmx` (apt gromacs), AmberTools (`tleap`, `sander` on PATH)
and `openmm` (pip). Usage: python tests/make_testdata.py [outdir]

Deviations from the CLAUDE.md recipe, forced by this build host (no GitHub downloads, GROMACS 2023.3 from apt):
  * GROMACS uses amber99sb-ildn + TIP3P shipped with the apt package, not a99SB-disp.
  * Amber writes a frame every 50 steps (60 frames, not 120) to keep the committed .nc at ~4 MB.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

SEQUENCE = "NTYR TYR ASP PRO GLU THR GLY THR TRP CTYR"   # CLN025 YYDPETGTWY with termini caps
ENV = {**os.environ, "GMX_MAXBACKUP": "-1"}

LEAP = f"""source leaprc.protein.ff14SB
source leaprc.water.tip3p
m = sequence {{ {SEQUENCE} }}
savepdb m solute.pdb
solvatebox m TIP3PBOX 8.0
addIons m Na+ 0
saveamberparm m cln025.prmtop cln025.inpcrd
quit
"""
SANDER_MIN = "min\n &cntrl\n  imin=1, maxcyc=200, ncyc=100, ntb=1, cut=8.0, ntpr=50,\n /\n"
SANDER_MD = """cln025 6 ps langevin npt
 &cntrl
  imin=0, irest=0, ntx=1, nstlim=3000, dt=0.002,
  ntc=2, ntf=2, cut=8.0, ntb=2, ntp=1, barostat=2, pres0=1.0, taup=2.0,
  ntt=3, gamma_ln=2.0, temp0=300.0, tempi=300.0, ig=12345,
  ntpr=25, ntwx=50, ntwr=3000, ioutfm=1,
 /
"""
MDP = """integrator = md
dt = 0.002
nsteps = 2500
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb = 1.0
rvdw = 1.0
DispCorr = EnerPres
constraints = h-bonds
tcoupl = v-rescale
tc-grps = System
tau-t = 0.5
ref-t = 300
pcoupl = C-rescale
tau-p = 2.0
ref-p = 1.0
compressibility = 4.5e-5
gen-vel = yes
gen-temp = 300
gen-seed = 7
nstxout-compressed = 50
nstenergy = 50
nstlog = 500
nstcalcenergy = 50
"""
EM_MDP = "integrator = steep\nnsteps = 500\nemtol = 500\ncutoff-scheme = Verlet\ncoulombtype = PME\n" \
         "rcoulomb = 1.0\nrvdw = 1.0\nconstraints = h-bonds\n"


def run(cmd: list[str], cwd: Path, stdin: str | None = None) -> None:
    """Run a command, raise with the output tail on failure."""
    r = subprocess.run(cmd, cwd=cwd, env=ENV, input=stdin, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)} failed in {cwd}:\n{(r.stdout + r.stderr)[-2000:]}")


def make_amber(work: Path) -> None:
    work.mkdir(parents=True, exist_ok=True)
    (work / "leap.in").write_text(LEAP)
    (work / "min.in").write_text(SANDER_MIN)
    (work / "md.in").write_text(SANDER_MD)
    run(["tleap", "-f", "leap.in"], work)
    run(["sander", "-O", "-i", "min.in", "-o", "min.out", "-p", "cln025.prmtop", "-c", "cln025.inpcrd",
         "-r", "min.rst7"], work)
    run(["sander", "-O", "-i", "md.in", "-o", "mdout", "-p", "cln025.prmtop", "-c", "min.rst7",
         "-r", "md.rst7", "-x", "md.nc", "-inf", "mdinfo"], work)


def make_gromacs(work: Path, solute_pdb: Path) -> None:
    """Parts: part0001 steps 0-1000, part0002 one frame at 1000 (ChainReader trap, duplicate boundary frame), part0003 to 2500."""
    work.mkdir(parents=True, exist_ok=True)
    shutil.copy(solute_pdb, work / "solute.pdb")
    (work / "md.mdp").write_text(MDP)
    (work / "em.mdp").write_text(EM_MDP)
    (work / "ions.mdp").write_text("integrator = steep\nnsteps = 0\n")
    run(["gmx", "pdb2gmx", "-f", "solute.pdb", "-o", "conf.gro", "-p", "topol.top", "-ff", "amber99sb-ildn",
         "-water", "tip3p", "-ignh"], work)
    run(["gmx", "editconf", "-f", "conf.gro", "-o", "box.gro", "-c", "-d", "0.8", "-bt", "dodecahedron"], work)
    run(["gmx", "solvate", "-cp", "box.gro", "-cs", "spc216.gro", "-p", "topol.top", "-o", "sol.gro"], work)
    run(["gmx", "grompp", "-f", "ions.mdp", "-c", "sol.gro", "-p", "topol.top", "-o", "ions.tpr", "-maxwarn", "2"], work)
    run(["gmx", "genion", "-s", "ions.tpr", "-o", "ion.gro", "-p", "topol.top", "-pname", "NA", "-nname", "CL",
         "-neutral"], work, stdin="SOL\n")
    run(["gmx", "grompp", "-f", "em.mdp", "-c", "ion.gro", "-p", "topol.top", "-o", "em.tpr"], work)
    run(["gmx", "mdrun", "-deffnm", "em", "-nt", "4"], work)
    run(["gmx", "grompp", "-f", "md.mdp", "-c", "em.gro", "-p", "topol.top", "-o", "md.tpr", "-po", "mdout.mdp"], work)
    run(["gmx", "mdrun", "-s", "md.tpr", "-deffnm", "md", "-noappend", "-nt", "4", "-cpt", "1", "-nsteps", "1000"], work)
    for nsteps in ("49", "1451"):   # -nsteps on a restart counts additional steps
        run(["gmx", "mdrun", "-deffnm", "md", "-cpi", "md.cpt", "-noappend", "-nt", "4", "-nsteps", nsteps], work)


def make_openmm(work: Path, prmtop: Path, inpcrd: Path) -> None:
    """LangevinMiddle 4 fs + HMR 3 amu + MonteCarloBarostat; DCD + StateDataReporter csv + system/integrator xml."""
    import openmm as mm
    from openmm import app, unit
    work.mkdir(parents=True, exist_ok=True)
    top, crd = app.AmberPrmtopFile(str(prmtop)), app.AmberInpcrdFile(str(inpcrd))
    system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=0.9 * unit.nanometer,
                              constraints=app.HBonds, hydrogenMass=3 * unit.amu)
    system.addForce(mm.MonteCarloBarostat(1 * unit.bar, 300 * unit.kelvin, 25))
    integrator = mm.LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.004 * unit.picoseconds)
    integrator.setRandomNumberSeed(12345)
    sim = app.Simulation(top.topology, system, integrator, mm.Platform.getPlatformByName("CPU"))
    sim.context.setPositions(crd.positions)
    sim.context.setPeriodicBoxVectors(*crd.boxVectors)
    sim.minimizeEnergy(maxIterations=200)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, 12345)
    with open(work / "topology.pdb", "w") as f:
        app.PDBFile.writeFile(sim.topology, sim.context.getState(getPositions=True).getPositions(), f)
    sim.reporters.append(app.DCDReporter(str(work / "traj.dcd"), 50))
    sim.reporters.append(app.StateDataReporter(
        str(work / "state.csv"), 25, step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, volume=True, density=True, separator=","))
    sim.step(1500)
    (work / "system.xml").write_text(mm.XmlSerializer.serialize(system))
    (work / "integrator.xml").write_text(mm.XmlSerializer.serialize(integrator))


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).parent / "data")
    make_amber(out / "amber")
    make_gromacs(out / "gromacs", out / "amber" / "solute.pdb")
    make_openmm(out / "openmm", out / "amber" / "cln025.prmtop", out / "amber" / "cln025.inpcrd")
