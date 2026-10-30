"""Readers on real CLN025 data. Expected values come from the inputs (mdp, echoed sander input, system.xml,
make_testdata.py recipe), not from reader output."""
import os

import MDAnalysis as mda
import numpy as np
import pytest

from trajinspector import readers as R


@pytest.fixture(scope="module")
def gmx(data_dir):
    return R.read_replica(data_dir / "gromacs")


@pytest.fixture(scope="module")
def amber(data_dir):
    return R.read_replica(data_dir / "amber")


@pytest.fixture(scope="module")
def omm(data_dir):
    return R.read_replica(data_dir / "openmm")


def test_engine_detection(data_dir, gmx, amber, omm):
    assert (gmx.engine, amber.engine, omm.engine) == ("gromacs", "amber", "openmm")


def test_gromacs_setup(gmx):
    s = gmx.setup
    assert (s.thermostat, s.barostat, s.dt_ps, s.constraints, s.cutoff_nm) == ("v-rescale", "c-rescale", 0.002, "h-bonds", 1.0)
    assert s.hmr_factor == pytest.approx(1.0, abs=1e-3)
    assert s.warnings_count == 0


def test_gromacs_energies_concatenated_without_duplicates(gmx):
    t = gmx.energies.time_ps
    assert np.all(np.diff(t) > 0)
    assert t[0] == 0.0 and t[-1] == pytest.approx(5.0)
    assert len(t) >= 51   # 0..5 ps every 0.1 ps; the one-frame part adds the final step
    assert gmx.energies.conserved is not None and not np.isnan(gmx.energies.conserved).any()
    assert 900 < gmx.energies.density_kg_m3.mean() < 1100


def test_gromacs_chain_with_one_frame_part(gmx):
    sizes = [mda.Universe(str(gmx.topology), str(p)).trajectory.n_frames for p in gmx.trajectories]
    assert sizes == [21, 1, 30]
    with pytest.raises(RuntimeError, match="at least two frames"):   # why open_universe does not use continuous=True
        mda.Universe(str(gmx.topology), [str(p) for p in gmx.trajectories], continuous=True)
    universe, frames = R.open_universe(gmx)
    times = np.array([universe.trajectory[i].time for i in frames])
    assert len(frames) == 51 and np.all(np.diff(times) > 0)
    assert len(universe.select_atoms("protein")) == 166


def test_amber_setup(amber):
    s = amber.setup
    assert (s.thermostat, s.barostat, s.dt_ps, s.constraints) == ("langevin", "monte-carlo", 0.002, "h-bonds")
    assert s.cutoff_nm == pytest.approx(0.8)   # cut = 8.0 A
    assert s.hmr_factor == pytest.approx(1.0, abs=1e-3) and s.warnings_count == 0


def test_amber_energies(amber):
    e = amber.energies
    assert len(e.time_ps) == 121 and e.time_ps[-1] == pytest.approx(6.0)   # nstlim 3000, ntpr 25
    assert e.conserved is None   # Langevin: Etot is not conserved
    assert e.density_kg_m3[0] == pytest.approx(712.3, abs=0.1)   # mdout: Density = 0.7123 g/cc
    assert e.volume[0] == pytest.approx(78.676, abs=1e-3)        # mdout: VOLUME = 78676.1190 A^3
    assert e.kinetic[0] == pytest.approx(5086.4908 * 4.184, rel=1e-6)   # first EKtot in mdout, kcal -> kJ


def test_amber_trajectory(amber):
    universe, frames = R.open_universe(amber)
    assert len(frames) == 60 and len(universe.select_atoms("protein")) == 166


def test_openmm_setup(omm):
    s = omm.setup
    assert (s.thermostat, s.barostat, s.dt_ps, s.constraints, s.cutoff_nm) == ("langevin", "monte-carlo", 0.004, "h-bonds", 0.9)
    assert s.hmr_factor == pytest.approx(3.0 / 1.008, rel=1e-3)   # hydrogenMass=3 amu
    assert s.warnings_count is None


def test_openmm_energies(omm):
    e = omm.energies
    assert len(e.time_ps) == 60 and e.time_ps[-1] == pytest.approx(6.0)   # 1500 steps x 4 fs, report every 25
    assert e.conserved is None
    assert e.density_kg_m3[0] == pytest.approx(712.308, abs=0.01)   # state.csv row 1: 0.7123084 g/mL


def test_openmm_universe_keeps_hmr_masses(omm):
    universe, frames = R.open_universe(omm)
    assert len(frames) == 30
    assert universe.select_atoms("protein and name H*").masses.mean() == pytest.approx(3.0)


@pytest.fixture
def gmx_dir(tmp_path, data_dir):
    for p in (data_dir / "gromacs").iterdir():
        if p.suffix != ".mdp":
            os.symlink(p, tmp_path / p.name)
    return tmp_path


@pytest.mark.parametrize("mdp, thermostat, barostat", [
    ("integrator = sd\ndt = 0.004\ntcoupl = no\npcoupl = parrinello-rahman\nconstraints = h-bonds\n", "langevin", "parrinello-rahman"),
    ("integrator = md\ndt = 0.002\ntcoupl = Berendsen\npcoupl = Berendsen\nrvdw = 1.2 ; comment\nrcoulomb = 0.9\n", "berendsen", "berendsen"),
])
def test_gromacs_mdp_mapping(gmx_dir, mdp, thermostat, barostat):
    (gmx_dir / "md.mdp").write_text(mdp)
    s = R.read_replica(gmx_dir).setup
    assert (s.thermostat, s.barostat) == (thermostat, barostat)


def test_gromacs_mdp_cutoff_is_larger_of_vdw_and_coulomb(gmx_dir):
    (gmx_dir / "md.mdp").write_text("dt = 0.002\nrvdw = 1.2 ; comment\nrcoulomb = 0.9\n")
    assert R.read_replica(gmx_dir).setup.cutoff_nm == 1.2


def test_gromacs_lincs_warnings_counted(gmx_dir):
    (gmx_dir / "md.mdp").write_text("dt = 0.002\n")
    (gmx_dir / "md.part0001.log").unlink()
    (gmx_dir / "md.part0001.log").write_text("Step 5  LINCS WARNING\nrelative constraint deviation\nStep 9  LINCS WARNING\n")
    assert R.read_replica(gmx_dir).setup.warnings_count == 2   # the 'LINCS:' citation line in other logs is not counted


def test_empty_part_is_skipped(gmx_dir):
    (gmx_dir / "md.mdp").write_text("dt = 0.002\n")
    (gmx_dir / "md.part0004.xtc").write_bytes(b"")
    assert len(R.read_replica(gmx_dir).trajectories) == 3


def test_unknown_folder_raises(tmp_path):
    with pytest.raises(R.ReaderError, match="no .tpr"):
        R.read_replica(tmp_path)
