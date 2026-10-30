from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA


@pytest.fixture(scope="session")
def small_rules():
    """Defaults, with the run-length gates relaxed: the 5-6 ps test runs would be N/A otherwise."""
    from trajinspector import rules as R
    rs = R.load_rules(None)
    params = {**rs.params, "msm_lag_ns": 0.0005, "drift_min_ns": 0.002, "ke_min_frames": 50}
    return R.RuleSet(rs.version, params, rs.rules)


@pytest.fixture(scope="session")
def campaign(data_dir, small_rules):
    from trajinspector import checks, readers
    reps = [readers.read_replica(data_dir / e) for e in ("gromacs", "amber", "openmm")]
    return checks.run_checks(reps, small_rules, reference=data_dir / "reference_compact.pdb")
