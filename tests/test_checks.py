import json

import numpy as np
import pytest

from trajinspector import checks as C
from trajinspector import rules as R
from trajinspector import readers


def table(campaign, rules):
    return {(r.result.rule_id, r.result.replica): r for r in C.evaluate_results(rules, campaign.results)}


def test_every_enabled_rule_has_a_row(campaign, small_rules):
    seen = {r.rule_id for r in campaign.results}
    assert seen == {r.id for r in small_rules.rules}


def test_gate1_statuses_follow_the_input_files(campaign, small_rules):
    t = table(campaign, small_rules)
    # gromacs md.mdp: v-rescale/c-rescale, 2 fs, h-bonds, 1.0 nm; amber echo: langevin/MC barostat, cut=8.0 A; openmm xml: 4 fs HMR, 0.9 nm
    assert t["thermostat", "gromacs"].status == R.PASS and t["barostat", "amber"].status == R.PASS
    assert t["cutoff_min_nm", "amber"].status == R.WARN and t["cutoff_min_nm", "openmm"].status == R.PASS
    assert t["hmr_rates", "openmm"].status == R.WARN and t["hmr_rates", "gromacs"].status == R.PASS
    assert t["dt_max_fs", "openmm"].status == R.PASS and t["dt_max_fs", "openmm"].result.variant == "hmr"
    assert all(t["forcefield_water", e].status == R.MANUAL for e in ("gromacs", "amber", "openmm"))


def test_gate2_na_reasons(campaign, small_rules):
    t = table(campaign, small_rules)
    assert t["drift", "amber"].status == R.NA and "langevin" in t["drift", "amber"].result.reason
    assert t["warnings", "openmm"].status == R.NA
    assert t["drift", "gromacs"].result.metric is not None   # gromacs has Conserved En.
    assert t["min_image", "gromacs"].result.metric > 0 and t["rmsd_jump", "gromacs"].result.metric < 0.2


def test_ke_width_needs_enough_frames(data_dir, small_rules):
    strict = R.RuleSet(small_rules.version, {**small_rules.params, "ke_min_frames": 200}, small_rules.rules)
    c = C.run_checks([readers.read_replica(data_dir / "amber")], strict, reference=data_dir / "reference_compact.pdb")
    row = table(c, strict)["ke_width", "amber"]
    assert row.status == R.NA and "needs >= 200" in row.result.reason


def test_rule_edit_changes_status_without_recomputing(campaign, small_rules):
    before = table(campaign, small_rules)["cutoff_min_nm", "amber"]
    assert before.status == R.WARN
    edited = R.RuleSet(small_rules.version, small_rules.params, tuple(
        R.replace(r, value=0.7) if r.id == "cutoff_min_nm" else r for r in small_rules.rules))
    after = table(campaign, edited)["cutoff_min_nm", "amber"]
    assert after.status == R.PASS and after.result is before.result   # same stored metric


def test_campaign_rows_and_notes(campaign, small_rules):
    t = table(campaign, small_rules)
    assert t["replicas_min", "campaign"].status == R.PASS
    assert t["transitions", "campaign"].result.metric == 0 and t["transitions", "campaign"].status == R.FAIL
    assert any("no match in the trajectory" in n for n in campaign.notes)   # gromacs C-terminal O1/O2 vs O/OXT


def test_extended_frame0_reference_stops(data_dir, small_rules):
    with pytest.raises(ValueError, match="no native contacts"):
        C.run_checks([readers.read_replica(data_dir / "amber")], small_rules)


def test_frame0_reference_is_flagged(data_dir, small_rules, monkeypatch):
    compact = C.reference_from_file(data_dir / "reference_compact.pdb", "protein")
    monkeypatch.setattr(C, "reference_from_frame0", lambda prep: C.Reference(compact.keys, compact.xyz_angstrom, "frame 0"))
    c = C.run_checks([readers.read_replica(data_dir / "amber")], small_rules)
    assert any("No reference structure given" in n for n in c.notes)


def test_replica_agreement():
    agree = {"a": {"p": 0.5, "n_eff": 100.0}, "b": {"p": 0.52, "n_eff": 100.0}}
    differ = {"a": {"p": 0.2, "n_eff": 100.0}, "b": {"p": 0.8, "n_eff": 100.0}}
    # se = sqrt(.25/100) = .05 each -> sqrt(se_i^2 + se_j^2) = .0707; |dp| = .02 -> z .28; .6 -> z 8.5
    assert C._agreement(agree).metric == pytest.approx(0.02 / np.hypot(0.05, 0.05), rel=1e-2)
    assert C._agreement(differ).metric > 8
    assert C._agreement({"a": agree["a"]}).metric is None


def test_q_results_on_synthetic_two_state_series():
    rng = np.random.default_rng(0)
    state = np.cumsum(rng.random(4000) < 0.01) % 2        # telegraph process, mean dwell 100 frames
    q = np.where(state == 1, 0.95, 0.05) + rng.normal(0, 0.02, 4000)
    s = C.Series(np.arange(4000.0), q, np.ones(4000), np.zeros((4000, 1)), 1.0, 10)
    results, stats = C._q_results("syn", s, {"q_unfolded": 0.3, "q_folded": 0.8, "min_blocks": 8})
    by = {r.rule_id: r for r in results}
    assert by["q_bimodal"].metric == 0.0
    assert by["equilibration"].metric > 3000   # pymbar cuts nearly everything on a two-state series ...
    assert stats["transitions"] == int(np.count_nonzero(np.diff(state)))   # ... but the transitions are counted on all frames
    assert stats["transitions"] > 10 and by["n_eff_min"].metric > 20
