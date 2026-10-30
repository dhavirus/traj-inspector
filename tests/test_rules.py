import pytest

from trajinspector import rules as R


def test_defaults_load_and_have_sources():
    rs = R.load_rules(None)
    assert len(rs.rules) == 28
    for rule in rs.rules:
        assert rule.rule_of_thumb or rule.source, f"{rule.id}: no source and not flagged rule_of_thumb"


@pytest.mark.parametrize("rule_id, metric, status", [
    ("thermostat", "berendsen", R.FAIL), ("thermostat", "v-rescale", R.PASS),
    ("barostat", "parrinello-rahman", R.WARN), ("barostat", "berendsen", R.FAIL),
    ("cutoff_min_nm", 0.9, R.PASS), ("cutoff_min_nm", 0.89, R.WARN),
    ("ke_width", 0.85, R.PASS), ("ke_width", 1.16, R.FAIL),
    ("hmr_rates", True, R.WARN), ("hmr_rates", False, R.PASS),
    ("forcefield_water", None, R.MANUAL), ("vamp2", 3.1, R.INFO), ("drift", None, R.NA),
    ("drift", float("nan"), R.NA),
])
def test_evaluate(rule_id, metric, status):
    assert R.evaluate(R.load_rules(None).get(rule_id), metric) == status


def test_user_override_changes_status_and_roundtrips(tmp_path):
    rs = R.load_rules(None)
    edited = R.RuleSet(rs.version, rs.params, tuple(
        R.replace(r, value=0.5, severity="fail") if r.id == "cutoff_min_nm" else r for r in rs.rules))
    path = tmp_path / "rules.yaml"
    R.save_rules(edited, path)
    again = R.load_rules(path)
    assert again.get("cutoff_min_nm").value == 0.5 and again.get("cutoff_min_nm").severity == "fail"
    assert R.evaluate(again.get("cutoff_min_nm"), 0.4) == R.FAIL
    assert again.get("barostat").extra == {"warn_in": ["parrinello-rahman"]}   # extra keys survive the round trip


def test_disabled_rule_is_na():
    rule = R.replace(R.load_rules(None).get("thermostat"), enabled=False)
    assert R.evaluate(rule, "berendsen") == R.NA


def test_malformed_rule_rejected(tmp_path):
    bad = tmp_path / "r.yaml"
    bad.write_text("version: 1\nparams: {}\nrules:\n- {id: x, gate: 1, check: c, kind: nope, value: 1, severity: fail,"
                   " source: s, tested_on: t, rule_of_thumb: false, enabled: true}\n")
    with pytest.raises(R.RulesError, match="unknown kind"):
        R.load_rules(bad)
