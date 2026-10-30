import json
import re

from trajinspector import report
from trajinspector import rules as R


def test_html_is_self_contained(campaign, small_rules):
    html = report.render_html(campaign, small_rules)
    assert not re.search(r"<script|<link|@import|src=[\"']http|url\(http", html)   # nothing fetched at runtime
    for title in report.GATES.values():
        assert title in html
    assert "rule of thumb" in html and "Not in v1" in html and "water model is not recorded" in html
    assert html.count("<svg") >= 4


def test_gate_card_counts_match_rows(campaign, small_rules):
    d = json.loads(report.export_qc_json(campaign, small_rules))
    assert len(d["rows"]) == len(campaign.results)
    for gate, title in report.GATES.items():
        n = sum(1 for r in d["rows"] if r["gate"] == gate)
        assert sum(d["summary"][title].values()) == n
    assert {r["status"] for r in d["rows"]} >= {"PASS", "WARN", "FAIL", "N/A", "MANUAL", "INFO"}


def test_editing_a_rule_changes_the_report(campaign, small_rules):
    edited = R.RuleSet(small_rules.version, small_rules.params, tuple(
        R.replace(r, enabled=False) if r.id == "thermostat" else r for r in small_rules.rules))
    a = json.loads(report.export_qc_json(campaign, small_rules))
    b = json.loads(report.export_qc_json(campaign, edited))
    pick = lambda d: [r["status"] for r in d["rows"] if r["rule_id"] == "thermostat"]   # noqa: E731
    assert pick(a) == ["PASS"] * 3 and pick(b) == ["N/A"] * 3


def test_qc_json_is_strict_json(campaign, small_rules):
    json.loads(report.export_qc_json(campaign, small_rules), parse_constant=lambda c: (_ for _ in ()).throw(ValueError(c)))
