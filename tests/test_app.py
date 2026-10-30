"""Headless window test (QT_QPA_PLATFORM=offscreen): run an analysis, edit a rule, check the report without recompute."""
import os
import time

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QTWEBENGINE_DISABLE_SANDBOX", "1")
pytest.importorskip("PySide6.QtWebEngineWidgets")

from PySide6.QtCore import QEventLoop  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from trajinspector import rules as R  # noqa: E402
from trajinspector.app import MainWindow, parse_value  # noqa: E402


def pump(app, until, timeout=120.0):
    end = time.time() + timeout
    while not until() and time.time() < end:
        app.processEvents(QEventLoop.AllEvents, 50)
        time.sleep(0.02)   # a busy loop would starve the analysis thread of the GIL
    assert until(), "timed out"


@pytest.fixture(scope="module")
def window(data_dir, small_rules):
    app = QApplication.instance() or QApplication([])
    w = MainWindow(small_rules)
    for engine in ("gromacs", "amber", "openmm"):
        w.add_folder(str(data_dir / engine))
    w.reference.setText(str(data_dir / "reference_compact.pdb"))
    w.run_analysis()
    pump(app, lambda: w.campaign is not None or "Failed" in w.status.text())
    assert w.campaign is not None, w.status.text()
    pump(app, lambda: w.thread.isFinished())
    yield w
    w.tmp.cleanup()


def test_run_fills_report_viewer_and_rules(window):
    html = (open(os.path.join(window.tmp.name, "report.html")).read())
    assert "Gate 3 - Sampling" in html and "Done" in window.status.text()
    assert os.path.exists(os.path.join(window.tmp.name, "viewer", "gromacs.pdb"))
    assert window.rules_panel.table.rowCount() == len(R.load_rules(None).rules)


def test_rule_edit_changes_report_without_recompute(window):
    app = QApplication.instance()
    campaign, thread = window.campaign, window.thread
    path = os.path.join(window.tmp.name, "report.html")
    before = open(path).read()
    row = [window.rules_panel.table.item(r, 2).text() for r in range(window.rules_panel.table.rowCount())].index("cutoff_min_nm")
    window.rules_panel.table.item(row, 5).setText("0.7")
    window.rules_panel.apply()
    pump(app, lambda: open(path).read() != before)
    assert window.campaign is campaign and window.thread is thread   # nothing was recomputed
    assert window.rules_panel.base.get("cutoff_min_nm").value == 0.7
    import json
    from trajinspector import report
    rows = [r for r in json.loads(report.export_qc_json(campaign, window.rules_panel.base))["rows"] if r["rule_id"] == "cutoff_min_nm"]
    assert [r["status"] for r in rows] == ["PASS"] * 3   # amber 0.8 nm was WARN against 0.9


@pytest.mark.parametrize("kind, text, ok", [("max", "2.5", True), ("max", "abc", False), ("range", "[0.8, 1.2]", True),
                                            ("range", "[1]", False), ("bool", "true", True), ("bool", "1", False),
                                            ("not_in", "[berendsen]", True)])
def test_parse_value(kind, text, ok):
    rule = R.replace(R.load_rules(None).get("thermostat"), kind=kind)
    if ok:
        parse_value(rule, text)
    else:
        with pytest.raises(ValueError, match="not valid"):
            parse_value(rule, text)
