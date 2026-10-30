"""Headless screenshots of the three tabs: QT_QPA_PLATFORM=offscreen python tests/make_screenshots.py OUTDIR
Looks at the real window with the test replicas loaded; not a pytest test (the images are for a human)."""
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QTWEBENGINE_DISABLE_SANDBOX", "1")

from PySide6.QtCore import QEventLoop, QTimer  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from trajinspector import rules as R  # noqa: E402
from trajinspector.app import MainWindow  # noqa: E402

DATA = Path(__file__).parent / "data"


def pump(app, seconds):
    end = time.time() + seconds
    while time.time() < end:
        app.processEvents(QEventLoop.AllEvents, 50)
        time.sleep(0.02)   # a busy loop would starve the analysis thread of the GIL


def main(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    app = QApplication([])
    rules = R.load_rules(None)
    window = MainWindow(R.RuleSet(rules.version, {**rules.params, "msm_lag_ns": 0.0005}, rules.rules))
    window.show()
    for engine in ("gromacs", "amber", "openmm"):
        window.add_folder(str(DATA / engine))
    window.reference.setText(str(DATA / "reference_compact.pdb"))
    window.run_analysis()
    while window.thread.isRunning():
        pump(app, 0.2)
    pump(app, 4)
    for index, name in enumerate(("report", "viewer", "rules")):
        window.tabs.setCurrentIndex(index)
        pump(app, 5 if name != "rules" else 1)
        window.grab().save(str(out / f"{name}.png"))
    print("status:", window.status.text())


if __name__ == "__main__":
    main(Path(sys.argv[1]))
