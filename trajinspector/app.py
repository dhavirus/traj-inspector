"""Qt window: replica list and Run on the left; Report, Viewer and Rules tabs on the right."""
from __future__ import annotations

import os
import sys
import tempfile
import traceback
from dataclasses import replace
from pathlib import Path

import yaml

# Software GL fallback so the NGL viewer also works on a PC without a usable GPU. Must be set before QtWebEngine starts.
os.environ.setdefault("QTWEBENGINE_CHROMIUM_FLAGS", "--ignore-gpu-blocklist --enable-unsafe-swiftshader")

from PySide6.QtCore import QEvent, QEventLoop, Qt, QThread, QUrl, Signal  # noqa: E402
from PySide6.QtGui import QAction  # noqa: E402
from PySide6.QtWebEngineCore import QWebEngineSettings  # noqa: E402
from PySide6.QtWebEngineWidgets import QWebEngineView  # noqa: E402
from PySide6.QtWidgets import (QApplication, QComboBox, QFileDialog, QHBoxLayout, QLabel, QLineEdit, QListWidget,  # noqa: E402
                               QMainWindow, QMessageBox, QProgressBar, QPushButton, QSplitter, QTableWidget,
                               QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from . import checks, readers, report, rules as R  # noqa: E402
from .viewer import export_viewer  # noqa: E402

EDITABLE_KINDS = ("max", "min", "range", "not_in", "bool")
RULE_COLUMNS = ("On", "Gate", "Id", "Check", "Kind", "Value", "Severity", "Source", "Rule of thumb")


def parse_value(rule: R.Rule, text: str):
    """Typed value for `rule` from the text typed in the table; ValueError with a message if it does not fit."""
    value = yaml.safe_load(text)
    numeric = lambda x: isinstance(x, (int, float)) and not isinstance(x, bool)   # noqa: E731
    ok = {"max": numeric(value), "min": numeric(value), "bool": isinstance(value, bool),
          "range": isinstance(value, list) and len(value) == 2 and all(map(numeric, value)),
          "not_in": isinstance(value, list) and all(isinstance(v, str) for v in value)}[rule.kind]
    if not ok:
        hint = {"max": "a number", "min": "a number", "bool": "true or false", "range": "[low, high]", "not_in": "[name, ...]"}[rule.kind]
        raise ValueError(f"{rule.id}: '{text}' is not valid for kind '{rule.kind}' (expected {hint})")
    return value


class AnalysisThread(QThread):
    progress = Signal(str)
    done = Signal(object, object)   # Campaign, viewer.html path
    failed = Signal(str)

    def __init__(self, folders, reference, selection, ruleset, out_dir):
        super().__init__()
        self.args = (folders, reference, selection, ruleset, out_dir)

    def run(self) -> None:
        folders, reference, selection, ruleset, out_dir = self.args
        try:
            replicas = []
            for folder in folders:
                self.progress.emit(f"reading {Path(folder).name}")
                replicas.append(readers.read_replica(folder, selection))
            campaign = checks.run_checks(replicas, ruleset, selection, reference, self.progress.emit)
            self.progress.emit("preparing the viewer")
            self.done.emit(campaign, export_viewer(campaign, selection, ruleset.params, Path(out_dir) / "viewer"))
        except Exception:   # reported in the window; a silent thread death would leave the UI stuck on "running"
            self.failed.emit(traceback.format_exc())


def web_view() -> QWebEngineView:
    view = QWebEngineView()
    settings = view.settings()
    settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessFileUrls, True)
    settings.setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls, False)
    settings.setAttribute(QWebEngineSettings.WebAttribute.WebGLEnabled, True)
    return view


class RulesPanel(QWidget):
    applied = Signal()

    def __init__(self, ruleset: R.RuleSet):
        super().__init__()
        self.base = ruleset
        self.table = QTableWidget(0, len(RULE_COLUMNS))
        self.table.setHorizontalHeaderLabels(RULE_COLUMNS)
        self.params = QTableWidget(0, 2)
        self.params.setHorizontalHeaderLabels(("Parameter (changes need Run)", "Value"))
        buttons = QHBoxLayout()
        for text, slot in (("Apply", self.apply), ("Save", self.save), ("Reset to defaults", self.reset)):
            button = QPushButton(text)
            button.clicked.connect(slot)
            buttons.addWidget(button)
        buttons.addStretch()
        layout = QVBoxLayout(self)
        layout.addLayout(buttons)
        layout.addWidget(self.table, 3)
        layout.addWidget(self.params, 1)
        self.fill(ruleset)

    def fill(self, ruleset: R.RuleSet) -> None:
        self.base = ruleset
        self.table.setRowCount(len(ruleset.rules))
        for row, rule in enumerate(ruleset.rules):
            on = QTableWidgetItem()
            on.setFlags(Qt.ItemIsUserCheckable | Qt.ItemIsEnabled)
            on.setCheckState(Qt.Checked if rule.enabled else Qt.Unchecked)
            cells = (on, str(rule.gate), rule.id, rule.check, rule.kind, "" if rule.value is None else yaml.safe_dump(rule.value, default_flow_style=True).strip().removesuffix("\n..."),
                     None, rule.source or "—", "yes" if rule.rule_of_thumb else "")
            for col, cell in enumerate(cells):
                if cell is None:
                    combo = QComboBox()
                    combo.addItems(R.SEVERITIES)
                    combo.setCurrentText(rule.severity)
                    self.table.setCellWidget(row, col, combo)
                    continue
                item = cell if isinstance(cell, QTableWidgetItem) else QTableWidgetItem(cell)
                if col not in (0, 5) or (col == 5 and rule.kind not in EDITABLE_KINDS):
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(row, col, item)
        self.table.resizeColumnsToContents()
        self.params.setRowCount(len(ruleset.params))
        for row, (key, value) in enumerate(ruleset.params.items()):
            name = QTableWidgetItem(key)
            name.setFlags(name.flags() & ~Qt.ItemIsEditable)
            self.params.setItem(row, 0, name)
            self.params.setItem(row, 1, QTableWidgetItem(str(value)))
        self.params.resizeColumnsToContents()

    def ruleset(self) -> R.RuleSet:
        """The rule set as shown in the tables. Raises ValueError on a value that does not fit its rule kind."""
        edited = []
        for row, rule in enumerate(self.base.rules):
            value = rule.value if rule.kind not in EDITABLE_KINDS else parse_value(rule, self.table.item(row, 5).text())
            edited.append(replace(rule, enabled=self.table.item(row, 0).checkState() == Qt.Checked, value=value,
                                  severity=self.table.cellWidget(row, 6).currentText()))
        params = {self.params.item(r, 0).text(): yaml.safe_load(self.params.item(r, 1).text()) for r in range(self.params.rowCount())}
        return R.RuleSet(self.base.version, params, tuple(edited))

    def _checked(self) -> R.RuleSet | None:
        try:
            return self.ruleset()
        except ValueError as exc:
            QMessageBox.warning(self, "Invalid value", str(exc))
            return None

    def apply(self) -> None:
        if (ruleset := self._checked()) is not None:
            self.base = ruleset
            self.applied.emit()

    def save(self) -> None:
        if (ruleset := self._checked()) is not None:
            R.save_rules(ruleset)
            self.base = ruleset
            QMessageBox.information(self, "Saved", f"Rules written to {R.USER_RULES}")

    def reset(self) -> None:
        if QMessageBox.question(self, "Reset", f"Restore the shipped defaults and remove {R.USER_RULES}?") == QMessageBox.Yes:
            R.USER_RULES.unlink(missing_ok=True)
            self.fill(R.load_rules(None))
            self.applied.emit()


class MainWindow(QMainWindow):
    def __init__(self, ruleset: R.RuleSet | None = None):
        super().__init__()
        self.setWindowTitle("Traj Inspector")
        self.resize(1400, 850)
        self.tmp = tempfile.TemporaryDirectory(prefix="trajinspector-")
        self.campaign: checks.Campaign | None = None
        self.thread: AnalysisThread | None = None
        self.folders = QListWidget()
        self.reference = QLineEdit()
        self.reference.setPlaceholderText("optional: folded reference (pdb/gro/cif); default frame 0 of replica 1 (WARN)")
        self.selection = QLineEdit("protein")
        self.run_button = QPushButton("Run")
        self.progress = QProgressBar()
        self.progress.setRange(0, 1)
        self.status = QLabel("Add replica folders, then Run.")
        self.status.setWordWrap(True)
        self.report_view, self.viewer_view = web_view(), web_view()
        self.rules_panel = RulesPanel(ruleset or R.load_rules())
        self.rules_panel.applied.connect(self.render_report)
        self._layout()
        self._menu()
        self.run_button.clicked.connect(self.run_analysis)

    def _layout(self) -> None:
        add, remove, browse = QPushButton("Add replica folder…"), QPushButton("Remove"), QPushButton("Browse…")
        add.clicked.connect(self.add_folder)
        remove.clicked.connect(lambda: [self.folders.takeItem(self.folders.row(i)) for i in self.folders.selectedItems()])
        browse.clicked.connect(self.pick_reference)
        left = QWidget()
        box = QVBoxLayout(left)
        row = QHBoxLayout()
        row.addWidget(add)
        row.addWidget(remove)
        ref_row = QHBoxLayout()
        ref_row.addWidget(self.reference)
        ref_row.addWidget(browse)
        for widget in (QLabel("Replicas"), self.folders):
            box.addWidget(widget)
        box.addLayout(row)
        box.addWidget(QLabel("Reference structure"))
        box.addLayout(ref_row)
        box.addWidget(QLabel("Atom selection"))
        for widget in (self.selection, self.run_button, self.progress, self.status):
            box.addWidget(widget)
        box.addStretch()
        tabs = QTabWidget()
        tabs.addTab(self.report_view, "Report")
        tabs.addTab(self.viewer_view, "Viewer")
        tabs.addTab(self.rules_panel, "Rules")
        self.tabs = tabs
        splitter = QSplitter()
        splitter.addWidget(left)
        splitter.addWidget(tabs)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([340, 1060])
        self.setCentralWidget(splitter)

    def _menu(self) -> None:
        menu = self.menuBar().addMenu("&File")
        for text, slot in (("Export HTML report…", self.export_html), ("Export qc.json…", self.export_json)):
            action = QAction(text, self)
            action.triggered.connect(slot)
            menu.addAction(action)

    # -------------------------------------------------------------- inputs
    def add_folder(self, folder: str | None = None) -> None:
        folder = folder or QFileDialog.getExistingDirectory(self, "Replica folder")
        if not folder:
            return
        try:
            engine = readers.detect_engine(Path(folder))
        except readers.ReaderError as exc:
            QMessageBox.warning(self, "Not a replica folder", str(exc))
            return
        self.folders.addItem(f"{folder}   [{engine}]")

    def pick_reference(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Reference structure", "", "Structures (*.pdb *.gro *.cif *.ent)")
        if path:
            self.reference.setText(path)

    def replica_paths(self) -> list[str]:
        return [self.folders.item(i).text().rsplit("   [", 1)[0] for i in range(self.folders.count())]

    # -------------------------------------------------------------- run
    def run_analysis(self) -> None:
        if not self.folders.count():
            self.status.setText("Add at least one replica folder.")
            return
        try:
            ruleset = self.rules_panel.ruleset()
        except ValueError as exc:
            self.status.setText(str(exc))
            return
        self.rules_panel.base = ruleset
        self.run_button.setEnabled(False)
        self.progress.setRange(0, 0)
        self.thread = AnalysisThread(self.replica_paths(), self.reference.text().strip() or None, self.selection.text().strip(),
                                     ruleset, self.tmp.name)
        self.thread.progress.connect(self.status.setText)
        self.thread.done.connect(self.on_done)
        self.thread.failed.connect(self.on_failed)
        self.thread.start()

    def _finish(self) -> None:
        self.run_button.setEnabled(True)
        self.progress.setRange(0, 1)

    def on_done(self, campaign: checks.Campaign, viewer_html: Path) -> None:
        self._finish()
        self.campaign = campaign
        self.render_report()
        self.viewer_view.load(QUrl.fromLocalFile(str(viewer_html)))
        self.status.setText(f"Done: {len(campaign.results)} results.")

    def on_failed(self, trace: str) -> None:
        self._finish()
        self.status.setText("Failed: " + trace.strip().splitlines()[-1])
        QMessageBox.critical(self, "Analysis failed", trace)

    def render_report(self) -> None:
        """Re-evaluate the stored metrics with the current rules; nothing is recomputed."""
        if self.campaign is None:
            return
        path = Path(self.tmp.name) / "report.html"
        path.write_text(report.render_html(self.campaign, self.rules_panel.base))   # file, not setHtml: setHtml is capped at 2 MB
        self.report_view.load(QUrl.fromLocalFile(str(path)))

    # -------------------------------------------------------------- export
    def _export(self, title: str, pattern: str, make) -> None:
        if self.campaign is None:
            QMessageBox.information(self, title, "Run an analysis first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, title, "", pattern)
        if path:
            Path(path).write_text(make(self.campaign, self.rules_panel.base))

    def export_html(self) -> None:
        self._export("Export HTML report", "HTML (*.html)", lambda c, r: report.render_html(c, r))

    def export_json(self) -> None:
        self._export("Export qc.json", "JSON (*.json)", lambda c, r: report.export_qc_json(c, r))


def self_check() -> None:
    """Build the window without showing it: proves Qt, QtWebEngine, the rules file and the bundled NGL are all present."""
    from . import __version__
    from .viewer import VIEWER_DIR
    window = MainWindow()
    print(f"trajinspector {__version__}: {len(window.rules_panel.base.rules)} rules, ngl.js {(VIEWER_DIR / 'ngl.js').stat().st_size} bytes, ok")
    for view in (window.report_view, window.viewer_view):   # pages must die before the WebEngine profile at interpreter exit
        view.page().deleteLater()
        view.deleteLater()
    window.deleteLater()
    QApplication.processEvents(QEventLoop.AllEvents)
    QApplication.sendPostedEvents(None, QEvent.DeferredDelete)


def main() -> None:
    os.environ.setdefault("QTWEBENGINE_DISABLE_SANDBOX", "1")   # also set by AppRun: the sandbox cannot be setuid inside an AppImage
    app = QApplication(sys.argv)
    if "--self-check" in sys.argv:
        self_check()
        return
    window = MainWindow()
    window.show()
    sys.exit(app.exec())
