"""Viewer export: selection-only multi-model PDB per replica plus per-frame time / Q / state for viewer.html."""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import MDAnalysis as mda
import numpy as np

from .. import analysis as A
from ..checks import Campaign, prepare

VIEWER_DIR = Path(__file__).parent
MAX_VIEWER_FRAMES = 500   # keeps the multi-model PDB small enough to load instantly
ASSETS = ("viewer.html", "ngl.js")


def _pick(n_frames: int) -> np.ndarray:
    return np.unique(np.linspace(0, n_frames - 1, min(n_frames, MAX_VIEWER_FRAMES)).astype(int))


def export_viewer(campaign: Campaign, selection: str, params: dict, out_dir: Path) -> Path:
    """Write viewer files into out_dir and return the path of viewer.html. Frames are centred on the selection."""
    out_dir.mkdir(parents=True, exist_ok=True)
    missing = [a for a in ASSETS if not (VIEWER_DIR / a).exists()]
    if missing:
        raise FileNotFoundError(f"{missing} missing in {VIEWER_DIR}; run packaging/fetch_ngl.sh")
    for asset in ASSETS:
        shutil.copy(VIEWER_DIR / asset, out_dir / asset)
    replicas = []
    for rep in campaign.replicas:
        series = campaign.series[rep.name]
        prep = prepare(rep, selection, params["max_frames"])
        states = A.count_transitions(series.q, params["q_unfolded"], params["q_folded"])["state"]
        picks = _pick(len(prep.frames))
        with mda.Writer(str(out_dir / f"{rep.name}.pdb"), len(prep.group), multiframe=True) as writer:
            for i in picks:
                prep.universe.trajectory[prep.frames[i]]
                prep.group.positions = prep.group.positions - prep.group.center_of_mass()
                writer.write(prep.group)
        replicas.append({"name": rep.name, "pdb": f"{rep.name}.pdb", "time_ps": series.time_ps[picks].tolist(),
                         "q": series.q[picks].tolist(), "state": states[picks].tolist()})
    (out_dir / "viewer_data.js").write_text("window.VIEWER_DATA = " + json.dumps({"replicas": replicas}) + ";\n")
    return out_dir / "viewer.html"
