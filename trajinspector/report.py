"""Self-contained HTML report (inline SVG from matplotlib, no external resources) and qc.json export."""
from __future__ import annotations

import dataclasses
import io
import json
import math
from datetime import datetime, timezone
from html import escape

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from . import __version__  # noqa: E402
from . import rules as R  # noqa: E402
from .checks import CAMPAIGN, Campaign, Row, evaluate_results  # noqa: E402

GATES = {1: "Gate 1 - Setup", 2: "Gate 2 - Artifacts", 3: "Gate 3 - Sampling", 4: "Gate 4 - Mechanism"}
STATUS_ORDER = (R.FAIL, R.WARN, R.PASS, R.NA, R.MANUAL, R.INFO)
PALETTE = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000")   # Okabe-Ito
NOT_IN_V1 = "DSSP per state, committor test, force distribution analysis, group temperatures, gate 5 (experiment comparison)"
FES_BINS = 24

CSS = """
:root{--bg:#fff;--fg:#1b1f24;--muted:#5a6572;--line:#d8dee4;--card:#f6f8fa;--pass:#1a7f37;--warn:#9a6700;--fail:#cf222e;--info:#0969da;--na:#6e7781}
@media (prefers-color-scheme:dark){:root{--bg:#0d1117;--fg:#e6edf3;--muted:#9da7b3;--line:#30363d;--card:#161b22;--pass:#3fb950;--warn:#d29922;--fail:#f85149;--info:#58a6ff;--na:#8b949e}}
body{font:14px/1.45 system-ui,sans-serif;background:var(--bg);color:var(--fg);margin:0;padding:16px 24px;max-width:1280px}
h1{font-size:20px}h2{font-size:16px;margin:28px 0 8px}.muted{color:var(--muted)}
.cards{display:flex;gap:12px;flex-wrap:wrap}.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px 14px;min-width:200px}
.card b{display:block;margin-bottom:4px}.card span{margin-right:10px}
table{border-collapse:collapse;width:100%}th,td{border-bottom:1px solid var(--line);padding:5px 8px;text-align:left;vertical-align:top}
th{font-size:12px;color:var(--muted)}td.num{font-variant-numeric:tabular-nums}
.chip{display:inline-block;min-width:54px;text-align:center;border-radius:10px;padding:1px 8px;font-size:12px;font-weight:600;color:#fff}
.PASS{background:var(--pass)}.WARN{background:var(--warn)}.FAIL{background:var(--fail)}.INFO{background:var(--info)}.NA,.MANUAL{background:var(--na)}
.rot{font-size:11px;border:1px solid var(--line);border-radius:4px;padding:0 4px;color:var(--muted)}
.note{background:var(--card);border-left:3px solid var(--warn);padding:6px 10px;margin:6px 0}
figure{margin:12px 0}figure svg{max-width:100%;height:auto}a{color:var(--info)}
"""


def threshold_text(rule: R.Rule, variant: str | None) -> str:
    value = rule.extra["value_hmr"] if variant == "hmr" else rule.value
    return {"max": f"≤ {value}", "min": f"≥ {value}", "range": f"{value[0]} – {value[1]}" if rule.kind == "range" else "",
            "not_in": f"not {', '.join(value)}" if rule.kind == "not_in" else "", "bool": f"must be {value}",
            "manual": "manual check", "info": "value only"}[rule.kind]


def _json_safe(x):
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return None if math.isnan(x) else ("inf" if math.isinf(x) else x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    return x


def summary_counts(rows: list[Row]) -> dict[int, dict[str, int]]:
    counts = {g: {s: 0 for s in STATUS_ORDER} for g in GATES}
    for row in rows:
        counts[row.rule.gate][row.status] += 1
    return counts


def qc_dict(campaign: Campaign, ruleset: R.RuleSet) -> dict:
    rows = evaluate_results(ruleset, campaign.results)
    return {
        "traj_inspector_version": __version__, "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": ruleset.params, "notes": campaign.notes, "not_in_v1": NOT_IN_V1,
        "replicas": [{"name": r.name, "engine": r.engine, "setup": dataclasses.asdict(r.setup)} for r in campaign.replicas],
        "summary": {GATES[g]: c for g, c in summary_counts(rows).items()},
        "rows": [{"rule_id": r.rule.id, "gate": r.rule.gate, "check": r.rule.check, "replica": r.result.replica, "status": r.status,
                  "metric": _json_safe(r.result.metric), "detail": r.result.detail, "reason": r.result.reason,
                  "threshold": threshold_text(r.rule, r.result.variant), "severity": r.rule.severity, "source": r.rule.source,
                  "source_url": r.rule.source_url, "tested_on": r.rule.tested_on, "rule_of_thumb": r.rule.rule_of_thumb}
                 for r in rows],
    }


def export_qc_json(campaign: Campaign, ruleset: R.RuleSet) -> str:
    return json.dumps(qc_dict(campaign, ruleset), indent=2, default=_json_safe)


# ------------------------------------------------------------------ plots

def _svg(fig) -> str:
    buffer = io.StringIO()
    fig.savefig(buffer, format="svg", bbox_inches="tight")
    plt.close(fig)
    text = buffer.getvalue()
    return text[text.index("<svg"):]


def _q_plot(campaign: Campaign, params: dict) -> str:
    fig, ax = plt.subplots(figsize=(7, 2.6))
    for color, (name, s) in zip(PALETTE, campaign.series.items()):
        ax.plot(s.time_ps / 1000.0, s.q, lw=1, color=color, label=name)
    for key in ("q_unfolded", "q_folded"):
        ax.axhline(params[key], color="gray", ls="--", lw=0.8)
    ax.set(xlabel="time (ns)", ylabel="Q", ylim=(0, 1.02))
    ax.legend(fontsize=7, frameon=False, ncol=3)
    return _svg(fig)


def _fes_plot(campaign: Campaign) -> str:
    """Free energy (kT) of Q x Rg: pooled and per replica. A projection: read basins, not barriers."""
    groups = {"pooled": list(campaign.series.values()), **{n: [s] for n, s in campaign.series.items()}}
    fig, axes = plt.subplots(1, len(groups), figsize=(3.0 * len(groups), 2.8), squeeze=False)
    q_all = np.concatenate([s.q for s in campaign.series.values()])
    rg_all = np.concatenate([s.rg_nm for s in campaign.series.values()])
    ranges = [[min(q_all.min(), 0.0), max(q_all.max(), 1.0)], [rg_all.min() - 1e-3, rg_all.max() + 1e-3]]
    for ax, (title, group) in zip(axes[0], groups.items()):
        h, xe, ye = np.histogram2d(np.concatenate([s.q for s in group]), np.concatenate([s.rg_nm for s in group]), FES_BINS, ranges)
        with np.errstate(divide="ignore"):
            f = -np.log(h / h.sum())
        f -= f[np.isfinite(f)].min()
        image = ax.pcolormesh(xe, ye, np.where(np.isfinite(f), f, np.nan).T, cmap="viridis")
        ax.set(title=title, xlabel="Q")
        ax.title.set_fontsize(8)
    axes[0][0].set_ylabel("Rg (nm)")
    fig.colorbar(image, ax=axes[0].tolist(), label="F / kT", shrink=0.8)
    return _svg(fig)


def _its_plot(campaign: Campaign) -> str:
    fig, ax = plt.subplots(figsize=(4, 2.6))
    lags = sorted(l for l in campaign.its_ns if campaign.its_ns[l] is not None)
    ax.plot(lags, [campaign.its_ns[l] for l in lags], "o-")
    ax.set(xlabel="lag (frames)", ylabel="slowest ITS (ns)")
    return _svg(fig)


def _energy_plot(campaign: Campaign) -> str:
    fig, (a, b) = plt.subplots(1, 2, figsize=(7, 2.6))
    for color, rep in zip(PALETTE, campaign.replicas):
        e = rep.energies
        a.plot(e.time_ps / 1000.0, e.temperature, lw=1, color=color, label=rep.name)
        b.plot(e.time_ps / 1000.0, e.density_kg_m3, lw=1, color=color)
    a.set(xlabel="time (ns)", ylabel="T (K)")
    b.set(xlabel="time (ns)", ylabel="density (kg/m3)")
    a.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    return _svg(fig)


# ------------------------------------------------------------------ html

def _chip(status: str) -> str:
    return f'<span class="chip {status.replace("/", "")}">{status}</span>'


def _row_html(row: Row) -> str:
    res, rule = row.result, row.rule
    value = f"{res.detail} ({res.reason})" if res.detail and res.reason else (res.detail or res.reason)
    source = escape(rule.source) or "—"
    if rule.source_url:
        source = f'<a href="{escape(rule.source_url)}">{source}</a>'
    rot = ' <span class="rot">rule of thumb</span>' if rule.rule_of_thumb else ""
    return (f"<tr><td>{_chip(row.status)}</td><td>{escape(rule.check)}{rot}</td><td>{escape(res.replica)}</td>"
            f'<td class="num">{escape(value)}</td><td>{escape(threshold_text(rule, res.variant))}</td>'
            f"<td>{source}</td><td>{escape(rule.tested_on) or chr(0x2014)}</td></tr>")


def _gate_html(gate: int, rows: list[Row], order: dict[str, int]) -> str:
    head = "<tr><th>Status</th><th>Check</th><th>Replica</th><th>Value</th><th>Threshold</th><th>Source</th><th>Tested on</th></tr>"
    body = "".join(_row_html(r) for r in sorted(rows, key=lambda r: (order[r.rule.id], r.result.replica == CAMPAIGN)) if r.rule.gate == gate)
    return f"<h2>{GATES[gate]}</h2><table>{head}{body}</table>"


def render_html(campaign: Campaign, ruleset: R.RuleSet) -> str:
    """Complete report. Statuses come from the rules passed in, so editing a rule only needs a new call."""
    rows = evaluate_results(ruleset, campaign.results)
    counts = summary_counts(rows)
    order = {r.id: i for i, r in enumerate(ruleset.rules)}
    cards = "".join(
        f'<div class="card"><b>{escape(GATES[g])}</b>' + "".join(f'<span>{_chip(s)} {c[s]}</span>' for s in (R.PASS, R.WARN, R.FAIL) )
        + f'<div class="muted">N/A {c[R.NA]} · manual {c[R.MANUAL]} · info {c[R.INFO]}</div></div>' for g, c in counts.items())
    notes = "".join(f'<div class="note">{escape(n)}</div>' for n in campaign.notes)
    replicas = ", ".join(f"{escape(r.name)} ({r.engine})" for r in campaign.replicas)
    figures = [("Temperature and density", _energy_plot(campaign)), ("Q(t) per replica", _q_plot(campaign, ruleset.params)),
               ("Free energy of Q x Rg (projection: no barriers)", _fes_plot(campaign))]
    if campaign.its_ns:
        figures.append(("Implied timescale vs lag", _its_plot(campaign)))
    figs = "".join(f"<figure><figcaption class='muted'>{escape(t)}</figcaption>{svg}</figure>" for t, svg in figures)
    return (f"<!doctype html><html lang='en'><meta charset='utf-8'><title>Traj Inspector report</title><style>{CSS}</style><body>"
            f"<h1>Trajectory evaluation</h1><p class='muted'>Replicas: {replicas}. Traj Inspector {__version__}.</p>{notes}"
            f"<div class='cards'>{cards}</div>" + "".join(_gate_html(g, rows, order) for g in GATES)
            + f"<h2>Plots</h2>{figs}<h2>Not in v1</h2><p class='muted'>{escape(NOT_IN_V1)}.</p></body></html>")
