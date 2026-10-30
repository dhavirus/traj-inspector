# Traj Inspector — CLAUDE.md

Offline desktop app for Ubuntu that evaluates MD trajectories against a five-gate checklist
(setup → artifacts → sampling → mechanism → experiment). Gates 1–4 run automatically; gate 5 is not
in scope for v1. Every threshold is a rule in an editable `rules.yaml` with its source paper, the
system that paper tested on, and a `rule_of_thumb` flag when no paper sets the number.

Owner: Bem (computational chemist, CADD). Reply in English; direct, technical style.
Review style: summary + evidence of what was verified and how.

## Decisions already made (do not re-ask)

| # | Decision |
|---|---|
| 1 | Input engines: **GROMACS, Amber and OpenMM** |
| 2 | Scope v1: **gates 1–4** (gate 5 experiment comparison = later) |
| 3 | Target: **Ubuntu 22.04 PC with no internet** → nothing fetched at runtime (no CDN, no fonts online, no update check) |
| 4 | **Rules editable inside the app** (settings table → writes user `rules.yaml`; report re-evaluates without recomputing) |
| — | Packaging: one **AppImage**, plus a `.tar.gz` of the AppDir as fallback |
| — | Stack: Python analysis core + **PySide6** window + **QtWebEngine** for report and 3D viewer (**NGL**, bundled locally) |

Rejected: Electron (second runtime, bigger), browser-only app (bad with large local files, not
self-contained), calling `gmx` (breaks self-contained; MDAnalysis + pyedr read tpr/xtc/edr), PyMOL/VMD
embedding (licence, size).

## Architecture

```
trajinspector/
  __main__.py      python -m trajinspector → launches the Qt window
  rules.yaml       shipped defaults (copied to ~/.config/trajinspector/rules.yaml on first edit)
  rules.py         load defaults + user override, save, evaluate one rule on a metric
  readers.py       engine detection by file extension, normalized Setup + Energies per replica
  checks.py        metrics for gates 1–4 → results table (rule × replica, plus campaign rows)
  report.py        self-contained HTML report (inline SVG plots from matplotlib), qc.json export
  viewer/          viewer.html + ngl.js (pinned, bundled) – multi-model PDB trajectory player
  app.py           Qt UI
tests/             pytest on real GROMACS/Amber/OpenMM test data + synthetic Q series
packaging/         build_appdir.sh, AppRun, glibc check, Colab build notebook
```

UI: left = replica list ("Add replica folder…", auto-classify files), optional reference structure,
atom selection (default `protein`), Run button with progress (QThread). Right tabs: **Report**
(gate cards with pass/warn/fail counts, rows with status chip, value, threshold, source, tested-on,
rule-of-thumb label), **Viewer** (replica picker, slider, play, per-frame time + Q + state badge),
**Rules** (table: enabled, value, severity editable; Save, Reset to defaults; Apply re-renders report).
Menu: Export HTML report, Export qc.json.

Statuses: PASS, WARN, FAIL, N/A (with reason), MANUAL (needs a human), INFO (value only).

## Normalized data model (engine-agnostic)

`Setup`: thermostat, barostat, dt_ps, constraints, cutoff_nm, hmr_factor, warnings_count.
`Energies`: time_ps, temperature, kinetic, potential, conserved (may be absent), density_kg_m3, volume.

### Reader facts (verified on real files)

**GROMACS** (`.tpr` + `.xtc/.trr` parts + `.edr` parts + `.log` parts + `md.mdp`/`mdout.mdp`)
- MDAnalysis reads the tpr (masses, bonds). Settings come from the mdp (tpr doesn't expose them).
- `-noappend` parts (`md.partNNNN.*`): **do not use `continuous=True`** — ChainReader fails when a part
  has <2 frames. Load with plain chaining and drop frames whose time ≤ previous time. Skip empty part files.
- pyedr keys: `Time, Potential, Kinetic En., Total Energy, Conserved En., Temperature, Pressure, Volume, Density`.
  Concatenate parts, dedupe on time.
- `integrator = sd` → thermostat `langevin`. `mass-repartition-factor` in the mdp (also visible in tpr H masses).
- Warnings: count `LINCS WARNING` in logs.

**Amber** (`.prmtop/.parm7` + `.nc/.mdcrd/.dcd` + `mdout`)
- Settings: parse `key = value` pairs from the echoed input ("Here is the input file:" section) of mdout,
  lower-cased. `ntt`: 1 berendsen, 2 andersen, 3 langevin. `ntp>0` = pressure coupling; barostat:
  text `Monte-Carlo Barostat` in mdout → monte-carlo, else berendsen (barostat=1 is the default).
  `ntc` 1 none / 2 h-bonds / 3 all-bonds. `dt` in ps. `cut` in Å → nm.
- Energies: blocks start `NSTEP = … TIME(PS) = … TEMP(K) = …`, then `Etot`, `EKtot`, `EPtot`, …,
  `VOLUME`, `Density` (g/cc). kcal/mol × 4.184. Stop at `A V E R A G E S`.
- Langevin (`ntt=3`) has no conserved energy → drift check = N/A with reason.
- Warnings: count `vlimit exceeded`.

**OpenMM** (`.pdb/.cif` (or prmtop) + `.dcd/.xtc` + StateDataReporter `.csv` + `system.xml` + `integrator.xml`)
- CSV header: `#"Step","Time (ps)","Potential Energy (kJ/mole)","Kinetic Energy (kJ/mole)","Total Energy (kJ/mole)","Temperature (K)","Box Volume (nm^3)","Density (g/mL)"`.
- integrator.xml: `<Integrator type="LangevinMiddleIntegrator" stepSize=".004" temperature=… friction=…>`.
- system.xml: `<Force type="MonteCarloBarostat">`, `<Force type="NonbondedForce" cutoff=".9">`, count `<Constraint `,
  `<Particle mass=…>` (use these masses — a PDB topology has none, MDAnalysis guesses H = 1.008).
- PDB topology has no bonds → `guess_bonds()` on the selection before unwrapping.
- No conserved energy for Langevin → drift N/A. Warnings: N/A.

**HMR detection (all engines):** mean mass of hydrogens in the selection / 1.008. Water H stays 1.008
(rigid water), so only look at the solute.

## Rules (shipped defaults)

Kinds: `max`, `min`, `range`, `not_in` (list), `bool`, `manual`, `info`. Fields: `id, gate, check, kind,
value, severity (fail|warn), source, source_url, tested_on, rule_of_thumb, enabled, note`.
Parameters block: `q_unfolded 0.3`, `q_folded 0.8`, `max_frames 5000`, `msm_lag_ns`, `image_stride`.

| Gate | id | Pass rule | Source → tested on | RoT |
|---|---|---|---|---|
| 1 | thermostat | not in [berendsen] | Bussi 2007 JCP; GROMACS manual → LJ argon, TIP4P water | |
| 1 | barostat | not in [berendsen]; warn Parrinello-Rahman from unequilibrated box | Bernetti & Bussi 2020 JCP → LJ fluid, water, membrane | |
| 1 | forcefield_water | MANUAL: water model matches force field | Robustelli 2018 PNAS (21 proteins); chignolin FF eval ChemRxiv 2024 | |
| 1 | dt_max_fs | ≤ 2 fs (≤ 4 fs with HMR) | GROMACS manual; Jung 2021 JCTC | |
| 1 | hmr_rates | WARN if HMR: no rates from this run | Sahil 2023 Biophys J (T4L, MopR, galectin-3); Bowers 2026 JCTC (Ala dipeptide) | |
| 1 | constraints | h-bonds when dt > 1 fs | GROMACS manual | |
| 1 | cutoff_min_nm | ≥ 0.9 nm | — | ✓ |
| 1 | box_padding_nm | frame-0 min image distance / 2 ≥ 1.2 nm | — | ✓ |
| 2 | drift | \|slope\| < 1e-3 kJ/mol/ps/atom over ≥ 1 ns (needs conserved energy) | — | ✓ |
| 2 | ke_width | sd(T)/mean(T) ÷ √(2/ndof) in 0.85–1.15, ≥ 200 frames; ndof = 2·KE/(kB·T) | Merz & Shirts 2018 PLoS ONE | ✓ (band) |
| 2 | warnings | constraint warnings = 0 | GROMACS manual | |
| 2 | min_image | min protein–own-image distance > cutoff, all analysed frames | — | ✓ |
| 2 | density_drift | \|slope\| ≤ 1 kg/m³/ns after equilibration | — | ✓ |
| 2 | rmsd_jump | frame-to-frame Cα RMSD < 0.5 nm | — | ✓ |
| 2 | group_temps | N/A in v1 (needs velocities) | Wong-ekkabut & Karttunen 2016 | |
| 3 | equilibration | INFO: t0 from pymbar `detect_equilibration` | Chodera 2016 JCTC → liquid argon, 500 runs | |
| 3 | n_eff_min | N_eff(Q) ≥ 20 | Grossfield 2018 LiveCoMS | ✓ (number) |
| 3 | block_plateau | SE(largest block, ≥ 8 blocks) / SE(half that) ≤ 1.25 | Grossfield 2018 | ✓ (number) |
| 3 | replicas_min | ≥ 3 replicas | Grossfield 2018 | ✓ (number) |
| 3 | replica_agreement | all pairs \|ff_i − ff_j\| ≤ 2·√(se_i²+se_j²), se = √(p(1−p)/N_eff) | Grossfield 2018 | ✓ |
| 3 | transitions | ≥ 5 (fail), ≥ 20 for Tm (warn), ≥ 50 for rates (warn); two-threshold Q states | — | ✓ |
| 4 | q_bimodal | fraction of frames with q_unfolded < Q < q_folded ≤ 0.2 | Best, Hummer & Eaton 2013 PNAS → Shaw fast folders | ✓ (number) |
| 4 | vamp2 | INFO: VAMP-2 score of features | Scherer et al. arXiv 1811.11714 → 12 fast folders incl. CLN025 | |
| 4 | its_flat | slowest ITS changes ≤ 15 % between lag and 2·lag | Prinz 2011 JCP → MR121-GSGS-W peptide | ✓ (number) |
| 4 | ck_test | max CK error ≤ 0.05 (deeptime `MarkovStateModel.ck_test`, 2 metastable sets) | Prinz 2011 | ✓ (number) |
| 4 | fes | plot only (Q × Rg, pooled + per replica); note: no barriers from projections | Grossfield 2018 | |

Source URLs and full "tested on" text: see the published explainer *Trajectory Evaluation Gates*
(claude.ai artifact) and the review doc *Modern MD Practice*. Copy them into `rules.yaml`.

Not in v1 (say so in the report): DSSP per state, committor test, force distribution analysis,
group temperatures, gate 5.

## Analysis details

- Q: Best–Hummer–Eaton. Heavy-atom pairs < 4.5 Å in the reference, |i−j| > 3 residues, β = 5 Å⁻¹, λ = 1.8.
  Reference = user file, else frame 0 of replica 1 (WARN). Match atoms by (residue index, atom name).
- Unwrap the selection before RMSD/Rg/Q (`transformations.unwrap`; guess bonds if none).
- Min image distance: unwrapped selection vs itself shifted by the 26 lattice vectors (triclinic,
  `triclinic_vectors`), `capped_distance`, every `image_stride` frames.
- Transitions: two thresholds; frames between keep the previous state (port `count_transitions` from
  `cln025.py`, it has a self-test).
- TICA/MSM: Cα–Cα distances (i, j ≥ i+2); if > 50 Cα, subsample Cα to ≤ 50 (mark as `shortcut:`).
  ITS at lag × [0.25, 0.5, 1, 2, 4]; catch `ValueError` for lags with too little data.
- Viewer: write selection-only multi-model PDB per replica (≤ 500 frames) + JSON (time, Q per frame)
  into a temp dir with `viewer.html` + `ngl.js`; QWebEngine with `LocalContentCanAccessFileUrls` on;
  NGL `stage.loadFile("traj.pdb", {asTrajectory: true})`.

## Packaging (Ubuntu 22.04, offline)

- Target glibc ≤ **2.35**. Build host may be newer → install wheels with
  `pip install --target … --only-binary :all: --platform manylinux_2_35_x86_64 --python-version 3.12`.
  Then scan every `.so` for the highest `GLIBC_x.y` symbol and **fail the build if > 2.35**.
- Python runtime: python-build-standalone (glibc 2.17) or conda-forge `python=3.12` via micromamba.
- Bundle `libxcb-cursor.so.0` (Qt ≥ 6.5 xcb plugin needs it; not installed on stock 22.04).
- AppRun sets `QTWEBENGINE_DISABLE_SANDBOX=1` (sandbox can't be setuid inside an AppImage mount) and runs
  `python -m trajinspector`.
- Use the **static type-2 runtime** (AppImage/appimagetool continuous) so libfuse2 is not needed.
  Document `--appimage-extract-and-run` as the fallback. Also ship `TrajInspector-x86_64.tar.gz` (AppDir).
- Pin and checksum `ngl.js` (npm `ngl`, fetched at build time only).
- Expected size 300–450 MB (estimate).
- **Preferred build route: GitHub Actions on an `ubuntu-22.04` runner** (`.github/workflows/appimage.yml`):
  native 22.04 build, appimagetool downloadable there, AppImage attached to a GitHub release on tag push.
  Bem downloads it once and carries it to the offline PC.
- Fallback: in the old sandbox, GitHub release downloads were blocked (pypi, npm, conda-forge worked). If
  Actions is not available, the final AppImage step runs in Colab via a build notebook that writes to Drive.
- Reference code: `reference/cln025.py` (copied from Drive) — port logic from it, do not import it.

## Test data and verification

Build tiny test systems for all three engines (CLN025, `YYDPETGTWY`):
- GROMACS: CPU `gmx` (apt `gromacs`), a99SB-disp (github.com/paulrobustelli/Force-Fields @ 25e729d) or any
  ff; a few ps with `-noappend` chunks (include a 1-frame part to test the ChainReader fix).
- Amber: tleap ff14SB + TIP3P, 8 Å box; rename C-terminal `OC1/OC2` → `O/OXT`; sander 6 ps,
  ntt=3, ntp=1, barostat=2, ntc=2, dt=0.002, ntwx=ntpr=25, ioutfm=1.
- OpenMM: same prmtop, LangevinMiddle 4 fs + `hydrogenMass=3 amu` + MonteCarloBarostat; DCD +
  StateDataReporter CSV + serialized system/integrator XML. The app must flag HMR here.

Done = all of:
1. pytest passes on the three engines + synthetic Q series (transitions, bimodality, replica agreement).
2. App starts headless (`QT_QPA_PLATFORM=offscreen`), runs an analysis, `widget.grab()` screenshots of
   Report, Viewer and Rules tabs look right (look at them).
3. Editing a rule value in the app changes the report status without recomputing.
4. AppDir runs via `./AppRun` on the build host; glibc scan ≤ 2.35; no network access at runtime
   (run once with networking blocked).
5. Reply to Bem with summary + evidence; say what could not be tested (real desktop GUI, the 22.04 machine).

## Conventions from Bem

- If anything goes to Google Drive: code under `googleColab_run/9.traj-inspector/` (`src/`, `notebooks/`);
  runs in `<yymm>/run_<yymmdd>_<hhmmss>/` with `config.json`; new versions as `name_v2`, old ones to `_archive/`.
- Notebooks: real line breaks per cell, `nbformat.validate()` before delivery; config dict at top;
  `importlib.reload`, never `%autoreload`.
- Big picture first; name rejected alternatives; evidence over claims.

## Status at hand-off (2026-10-30)

- Decisions 1–4 confirmed. Design above agreed. **No app code written yet.**
- Verified so far: MDAnalysis 2.11 reads the GROMACS tpr + xtc parts, Amber prmtop + nc, OpenMM pdb + dcd;
  pyedr keys listed above; deeptime 0.4.5 has `MarkovStateModel.ck_test`; Amber mdout and OpenMM XML formats above.
- Test data from the old sandbox is not carried over — regenerate it with the recipe above.
- Related existing code: `cln025.py` (Drive `googleColab_run/7.protein-folding/src/`) has
  `count_transitions`, `best_hummer_q`, `stat_inefficiency`, `tica_msm`, QC helpers — reuse logic.

## Status after the v1 build (2026-10-30)

Built: rules.yaml (28 rules) + rules.py, readers.py, analysis.py, checks.py, report.py, viewer/, app.py, packaging/, appimage.yml.
Tests: `pytest` (69) incl. offscreen window test; `tests/make_screenshots.py OUT` for the three tab screenshots;
`tests/make_testdata.py` regenerates `tests/data` (needs gmx, AmberTools, openmm).

Corrections to the facts above, found while building:
- `cln025.py` was at the repo root; it now lives in `reference/`.
- `pip install --platform manylinux_2_35_x86_64` alone matches nothing (pip 26.2 does not expand to older manylinux tags);
  `build_appdir.sh` lists every tag up to 2_35.
- pymbar `detect_equilibration` cuts almost all frames of a two-state Q series (3835 of 4000 in a synthetic telegraph series).
  N_eff, block SE and transitions therefore use all frames; t0 is INFO only (see `analysis.stat_inefficiency`).
- Test data: GROMACS with amber99sb-ildn + TIP3P (apt gromacs 2023.3; GitHub downloads of a99SB-disp were blocked), so
  `mass-repartition-factor` (GROMACS >= 2024) was not exercised; HMR is detected from masses.
  In a `-noappend` restart, `-nsteps` counts additional steps; GROMACS repeats the checkpoint frame in the next part.
- Q needs a compact reference: an extended chain has no native contacts and `run_checks` stops with a ValueError.
- Citations in rules.yaml: only Bussi 2007, Bernetti & Bussi 2020 and Sahil 2023 were confirmed in PubMed (source_url set);
  the others are as given above and marked "not verified". The explainer and review doc were not available in this session.
