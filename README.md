# Traj Inspector

Offline desktop app (Ubuntu 22.04) that evaluates MD trajectories (GROMACS, Amber, OpenMM) against gates 1-4 of a
five-gate checklist. Thresholds live in `trajinspector/rules.yaml` and can be edited in the Rules tab.
See CLAUDE.md for the design.

## Run from source

```bash
pip install -e ".[test]"
packaging/fetch_ngl.sh          # bundles NGL (build time only; the app never fetches anything)
python -m trajinspector
```

## AppImage

```bash
chmod +x TrajInspector-x86_64.AppImage
./TrajInspector-x86_64.AppImage
```

The static type-2 runtime needs no libfuse2. If FUSE is unavailable: `./TrajInspector-x86_64.AppImage --appimage-extract-and-run`.
Fallback: unpack `TrajInspector-x86_64.tar.gz` and run `AppDir/AppRun`.
Expected from the host: the usual desktop GL/X11 stack (libEGL, libGL, libX11 family, libdbus, fontconfig, freetype, glib).
Everything else, including libxcb-cursor, is inside the AppDir.

## Replica folder contents

| Engine | Files |
|---|---|
| GROMACS | `*.tpr`, `*.xtc`/`*.trr` parts, `*.edr` parts, `*.log` parts, `mdout.mdp` or `md.mdp` |
| Amber | `*.prmtop`/`*.parm7`, `*.nc`/`*.mdcrd`, `mdout*` (MD runs only; minimization outputs are skipped) |
| OpenMM | `topology.pdb`/`*.cif`/prmtop, `*.dcd`/`*.xtc`, `state*.csv`, `system.xml`, `integrator.xml` |
