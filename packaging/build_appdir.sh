#!/usr/bin/env bash
# Build build/AppDir (Python 3.12 from conda-forge, glibc 2.17 sysroot) + TrajInspector-x86_64.tar.gz.
# Works on a host newer than Ubuntu 22.04 because wheels are chosen for manylinux tags <= 2_35 and the result is glibc-scanned.
# usage: packaging/build_appdir.sh      (OUT=dir to change the output directory; it must not exist yet)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${OUT:-$ROOT/build}"
APPDIR="$OUT/AppDir"
MICROMAMBA_BUILD=2.10.0-0
PYTHON_VERSION=3.12
GLIBC_MAX=2.35
# Libraries Qt's xcb plugin and QtWebEngine want that a stock Ubuntu 22.04 may lack (libxcb-cursor0 is the known one).
CONDA_LIBS="xcb-util-cursor xcb-util-image xcb-util-keysyms xcb-util-renderutil xcb-util-wm xcb-util libxkbcommon nss nspr alsa-lib"
# Qt modules this app never imports; pruned to save space. Verified afterwards by --self-check.
PRUNE_QT="Qt3D Charts DataVisualization Graphs Multimedia Quick3D Location Sensors Bluetooth Nfc SerialPort SerialBus RemoteObjects Scxml SpatialAudio TextToSpeech HttpServer Designer UiTools Help"

[ ! -e "$OUT" ] || { echo "$OUT exists; remove it first" >&2; exit 1; }
[ -f "$ROOT/trajinspector/viewer/ngl.js" ] || "$ROOT/packaging/fetch_ngl.sh"
mkdir -p "$OUT/mm"

curl -fsSL "https://conda.anaconda.org/conda-forge/linux-64/micromamba-${MICROMAMBA_BUILD}.tar.bz2" | tar xj -C "$OUT/mm" bin/micromamba
MAMBA_ROOT_PREFIX="$OUT/mamba" "$OUT/mm/bin/micromamba" create -y -q -p "$APPDIR/usr" -c conda-forge --override-channels \
    "python=${PYTHON_VERSION}" $CONDA_LIBS

PY="$APPDIR/usr/bin/python${PYTHON_VERSION}"
SITE="$APPDIR/usr/lib/python${PYTHON_VERSION}/site-packages"
# pip does not expand --platform manylinux_2_35_x86_64 to older manylinux tags, so list every tag up to glibc 2.35 explicitly.
PLATFORMS=()
for tag in manylinux_2_35 manylinux_2_34 manylinux_2_31 manylinux_2_28 manylinux_2_27 manylinux_2_24 manylinux_2_17 manylinux2014 manylinux_2_12 manylinux2010 manylinux_2_5 manylinux1; do
    PLATFORMS+=(--platform "${tag}_x86_64")
done
"$PY" -m pip install --no-compile --only-binary :all: --target "$SITE" --python-version "$PYTHON_VERSION" \
    --implementation cp --abi "cp${PYTHON_VERSION/./}" "${PLATFORMS[@]}" -r "$ROOT/packaging/requirements.txt"
cp -r "$ROOT/trajinspector" "$SITE/"
mkdir -p "$APPDIR/usr/share/trajinspector"
"$PY" -m pip list --path "$SITE" --format=freeze > "$APPDIR/usr/share/trajinspector/requirements-frozen.txt"

# prune: unused Qt modules, translations, dev files, caches, test suites
for m in $PRUNE_QT; do
    rm -rf "$SITE"/PySide6/Qt*"$m"*.abi3.so "$SITE"/PySide6/Qt/lib/libQt6*"$m"*.so* "$SITE"/PySide6/Qt/qml/Qt*"$m"* "$SITE"/PySide6/Qt/plugins/*/*"$m"* 2>/dev/null || true
done
rm -rf "$SITE"/PySide6/{examples,include,typesystems,glue,scripts} "$SITE"/PySide6/Qt/metatypes \
       "$SITE"/PySide6/Qt/libexec/{assistant,designer,linguist,lupdate,lrelease,qmlformat,qmllint,qmlls} "$APPDIR"/usr/include "$APPDIR"/usr/share/man
find "$SITE/PySide6/Qt/translations" -maxdepth 1 -name "*.qm" -delete   # keep qtwebengine_locales: Chromium needs its .pak files
find "$APPDIR/usr" -name __pycache__ -type d -prune -exec rm -rf {} +
find "$SITE" -type d \( -name tests -o -name test \) -prune -exec rm -rf {} +

chmod -R a+rX "$APPDIR"   # some wheel directories are extracted 0700; appimagetool refuses those
install -m 755 "$ROOT/packaging/AppRun" "$APPDIR/AppRun"
install -m 644 "$ROOT/packaging/trajinspector.desktop" "$APPDIR/trajinspector.desktop"
"$PY" "$ROOT/packaging/make_icon.py" "$APPDIR/trajinspector.png"
cp "$APPDIR/trajinspector.png" "$APPDIR/.DirIcon"

echo "== libraries still expected from the host (GL/EGL/X11 stack)"
for f in "$SITE"/PySide6/Qt/libexec/QtWebEngineProcess "$SITE"/PySide6/Qt/lib/libQt6WebEngineCore.so.6 "$SITE"/PySide6/Qt/plugins/platforms/libqxcb.so; do
    LD_LIBRARY_PATH="$APPDIR/usr/lib" ldd "$f" | awk '/not found/ {print $1}'
done | sort -u
echo "== glibc scan (limit $GLIBC_MAX)"
"$PY" "$ROOT/packaging/glibc_check.py" "$APPDIR" "$GLIBC_MAX"

tar -C "$OUT" -czf "$OUT/TrajInspector-x86_64.tar.gz" AppDir
du -sh "$APPDIR" "$OUT/TrajInspector-x86_64.tar.gz"
