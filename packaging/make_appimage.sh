#!/usr/bin/env bash
# Pack build/AppDir into TrajInspector-x86_64.AppImage with the static type-2 runtime (no libfuse2 needed to run it).
# The build host runs appimagetool without FUSE (APPIMAGE_EXTRACT_AND_RUN). Fallback for users: ./TrajInspector-x86_64.AppImage --appimage-extract-and-run
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${OUT:-$ROOT/build}"
BASE=https://github.com/AppImage
curl -fsSL -o "$OUT/appimagetool" "$BASE/appimagetool/releases/download/continuous/appimagetool-x86_64.AppImage"
curl -fsSL -o "$OUT/runtime-x86_64" "$BASE/type2-runtime/releases/download/continuous/runtime-x86_64"
chmod +x "$OUT/appimagetool"
sha256sum "$OUT/appimagetool" "$OUT/runtime-x86_64"   # 'continuous' moves; the build log is the record of what was used
ARCH=x86_64 APPIMAGE_EXTRACT_AND_RUN=1 "$OUT/appimagetool" --runtime-file "$OUT/runtime-x86_64" "$OUT/AppDir" "$OUT/TrajInspector-x86_64.AppImage"
ls -la "$OUT/TrajInspector-x86_64.AppImage"
