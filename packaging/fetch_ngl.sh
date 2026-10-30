#!/usr/bin/env bash
# Fetch the pinned NGL build from the npm registry into trajinspector/viewer/. Build time only: the app never fetches.
set -euo pipefail
VERSION=2.5.0
TARBALL_SRI=sha512-7pvmhCxBTBb9dxMVrqDMR4Lb2ijrXFckokM0z16vZ0pO36x+j5wV8fGJJQSAK4hpQEzpawk+BHu0bJAFSkXOCA==
DEST="$(cd "$(dirname "$0")/.." && pwd)/trajinspector/viewer"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
curl -fsSL "https://registry.npmjs.org/ngl/-/ngl-${VERSION}.tgz" -o "$WORK/ngl.tgz"
ACTUAL="sha512-$(openssl dgst -sha512 -binary "$WORK/ngl.tgz" | base64 -w0)"
[ "$ACTUAL" = "$TARBALL_SRI" ] || { echo "ngl-${VERSION}.tgz checksum mismatch: $ACTUAL" >&2; exit 1; }
tar xzf "$WORK/ngl.tgz" -C "$WORK"
cp "$WORK/package/dist/ngl.js" "$DEST/ngl.js"
cp "$WORK/package/LICENSE" "$DEST/NGL-LICENSE"
echo "ngl ${VERSION} -> $DEST/ngl.js ($(sha256sum "$DEST/ngl.js" | cut -d' ' -f1))"
