"""Fail if any ELF file under a directory needs a GLIBC symbol version above the limit (default 2.35, Ubuntu 22.04).

usage: python glibc_check.py DIR [MAX_VERSION]
"""
import re
import subprocess
import sys
from pathlib import Path

ELF_MAGIC = b"\x7fELF"
VERSION = re.compile(r"Name: GLIBC_(\d+)\.(\d+)(?:\.(\d+))?")


def needed_glibc(path: Path) -> tuple[int, ...]:
    """Highest GLIBC_x.y[.z] version a file requires, from its version-needs table; () if none."""
    out = subprocess.run(["readelf", "-V", str(path)], capture_output=True, text=True).stdout
    return max((tuple(int(g) for g in m.groups() if g is not None) for m in VERSION.finditer(out)), default=())


def main(root: Path, limit: tuple[int, ...]) -> int:
    worst, offenders, scanned = (), [], 0
    for path in sorted(root.rglob("*")):
        if path.is_symlink() or not path.is_file() or path.stat().st_size < 64:
            continue
        with open(path, "rb") as f:
            if f.read(4) != ELF_MAGIC:
                continue
        scanned += 1
        need = needed_glibc(path)
        worst = max(worst, need)
        if need > limit:
            offenders.append((need, path))
    print(f"scanned {scanned} ELF files under {root}; highest GLIBC needed: {'.'.join(map(str, worst)) or 'none'} (limit {'.'.join(map(str, limit))})")
    for need, path in sorted(offenders, reverse=True)[:25]:
        print(f"  TOO NEW GLIBC_{'.'.join(map(str, need))}: {path.relative_to(root)}")
    return 1 if offenders else 0


if __name__ == "__main__":
    limit = tuple(int(p) for p in (sys.argv[2] if len(sys.argv) > 2 else "2.35").split("."))
    sys.exit(main(Path(sys.argv[1]), limit))
