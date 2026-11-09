#!/usr/bin/env python3
"""Check (or fetch) the SOURCES 4A code and libraries used for the 2-Flux calculations.

The code and libraries are not distributed with this repository. Obtain SOURCES from the Oak Ridge RSICC package and
the modified libraries from the University of Sheffield group (see 2-Flux/SOURCES/README.md). This script checks that a
directory holds exactly the files used here.

  python scripts/get_sources.py                 # check $SOURCES_DIR or 2-Flux/SOURCES/code
  python scripts/get_sources.py --dir PATH      # check PATH
  python scripts/get_sources.py --clone         # clone the private copy (maintainers only) and check it
"""
import argparse
import hashlib
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_DIR = REPO / "2-Flux" / "SOURCES" / "code"
PRIVATE_URL = "https://github.com/villano-lab/sources4a-local.git"

# SHA-256 of the files used for the results in this repository (relative to the SOURCES directory).
EXPECTED = {
    "src/sources_v20240312.f": "0c2881202031e8e9e57ebd55c1180a6a87e0042eebd5c228a100b14eae5576bd",
    "libs/tape2": "a228bc33e34465cca0f28ad45506adcd1a4304cbec719f1974413e5a9c1c823f",
    "libs/tape3": "4cfdb1fe900c886d0a51f571c21a7a055509b249e40d6e84f68cc6cfed39cb43",
    "libs/tape4": "cead67bcb335108cb3443f757cc0584b12cd412238d7382e9b2aa3101fbef261",
    "libs/tape5": "a54d9c08697392fab55218f1a8ef6a5fd249268245284d78dbebd99f9fea73d4",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check(directory):
    ok = True
    for rel, want in EXPECTED.items():
        p = directory / rel
        if not p.is_file():
            print(f"missing   {rel}")
            ok = False
        elif sha256(p) != want:
            print(f"different {rel} (checksum does not match the version used here)")
            ok = False
        else:
            print(f"ok        {rel}")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", help="directory holding src/ and libs/")
    ap.add_argument("--clone", action="store_true", help="clone the private copy into the default directory")
    args = ap.parse_args()

    directory = Path(args.dir or os.environ.get("SOURCES_DIR") or DEFAULT_DIR)
    if args.clone:
        if directory.exists() and any(directory.iterdir()):
            sys.exit(f"{directory} exists and is not empty; not cloning")
        if subprocess.call(["git", "clone", PRIVATE_URL, str(directory)]) != 0:
            sys.exit("clone failed: you need access to the private repository (see 2-Flux/SOURCES/README.md)")
    if not directory.is_dir():
        sys.exit(f"{directory} not found. See 2-Flux/SOURCES/README.md for how to obtain SOURCES.")
    print(f"checking {directory}")
    sys.exit(0 if check(directory) else 1)


if __name__ == "__main__":
    main()
