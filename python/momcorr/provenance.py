"""Small provenance utilities shared by the optional Python workflows."""

import hashlib
from pathlib import Path
import subprocess


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def code_state():
    root = Path(__file__).resolve().parents[2]
    return {
        "commit": subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"], text=True
        ).strip(),
        "working_tree_status": subprocess.check_output(
            ["git", "-C", str(root), "status", "--porcelain"], text=True
        ),
        "modules": {
            p.name: sha256_file(p)
            for p in sorted((root / "python/momcorr").glob("*.py"))
        },
    }
