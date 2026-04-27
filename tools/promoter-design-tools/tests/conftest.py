from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def run_cli(args: list[str], cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    completed = subprocess.run([sys.executable, *args], cwd=cwd or ROOT, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    return completed


def example_fasta() -> Path:
    return ROOT / "tools" / "promoter_window_normalizer" / "test-data" / "example_promoters.fasta"

