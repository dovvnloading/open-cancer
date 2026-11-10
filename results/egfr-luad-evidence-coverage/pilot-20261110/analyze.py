#!/usr/bin/env python3
"""Run the pinned Open Cancer EGFR/LUAD pilot from the repository root."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from scripts.egfr_luad_pilot import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
