"""Build aggregate paper tables from per-well summaries."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Sequence

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from cemdisp.data.loaders.paper import iter_paper_wells  # noqa: E402
from scripts.paper_data.run_paper_main_results import RESULTS_ROOT, write_all_wells_main_results  # noqa: E402


def build_paper_tables(root: Path = RESULTS_ROOT) -> list[dict[str, object]]:
    """Rebuild paper main result tables from per-well summaries."""

    well_ids = [record.well_id for record in iter_paper_wells(include_pending=True)]
    return write_all_wells_main_results(well_ids, root=root)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Build aggregate paper tables")
    parser.add_argument("--output-root", type=Path, default=RESULTS_ROOT)
    args = parser.parse_args(argv)
    rows = build_paper_tables(args.output_root)
    print(f"Wrote main result rows: {len(rows)}")


if __name__ == "__main__":
    main()
