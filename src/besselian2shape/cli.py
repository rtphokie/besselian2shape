"""Command-line interface for besselian2shape."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .cache import download_besselian_csv
from .elements import find_by_date
from .export import generate_eclipse_kml, generate_eclipse_shapefiles

_ALL_FORMATS = ("shp", "kml", "kmz")


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="besselian2shape",
        description="Generate ESRI shapefiles and/or KML/KMZ files for a solar eclipse's path.",
    )
    parser.add_argument(
        "year",
        type=int,
        help="eclipse year, astronomical numbering (1 BCE = 0, 2 BCE = -1, ...)",
    )
    parser.add_argument("month", type=int, help="eclipse month (1-12)")
    parser.add_argument("day", type=int, help="eclipse day of month")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        metavar="DIR",
        help="output directory (default: ./eclipse_YYYY-MM-DD)",
    )
    parser.add_argument(
        "-f",
        "--format",
        dest="formats",
        action="append",
        choices=(*_ALL_FORMATS, "all"),
        metavar="{shp,kml,kmz,all}",
        help="output format to generate; may be repeated (default: shp)",
    )
    parser.add_argument(
        "--step-minutes",
        type=float,
        default=0.5,
        metavar="MINUTES",
        help="time resolution used to sample the eclipse path (default: 0.5)",
    )
    parser.add_argument(
        "--penumbral-resolution",
        type=float,
        default=0.25,
        metavar="DEGREES",
        help="grid resolution for the penumbral visibility boundary (default: 0.25)",
    )
    parser.add_argument(
        "--penumbral-step-minutes",
        type=float,
        default=1.5,
        metavar="MINUTES",
        help="time step for the penumbral visibility raster scan (default: 1.5)",
    )
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="re-download the Besselian elements dataset even if already cached",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    formats = list(dict.fromkeys(args.formats or ["shp"]))  # de-dupe, preserve order
    if "all" in formats:
        formats = list(_ALL_FORMATS)

    if args.refresh_cache:
        print("Refreshing cached Besselian elements dataset...", file=sys.stderr)
        download_besselian_csv(force=True)

    try:
        elements = find_by_date(args.year, args.month, args.day)
    except LookupError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    output_dir = args.output or Path(f"eclipse_{args.year:04d}-{args.month:02d}-{args.day:02d}")

    if "shp" in formats:
        written = generate_eclipse_shapefiles(
            args.year,
            args.month,
            args.day,
            output_dir,
            step_minutes=args.step_minutes,
            elements=elements,
            penumbral_resolution_deg=args.penumbral_resolution,
            penumbral_step_minutes=args.penumbral_step_minutes,
        )
        for name, path in written.items():
            print(f"wrote {name}: {path}")

    for fmt in ("kml", "kmz"):
        if fmt not in formats:
            continue
        out_path = output_dir / f"eclipse.{fmt}"
        written_path = generate_eclipse_kml(
            args.year,
            args.month,
            args.day,
            out_path,
            step_minutes=args.step_minutes,
            elements=elements,
            penumbral_resolution_deg=args.penumbral_resolution,
            penumbral_step_minutes=args.penumbral_step_minutes,
        )
        print(f"wrote {fmt}: {written_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
