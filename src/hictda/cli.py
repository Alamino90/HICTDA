"""Command-line interface for HiCTDA."""

from __future__ import annotations

import argparse

from .pipeline import analyze_hic


def build_parser():
    parser = argparse.ArgumentParser(
        prog="hictda",
        description="Persistent-homology analysis of Hi-C .hic data.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    analyze = sub.add_parser("analyze", help="Run the complete analysis.")
    analyze.add_argument("hic_file", help="Input .hic file or supported URL.")
    analyze.add_argument("--chromosome", required=True)
    analyze.add_argument("--start", type=int, default=None)
    analyze.add_argument("--end", type=int, default=None)
    analyze.add_argument("--resolution", type=int, default=50_000)
    analyze.add_argument("--normalization", default="KR")
    analyze.add_argument("--maxdim", type=int, default=1)
    analyze.add_argument("--distance-method", default="inverse")
    analyze.add_argument("--epsilon", type=float, default=1e-5)
    analyze.add_argument("--output", default="hictda_results")
    analyze.add_argument(
        "--no-arrays",
        action="store_true",
        help="Do not save NumPy arrays.",
    )
    return parser


def main():
    args = build_parser().parse_args()

    if args.command == "analyze":
        result = analyze_hic(
            hic_file=args.hic_file,
            chromosome=args.chromosome,
            start=args.start,
            end=args.end,
            resolution=args.resolution,
            normalization=args.normalization,
            maxdim=args.maxdim,
            distance_method=args.distance_method,
            epsilon=args.epsilon,
            output_dir=args.output,
            save_arrays=not args.no_arrays,
        )

        print("\nHiCTDA analysis completed.")
        print(f"Chromosome        : {result.metadata['chromosome']}")
        print(
            f"Region            : "
            f"{result.metadata['start']:,}–{result.metadata['end']:,} bp"
        )
        print(f"Resolution        : {result.metadata['resolution']:,} bp")
        print(f"Normalization     : {result.metadata['normalization']}")
        print(f"H1 features       : {result.statistics['n_h1']}")
        print(
            f"Mean persistence  : "
            f"{result.statistics['mean_persistence']:.6g}"
        )
        print(
            f"Median persistence: "
            f"{result.statistics['median_persistence']:.6g}"
        )
        print(
            f"Std persistence   : "
            f"{result.statistics['std_persistence']:.6g}"
        )
        print(f"Output directory  : {args.output}")


if __name__ == "__main__":
    main()
