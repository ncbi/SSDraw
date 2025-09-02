import sys
import argparse
import typing as T

from SSDraw.core import SSDraw
from SSDraw.multi import run_multiple_pdbs_on_one_msa


def get_args(
    argv: T.Optional[T.List[str]] = None,
) -> T.Tuple[argparse.Namespace, argparse.ArgumentParser]:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="SSDraw is a program that generates publication-quality protein secondary structure diagrams from three-dimensional protein structures. To depict relationships between secondary structure and other protein features, diagrams can be colored by conservation score, B-factor, or custom scoring.",
        epilog="",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    single_parser = subparsers.add_parser(
        "single", help="Run SSDraw on a single PDB/MSA"
    )
    single_parser.add_argument(
        "-f",
        "--fasta",
        help="(required) sequence/alignment file in fasta format",
    )
    single_parser.add_argument("-p", "--pdb", help="(required) pdb file")
    single_parser.add_argument(
        "-n",
        "--name",
        help="(required) id of the protein in the alignment file",
    )
    single_parser.add_argument(
        "-o", "--output", help="(required) name for output file"
    )
    single_parser.add_argument(
        "--SS",
        default=None,
        help="secondary structure annotation in DSSP or .horiz format. If this option is not provided, SSDraw will compute secondary structure from the given PDB file with DSSP.",
    )
    single_parser.add_argument(
        "--chain_id",
        default="A",
        help="chain id to use in pdb. Defaults to chain A.",
    )
    single_parser.add_argument(
        "--color_map",
        default=["inferno"],
        nargs="*",
        help="color map to use for heat map",
    )
    single_parser.add_argument(
        "--scoring_file",
        default=None,
        help="custom scoring file for alignment",
    )
    single_parser.add_argument(
        "--color",
        default="white",
        help="color for the image. Can be a color name (eg. white, black, green), or a hex code",
    )
    single_parser.add_argument(
        "-conservation_score",
        action="store_true",
        help="score alignment by conservation score",
    )
    single_parser.add_argument(
        "--output_file_type",
        default="png",
        help="output file type. Options: png, ps, eps, tif, svg",
    )
    single_parser.add_argument(
        "-bfactor", action="store_true", help="score by B-factor"
    )
    single_parser.add_argument(
        "-mview", action="store_true", help="color by mview color map"
    )
    single_parser.add_argument(
        "--dpi", default=600, type=int, help="dpi to use for final plot"
    )
    single_parser.add_argument(
        "--ticks", default=0, type=int, help="set ticks at every nth position"
    )
    single_parser.add_argument("--start", default=0, type=int)
    single_parser.add_argument("--end", default=0, type=int)
    single_parser.add_argument(
        "--dssp_exe",
        default="mkdssp",
        help="The path to your dssp executable. Default: mkdssp",
    )
    single_parser.add_argument(
        "--consurf",
        default="",
        help="consurf or rate4site file to color image with. If rate4site file is given, SSDraw will convert raw scores to grades.",
    )
    single_parser.add_argument(
        "--fontsize",
        default=12,
        type=int,
        help="font size for residue numbers",
    )
    single_parser.add_argument(
        "--fontcolor",
        default="black",
        type=str,
        help="font color for residue numbers",
    )
    single_parser.set_defaults(func=lambda args: SSDraw(args))

    multi_parser = subparsers.add_parser(
        "multi", help="Run SSDraw for multiple PDBs from one MSA"
    )
    multi_parser.add_argument(
        "-i", "--input", required=True, help="Name of input script"
    )
    multi_parser.add_argument(
        "-o", "--output", required=True, help="Name of output directory"
    )
    multi_parser.set_defaults(
        func=lambda args: run_multiple_pdbs_on_one_msa(args)
    )

    return parser.parse_args(argv)


def main(argv=None):
    args = get_args(argv)

    if args.command == "single":
        if args.start > args.end:
            sys.exit("--start cannot be greater than --end")
        args.func(args)
    elif args.command == "multi":
        args.func(args)
    else:
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
