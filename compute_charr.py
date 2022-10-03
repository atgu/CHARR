import argparse
from .utils import *


def main(args):
    hl.init(default_reference="GRCh38")
    hl._set_flags(grouped_aggregate_buffer_size="5")
    hl._set_flags(use_new_shuffle="1")

    run_charr(
        format=args.input_file_format,
        path=args.input_file_path,
        output_dir=args.output_dir,
        output_name=args.output_name,
        min_af=args.min_af,
        max_af=args.max_af,
        min_dp=args.min_dp,
        max_dp=args.max_dp,
        min_gq=args.min_gq,
        ref_AF_field=args.ref_AF_field,
        data_type=args.data_type,
        num_partitions=args.n_partitions,
        write_charr=True,
        extension=args.extension,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-file-format",
        type=str,
        nargs="?",
        default="vcf",
        choices=["gvcf", "mt", "vcf", "vds"],
        help="Format of the input genotype data.",
    )
    parser.add_argument(
        "--input-file-path", type=str, help="Path to import the genotype data from"
    )
    parser.add_argument(
        "--output-dir", type=str, help="Directory to write charr table to"
    )
    parser.add_argument("--output-name", type=str, help="Filename to write charr table")
    parser.add_argument(
        "--min-af",
        type=float,
        nargs="?",
        help="Filter to variants with reference allele frequency above this value",
        default=0.05,
    )
    parser.add_argument(
        "--max-af",
        type=float,
        nargs="?",
        help="Filter to variants with reference allele frequency below this value",
        default=0.95,
    )
    parser.add_argument(
        "--min-dp",
        type=float,
        nargs="?",
        help="Filter to variants with DP above this value",
        default=20,
    )
    parser.add_argument(
        "--max-dp",
        type=float,
        nargs="?",
        help="Filter to variants with DP below this value",
        default=100,
    )
    parser.add_argument(
        "--min-gq",
        type=float,
        nargs="?",
        help="Filter to variants with GQ above this value",
        default=20,
    )
    parser.add_argument("--ref-AF-field", type=str, default=None)
    parser.add_argument(
        "--data_type",
        type=str,
        nargs="?",
        default="genomes",
        choices=["genomes", "exomes"],
    )
    parser.add_argument(
        "--extension",
        type=str,
        nargs="?",
        default="ht",
        choices=["ht", "tsv"],
    )
    parser.add_argument(
        "--n_partitions",
        help="Number of desired partitions for import VCF",
        default=1000,
        type=int,
    )
    parser.add_argument("--overwrite", help="Overwrite", action="store_true")
    args = parser.parse_args()
    print(args)
    main(args)
