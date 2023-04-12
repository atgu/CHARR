import hail as hl
import argparse

BILLING_PROJECT = "gnomad-production"
MY_BUCKET = "gs://gnomad-wenhan/charr_simulation"
TMP_DIR = "gs://gnomad-wenhan/tmp/contam_free_file/"
SAMTOOLS_IMAGE = "staphb/samtools:latest"
# SAMTOOLS_IMAGE = 'us-central1-docker.pkg.dev/broad-mpg-gnomad/wlu/hail/samtools:1.0'

REF_FASTA_PATH = (
    "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
)
REF_FASTA_INDEX = (
    "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
)
REF_DICT = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"


CALLING_INTERVAL_LIST = "gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list"
EVALUATION_INTERVAL_LIST = "gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list"
DBSNP_VCF = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
DBSNP_VCF_INDEX = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx"

ALLELES = ("A", "T", "C", "G")

chromosomes = list(map(str, range(1, 23))) + ["X"]
CHROMOSOMES = [f"chr{chrom}" for chrom in chromosomes]


reference = "GRCh38"
CHROM_LENGTHS = hl.get_reference(reference).lengths
CONTAM_RATES = [0.005, 0.01, 0.02, 0.05, 0.1]

POPs = ["afr", "amr", "eas", "mid", "nfe", "sas"]


def main():
    gvcfs_to_combine = hl.utils.hadoop_ls(f'{args.input_gvcf_bucket}/merged-gvcf/*.vcf.gz')
    gvcfs_list = []
    gvcf_names = []
    for file in gvcfs_to_combine:
        if hl.hadoop_exists(f"{file['path']}.tbi"):
            gvcfs_list.append(file['path'])
            gvcf_names.append(file['path'].split("/")[-1][:-9])
    print(gvcf_names)

    combiner = hl.vds.new_combiner(
        output_path=f"{MY_BUCKET}/{args.output_vds_name}.vds",
        temp_path=TMP_DIR,
        gvcf_paths=gvcfs_list,
        gvcf_sample_names=gvcf_names,
        gvcf_external_header="gs://gnomad-wenhan/charr_simulation/header.txt",
        use_genome_default_intervals=True,
        reference_genome="GRCh38",
    )

    ## Run combiner
    combiner.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample_ids",
        help="Path to a list of sample IDs to merge",
        nargs="?",
    )
    parser.add_argument(
        "--input-gvcf-bucket",
        help="Input gvcf bucket",
        nargs="?",
        default="gs://gnomad-wenhan/charr_simulation/contam_free/"
    )
    parser.add_argument(
        "--output-vds-name",
        help="Output VDS name",
        nargs="?",
    )
    args = parser.parse_args()
    print(args)
    main()

# hailctl dataproc submit wlu run_vds_combiner.py --output-vds-name hgdp_decontaminated_28_samples
# hailctl dataproc submit wlu run_vds_combiner.py --input-gvcf-bucket gs://gnomad-wenhan/charr_simulation/mixed_samples/ --output-vds-name hgdp_mixed_samples_full