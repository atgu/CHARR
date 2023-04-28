import hail as hl
import argparse
import pandas as pd

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

def check_contam(
    contamination_underestimation_factor: float,
    output_path: str,
    output_prefix: str,
):
    # used to read from the selfSM file and calculate contamination, which gets printed out
    import csv
    import sys
    import pickle
    import subprocess
    local = subprocess.check_output("pwd", shell=True)
    local = local.decode("utf-8").split("\n")[0]
    hl.hadoop_copy(
        f"{output_path}{output_prefix}.selfSM",
        f"file://{local}/contam_temp.selfSM",
    )
    with open(f"{local}/contam_temp.selfSM", 'rt') as selfSM:
        reader = csv.DictReader(selfSM, delimiter='\t')
        i = 0
        for row in reader:
            if float(row["FREELK0"]) == 0 and float(row["FREELK1"]) == 0:
                # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
                # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
                # vcf and bam.
                sys.stderr.write(
                    "Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
                sys.exit(1)
            print(f'{output_prefix} contamination rate: {float(row["FREEMIX"]) /contamination_underestimation_factor}')
            c = float(row["FREEMIX"]) /contamination_underestimation_factor
            i = i + 1
            # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
            # and the results are not reliable.
            if i != 1:
                sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error" % (i))
                sys.exit(2)
    return c


def main():

    if args.run_freemix_sum_table:
        mixed_labels = ['HGDP00021_0.5', 'HGDP00105_0.5', 'HGDP00118_0.5', 'HGDP00281_0.5', 'HGDP00341_0.5',
                        'HGDP00529_0.5', 'HGDP00610_0.5', 'HGDP00669_0.5', 'HGDP00688_0.5', 'HGDP00689_0.5',
                        'HGDP00690_0.5', 'HGDP00703_0.5', 'HGDP00714_0.5', 'HGDP00726_0.5', 'HGDP00768_0.5',
                        'HGDP00845_0.5', 'HGDP00875_0.5', 'HGDP00905_0.5', 'HGDP00931_0.5', 'HGDP00944_0.5',
                        'HGDP01009_0.5', 'HGDP01050_0.5', 'HGDP01092_0.5', 'HGDP01096_0.5', 'HGDP01190_0.5',
                        'HGDP01243_0.5', 'HGDP01371_0.5', 'HGDP01385_0.5', 'HGDP01396_0.5', 'HGDP01406_0.5',
                        'HGDP00021_1.0', 'HGDP00105_1.0', 'HGDP00118_1.0', 'HGDP00281_1.0', 'HGDP00341_1.0',
                        'HGDP00529_1.0', 'HGDP00610_1.0', 'HGDP00669_1.0', 'HGDP00688_1.0', 'HGDP00689_1.0',
                        'HGDP00690_1.0', 'HGDP00703_1.0', 'HGDP00714_1.0', 'HGDP00726_1.0', 'HGDP00768_1.0',
                        'HGDP00845_1.0', 'HGDP00875_1.0', 'HGDP00905_1.0', 'HGDP00931_1.0', 'HGDP00944_1.0',
                        'HGDP01009_1.0', 'HGDP01050_1.0', 'HGDP01092_1.0', 'HGDP01096_1.0', 'HGDP01190_1.0',
                        'HGDP01243_1.0', 'HGDP01371_1.0', 'HGDP01385_1.0', 'HGDP01396_1.0', 'HGDP01406_1.0',
                        'HGDP00021_2.0', 'HGDP00105_2.0', 'HGDP00118_2.0', 'HGDP00281_2.0', 'HGDP00341_2.0',
                        'HGDP00529_2.0', 'HGDP00610_2.0', 'HGDP00669_2.0', 'HGDP00688_2.0', 'HGDP00689_2.0',
                        'HGDP00690_2.0', 'HGDP00703_2.0', 'HGDP00714_2.0', 'HGDP00726_2.0', 'HGDP00768_2.0',
                        'HGDP00845_2.0', 'HGDP00875_2.0', 'HGDP00905_2.0', 'HGDP00931_2.0', 'HGDP00944_2.0',
                        'HGDP01009_2.0', 'HGDP01050_2.0', 'HGDP01092_2.0', 'HGDP01096_2.0', 'HGDP01190_2.0',
                        'HGDP01243_2.0', 'HGDP01371_2.0', 'HGDP01385_2.0', 'HGDP01396_2.0', 'HGDP01406_2.0',
                        'HGDP00021_5.0', 'HGDP00105_5.0', 'HGDP00118_5.0', 'HGDP00281_5.0', 'HGDP00341_5.0',
                        'HGDP00529_5.0', 'HGDP00610_5.0', 'HGDP00669_5.0', 'HGDP00688_5.0', 'HGDP00689_5.0',
                        'HGDP00690_5.0', 'HGDP00703_5.0', 'HGDP00714_5.0', 'HGDP00726_5.0', 'HGDP00768_5.0',
                        'HGDP00845_5.0', 'HGDP00875_5.0', 'HGDP00905_5.0', 'HGDP00931_5.0', 'HGDP00944_5.0',
                        'HGDP01009_5.0', 'HGDP01050_5.0', 'HGDP01092_5.0', 'HGDP01096_5.0', 'HGDP01190_5.0',
                        'HGDP01243_5.0', 'HGDP01371_5.0', 'HGDP01385_5.0', 'HGDP01396_5.0', 'HGDP01406_5.0',
                        'HGDP00021_10.0', 'HGDP00105_10.0', 'HGDP00118_10.0', 'HGDP00281_10.0', 'HGDP00341_10.0',
                        'HGDP00529_10.0', 'HGDP00610_10.0', 'HGDP00669_10.0', 'HGDP00688_10.0', 'HGDP00689_10.0',
                        'HGDP00690_10.0', 'HGDP00703_10.0', 'HGDP00714_10.0', 'HGDP00726_10.0', 'HGDP00768_10.0',
                        'HGDP00845_10.0', 'HGDP00875_10.0', 'HGDP00905_10.0', 'HGDP00931_10.0', 'HGDP00944_10.0',
                        'HGDP01009_10.0', 'HGDP01050_10.0', 'HGDP01092_10.0', 'HGDP01096_10.0', 'HGDP01190_10.0',
                        'HGDP01243_10.0', 'HGDP01371_10.0', 'HGDP01385_10.0', 'HGDP01396_10.0', 'HGDP01406_10.0']

        contam_est = []
        split_labels = pd.DataFrame(
            [i.split("_") for i in mixed_labels],
            columns=["original", "contam_rate"],
        )
        for s in mixed_labels:
            contam_est.append(
                check_contam(1, f"{MY_BUCKET}/mixed_samples/v2/verifybamid/", s)
            )
        split_labels["freemix_score"] = contam_est
        mixed_table = hl.Table.from_pandas(split_labels)
        mixed_table.write(f"{MY_BUCKET}/hgdp_n_way_mixed_sample_freemix_score.ht")

    if args.run_vds_combiner:
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
        "--run-freemix-sum-table",
        help="Whether to run freemix summary table",
        action="store_true",
    )
    parser.add_argument(
        "--disable-sanity-check",
        help="Whether to use sanity check in verifybamID",
        action="store_true",
    )
    parser.add_argument(
        "--sample_ids",
        help="Path to a list of sample IDs to merge",
        nargs="?",
    )
    parser.add_argument(
        "--run-vds-combiner",
        help="Whether to run freemix summary table",
        action="store_true",
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
# hailctl dataproc submit wlu run_vds_combiner.py --input-gvcf-bucket gs://gnomad-wenhan/charr_simulation/mixed_samples/v2/ --output-vds-name hgdp_n_way_mixed_samples_115