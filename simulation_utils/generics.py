import hail as hl
import pipes
import hailtop.batch as hb
import numpy as np
import random
import pysam
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