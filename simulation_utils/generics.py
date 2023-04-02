import hail as hl
import pipes
import hailtop.batch as hb
import numpy as np
import random
import pysam
import pandas as pd

BILLING_PROJECT = "gnomad-production"
TMP_DIR = "gs://gnomad-wenhan/tmp/contam_free_file/"
REF_FASTA_PATH = (
    "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
)
REF_FASTA_INDEX = (
    "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
)
REF_DICT = ('gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict')
MY_BUCKET = "gs://gnomad-wenhan/charr_simulation"
SAMTOOLS_IMAGE = 'staphb/samtools:latest'
ALLELES = ("A", "T", "C", "G")

chromosomes = list(map(str, range(1, 23))) + ['X']
chromosomes = [f'chr{chrom}' for chrom in chromosomes]

CONTAM_RATES = [0.005, 0.01, 0.02, 0.05, 0.1]
reference = 'GRCh38'
CHROM_LENGTHS = hl.get_reference(reference).lengths

POPs = ['afr', 'amr', 'eas', 'mid', 'nfe', 'sas']