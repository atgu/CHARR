"""
Create job to run verifyBamID
"""
import os
import logging
import argparse
from math import ceil
import hail as hl
import hailtop.batch as hb
from hailtop.batch.job import Job
from typing import List, Optional, Dict
from google.cloud import storage
import pandas as pd

# ResourceGroup: set of files with multiple extensions
IMAGE="us.gcr.io/broad-gotc-prod/verify-bam-id:1.0.1-c1cba76e979904eb69c31520a0d7f5be63c72253-1639071840"
REF_FASTA="gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
CONTAM_SITES="gs://gcp-public-data--broad-references/hg38/v0/contamination-resources/1000g/1000g.phase3.100k.b38.vcf.gz.dat"
BILLING_PROJECT = 'gnomad-production'

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("VerifyBamIDn")
logger.setLevel(logging.INFO)


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


def run_verifybamid(
        b,
        input_cram_path: str,
        input_crai_path: str,
        cram_project_id: str,
        ref_fasta_path: str,
        contamination_sites_path: str,
        output_path: str,
        output_prefix: str,
        depend_on: hb.batch.job=None,
        disable_sanity_check: bool=True,
):
    if hl.hadoop_exists(input_cram_path):
        bucket_name = input_cram_path.split('/')[2]
        client: storage.client.Client = storage.Client(cram_project_id)
        bucket: storage.bucket.Bucket = client.get_bucket(bucket_name)
        bam_size = bucket.get_blob(input_cram_path.replace(f'gs://{bucket_name}/', '')).size
        logger.info(
            f'Cram file size: {bam_size/(1024**3)} GiB ... '
        )
    else:
        raise ValueError(
            f'Input cram file {input_cram_path} does not exist!'
        )

    if hl.hadoop_exists(ref_fasta_path):
        client: storage.client.Client = storage.Client('bigquery-public-data')
        bucket: storage.bucket.Bucket = client.get_bucket('gcp-public-data--broad-references')
        ref_size = bucket.get_blob('hg38/v0/Homo_sapiens_assembly38.fasta').size
        logger.info(
            f'Reference file size: {ref_size/(1024**3)} GiB ... '
        )
    else:
        raise ValueError(
            f'Input reference file {ref_fasta_path} does not exist!'
        )


    ncpu = 8  # ~ 8G/core ~ 64G
    disk_size = ceil((bam_size + ref_size)/ (1024**3)) + 30
    # disk_size = 40
    # path, name = input_cram_path.rsplit('/', 1)
    # output_prefix = name[:-5].split('.')[0]

    j = b.new_job(f'Run_VerifyBamID_{output_prefix}')
    j.image(IMAGE).cpu(ncpu).storage(f'{disk_size}G').memory('highmem')
    if depend_on is not None:
        j.depends_on(depend_on)

    logger.info(
        f'Requesting storage: {disk_size} GiB ... '
    )
    j.declare_resource_group(
        **{f'{output_prefix}': {'selfSM': f'{output_prefix}.selfSM'}}
    )

    input_cram = b.read_input_group(
            **{
                'file': input_cram_path,
                'index': input_crai_path,
            }
        )

    ref_fasta = b.read_input_group(
            **{
                'file': ref_fasta_path,
                'index': ref_fasta_path + '.fai',
            }
        )

    contamination_sites = b.read_input_group(
            **{
                'ud': contamination_sites_path + '.UD',
                'mu': contamination_sites_path + '.mu',
                'bed': contamination_sites_path + '.bed'
            }
        )

    j.command(
        f"""set -euo pipefail
        ./VerifyBamID \\
        --Verbose \\
        --NumPC 4 \\
        --Output {j[output_prefix]} \\
        --BamFile {input_cram.file} \\
        --Reference {ref_fasta['file']} \\
        --UDPath {contamination_sites['ud']} \\
        --MeanPath {contamination_sites['mu']} \\
        --BedPath {contamination_sites['bed']} \\
        {"--DisableSanityCheck" if disable_sanity_check else ""} ;
        """
    )

    if output_path:
        b.write_output(j[output_prefix], f'{output_path}{output_prefix}')
    return j[output_prefix]


def main(args):
    hl.init(default_reference="GRCh38")

    try:
        tmp_bucket = f'gs://{args.bucket}/verifybamid/'


        if args.test:
            logger.info(
                f'Starting hail Batch with the project {BILLING_PROJECT}, '
                f'bucket {tmp_bucket}'
            )
            backend = hb.ServiceBackend(
                billing_project=BILLING_PROJECT,
                remote_tmpdir=tmp_bucket,
            )

            b = hb.Batch(
                'VerifyBamID',
                backend=backend,
            )
            run_verifybamid(
                b=b,
                input_cram_path=args.input_cram_path,
                input_crai_path=args.input_crai_path,
                ref_fasta_path=REF_FASTA,
                contamination_sites_path=CONTAM_SITES,
                output_path=tmp_bucket,
                disable_sanity_check=args.disable_sanity_check,
            )

            b.run()
        else:
            ht = hl.import_table('gs://gnomad-wenhan/verifybamid/gnomad_v3_hgdp_low_contam.txt').key_by('s')
            cram_data = ht.aggregate(hl.struct(cram_path=hl.agg.collect(ht.final_cram_path), crai_path=hl.agg.collect(ht.final_crai_path), sample_ids=hl.agg.collect(ht.s)))
            if args.run_verifybamid:
                logger.info(
                    f'Starting hail Batch with the project {BILLING_PROJECT}, '
                    f'bucket {tmp_bucket}'
                )
                backend = hb.ServiceBackend(
                    billing_project=BILLING_PROJECT,
                    remote_tmpdir=tmp_bucket,
                )

                b = hb.Batch(
                    'VerifyBamID',
                    backend=backend,
                )

                for cram_path, crai_path in zip(cram_data.cram_path, cram_data.crai_path):
                    run_verifybamid(
                        b=b,
                        input_cram_path=cram_path,
                        input_crai_path=crai_path,
                        ref_fasta_path=REF_FASTA,
                        contamination_sites_path=CONTAM_SITES,
                        output_path=f'{tmp_bucket}contam_selSM/',
                        disable_sanity_check=args.disable_sanity_check,
                    )

                b.run()
            elif args.get_sum_table:
                contam_names = [x.split('/')[-1][:-5] for x in cram_data.cram_path]
                contam_paths = [tmp_bucket + x.split('/')[-1][:-5] + '.selfSM' for x in cram_data.cram_path]
                contam_est = []
                for name in contam_names:
                    contam_est.append(check_contam(args.contamination_underestimation_factor, tmp_bucket, name))
                cram_data = cram_data.annotate(recomputed_contam_rate = contam_est,
                                               contam_path = contam_paths)
                cram_data.write(f'{tmp_bucket}/{args.output_contam_table_name}.ht', overwrite=True)

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(f"{tmp_bucket}log/")



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-cram-path')
    parser.add_argument('--input-crai-path')
    parser.add_argument('--cram-project-id', default="broad-mpg-gnomad")
    parser.add_argument('--bucket', default="gnomad-wenhan")
    parser.add_argument('--disable-sanity-check', action="store_true")
    parser.add_argument('--run-verifybamid', action="store_true")
    parser.add_argument('--get-sum-table', action="store_true")
    parser.add_argument('--test', action="store_true")
    parser.add_argument('--contamination-underestimation-factor', type=float, default=1)
    parser.add_argument('--output_contam_table_name')


    args = parser.parse_args()
    main(args)


# python3 batch_verifybamid.py \
# --test \
# --input-cram-path gs://gnomad-wenhan/verifybamid/ELGH7934070_30702_1_62.cram \
# --input-crai-path gs://gnomad-wenhan/verifybamid/ELGH7934070_30702_1_62.cram.crai


# python3 batch_verifybamid.py \
# --test \
# --input-cram-path gs://gnomad-wenhan/verifybamid/cram/HGDP00967.alt_bwamem_GRCh38DH.20181023.Yakut.cram \
# --input-crai-path gs://gnomad-wenhan/verifybamid/cram/HGDP00967.alt_bwamem_GRCh38DH.20181023.Yakut.cram.crai

# python3 batch_verifybamid.py --run-verifybamid