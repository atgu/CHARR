import pysam
import pipes
import copy
import random
from random import choice
import hail as hl
import hailtop.batch as hb
import logging
import pandas as pd
import argparse
import numpy as np
from typing import List
from functools import reduce

from variant_calling.haplotype_caller import haplotype_caller_gatk
from variant_calling.merge_gvcfs import merge_vcf
from variant_calling.gvcf_index import index_gvcf
from variant_calling.variant_calling_pipeline import var_call_pipeline
from variant_calling.get_file_size import bytes_to_gb
from batch_verifybamid import *

reference = "GRCh38"
CHROM_LENGTHS = hl.get_reference(reference).lengths
BILLING_PROJECT = "gnomad-production"
MY_BUCKET = "gs://gnomad-wenhan/charr_simulation"
TMP_DIR = "gs://gnomad-wenhan/tmp/contam_free_file/"
REF_FASTA_PATH = (
    "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
)
REF_FASTA_INDEX = (
    "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"
)
REF_DICT = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict"
chromosomes = list(map(str, range(1, 23))) + ["X"]
CHROMOSOMES = [f"chr{chrom}" for chrom in chromosomes]
CONTAM_RATES = [0.005, 0.01, 0.02, 0.05, 0.1]
SAMTOOLS_IMAGE = "staphb/samtools:latest"

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("CHARR simulation pipeline")
logger.setLevel(logging.INFO)

def open_pipes_output(ref_fasta, output_name):
    pipe = pipes.Template()
    pipe.append(
        f"samtools view -C -T {ref_fasta} -h -o {output_name}  2> /dev/null",
        "-.",
    )
    f = pipe.open(
        f"temp.sam", "w"
    )
    return f

def get_read_groups(file):
    read = next(file.fetch(until_eof=True))
    rg_ind = [i for i, tuple in enumerate(read.tags) if 'RG' in tuple]
    main_rg = read.tags[rg_ind[0]][-1]
    return main_rg

def edit_read_group(read, rg_name='ERR1349727'):
    rg_ind = [i for i, tuple in enumerate(read.tags) if 'RG' in tuple]
    if len(rg_ind)>0:
        temp_tag = [list(ele) for ele in read.tags]
        temp_tag[rg_ind[0]][-1] = rg_name
        read.tags = [tuple(ele) for ele in temp_tag]
    return read

# def mixing_many_samples(
#         input_list: List,
#         output_list: List,
#         sample_list: List,
#         input_ref_fasta: str,
#         chromosome: str,
#         contam_rate: float,
# ):
#     save_pos=[]
#     save_next_reads_len=[]
#     save_current_reads_len=[]
#     save_current_reads_actually_at_position_len=[]
#
#     inputs = [pysam.AlignmentFile(cram_input["cram"], mode="rc", reference_filename=input_ref_fasta) for
#               cram_input in input_list]
#     outputs = [open_pipes_output(ref_fasta=input_ref_fasta, output_name=output_list[i]) for i in
#                range(len(sample_list))]
#     # main_rgs = [get_read_groups(input) for input in inputs]
#     # print(main_rgs)
#     next_reads = [[edit_read_group(next(inputs[i].fetch(chromosome, until_eof=True)))] for i in range(len(inputs))]
#     current_pos = [next_read[0].pos for next_read in next_reads]
#     last_reads = [False for input in inputs]
#
#     for pos in range(CHROM_LENGTHS[chromosome]):
#         save_pos.append(pos)
#         if (not any(x == pos for x in current_pos)):
#             continue
#         if any(last_reads):
#             break
#         next_reads = [[read for read in reads if read.pos >= pos] for reads in next_reads]
#         save_next_reads_len.append(np.mean([len(x) for x in next_reads]))
#         current_reads = copy.deepcopy(next_reads)
#         for sample in range(len(inputs)):
#             if len(next_reads[sample]) > 0 and next_reads[sample][0].pos != pos:
#                 continue
#             while True:
#                 # print(f'current sample: {samples[sample]}')
#                 # print(f'current read: {next_read}')
#                 # print(f'current position: {pos}')
#                 # print(f'current read position: {next_read.pos}')
#                 try:
#                     next_read = edit_read_group(next(inputs[sample]))
#                 except StopIteration:
#                     last_reads[sample] = True
#                     break
#
#                 if next_read.pos == pos and (not last_reads[sample]):
#                     current_reads[sample].append(next_read)
#                 elif not last_reads[sample]:
#                     next_reads[sample].append(next_read)
#                     break
#         save_current_reads_len.append(np.mean([len(x) for x in current_reads]))
#         current_pos = [next_read[-1].pos for next_read in next_reads]
#
#         # contamination
#         current_reads_actually_at_position = [[read for read in reads if read.pos == pos] for reads in current_reads]
#         save_current_reads_actually_at_position_len.append(np.mean([len(x) for x in current_reads_actually_at_position]))
#         for sample in range(len(inputs)):
#             for read in current_reads_actually_at_position[sample]:
#                 contaminate = np.random.binomial(1, contam_rate)
#                 if not contaminate:
#                     outputs[sample].write(read.tostring() + "\n")
#                 else:
#                     indexes = list(
#                         filter(lambda x: x not in [sample] and len(current_reads_actually_at_position[x]) > 0,
#                                range(0, len(inputs))))
#                     if len(indexes)==0:
#                         continue
#                     random_ind = choice(indexes)
#                     random_sample_reads = current_reads_actually_at_position[random_ind]
#                     random_read_ind = random.randint(0, len(random_sample_reads) - 1)
#                     random_read = random_sample_reads[random_read_ind]
#                     outputs[sample].write(random_read.tostring() + "\n")
#     for output in outputs:
#         output.close()
#     for input in inputs:
#         input.close()
#     df = pd.DataFrame(data=zip(save_pos, save_current_reads_len, save_next_reads_len, save_current_reads_actually_at_position_len),
#                       columns=['Position', 'Current_reads_length', 'Next_reads_length', 'Current_reads_actually_at_position_length'])
#     with hl.hadoop_open(f'{MY_BUCKET}/mixed_samples/v2/reads_summary_{chromosome}.csv', 'w') as f:
#         df.to_csv(f)

def mixing_many_samples(
        input_list: List,
        output_list: List,
        sample_list: List,
        input_ref_fasta: str,
        chromosome: str,
        contam_rate: float,
):
    inputs = [pysam.AlignmentFile(cram_input["cram"], mode="rc", reference_filename=input_ref_fasta) for
              cram_input in input_list]
    outputs = [open_pipes_output(ref_fasta=input_ref_fasta, output_name=output_list[i]) for i in
               range(len(sample_list))]
    # main_rgs = [get_read_groups(input) for input in inputs]
    # print(main_rgs)
    next_reads = [[edit_read_group(next(inputs[i].fetch(chromosome, until_eof=True)))] for i in range(len(inputs))]
    current_pos = [next_read[0].pos for next_read in next_reads]
    last_reads = [False for input in inputs]

    for pos in range(CHROM_LENGTHS[chromosome]):
        if (not any(x == pos for x in current_pos)):
            continue
        if any(last_reads):
            break
        next_reads = [[read for read in reads if read.pos >= pos] for reads in next_reads]
        current_reads = copy.deepcopy(next_reads)
        for sample in range(len(inputs)):
            if len(next_reads[sample]) > 0 and next_reads[sample][0].pos != pos:
                continue
            while True:
                # print(f'current sample: {samples[sample]}')
                # print(f'current read: {next_read}')
                # print(f'current position: {pos}')
                # print(f'current read position: {next_read.pos}')
                try:
                    next_read = edit_read_group(next(inputs[sample]))
                except StopIteration:
                    last_reads[sample] = True
                    break

                if next_read.pos == pos and (not last_reads[sample]):
                    current_reads[sample].append(next_read)
                elif not last_reads[sample]:
                    next_reads[sample].append(next_read)
                    break
        current_pos = [next_read[-1].pos for next_read in next_reads]

        # contamination
        current_reads_actually_at_position = [[read for read in reads if read.pos == pos] for reads in current_reads]
        for sample in range(len(inputs)):
            for read in current_reads_actually_at_position[sample]:
                contaminate = np.random.binomial(1, contam_rate)
                if not contaminate:
                    outputs[sample].write(read.tostring() + "\n")
                else:
                    indexes = list(
                        filter(lambda x: x not in [sample] and len(current_reads_actually_at_position[x]) > 0,
                               range(0, len(inputs))))
                    if len(indexes)==0:
                        continue
                    random_ind = choice(indexes)
                    random_sample_reads = current_reads_actually_at_position[random_ind]
                    random_read_ind = random.randint(0, len(random_sample_reads) - 1)
                    random_read = random_sample_reads[random_read_ind]
                    outputs[sample].write(random_read.tostring() + "\n")
    for output in outputs:
        output.close()
    for input in inputs:
        input.close()


def main():
    backend = hb.ServiceBackend(
        billing_project=BILLING_PROJECT,
        remote_tmpdir=TMP_DIR,
    )
    b = hb.Batch(
        name=f"Decontam-Mixing-Cram-Files",
        requester_pays_project="daly-ibd",
        default_python_image="us-central1-docker.pkg.dev/broad-mpg-gnomad/wlu/hail/hail-pysam-samtools:4.16.0",
        backend=backend,
    )

    input_ref_fasta = b.read_input_group(
        **{
            "fasta": REF_FASTA_PATH,
            "fasta.fai": REF_FASTA_INDEX,
            "dict": REF_DICT,
        }
    )

    logger.info("Loading sample info...")
    sample_paths = pd.read_csv(args.input_sample_info, sep=",")
    sample_ids = pd.read_csv(args.selected_samples, sep=",")
    sample_paths = sample_paths[sample_paths["s"].isin(sample_ids["hgdp_id"])]
    samples = []
    cram_files = {}
    for sample in sample_paths.iterrows():
        samples.append(sample[1][0])
        cram_files[sample[1][0]] = (f"{sample[1][1]}", f"{sample[1][1]}.crai")

    # if args.test:
    #     samples = samples[:4]

    mixed_labels = []
    for contam_rate in CONTAM_RATES:
        print(contam_rate)
        mixed_crams = {}
        mixed_gvcfs = {}
        mix_gvcf_job_depend_list = []
        mix_cram_job_depend_list = []
        for sample in samples:
            mixing_samples_label = f"{sample}_{contam_rate * 100}"
            mixed_crams[mixing_samples_label] = {}
            mixed_gvcfs[mixing_samples_label] = {}

        if args.run_mixing:
            for chromosome in CHROMOSOMES:
                if args.test:
                    chromosome = 'chr21'
                output_mix_cram_paths = [f'{MY_BUCKET}/mixed_samples/v2/crams/cram_by_chrom/{sample}/contam_rate_{contam_rate*100}/{sample}_{chromosome}_{contam_rate*100}.cram' for sample in samples]
                output_mix_vcf_paths = [
                    f'{MY_BUCKET}/mixed_samples/v2/variant-calling/{sample}_{contam_rate * 100}/{sample}_{contam_rate * 100}_{chromosome}.g.vcf'
                    for sample in samples]
                j_reheader_depend_on = None
                if any([not hl.hadoop_exists(path) for path in output_mix_cram_paths]):
                    logger.info(
                        f"Mixing contamination free crams: {chromosome}_{contam_rate*100}_percent_contamination..."
                    )
                    j_mix = b.new_python_job(
                        name=f"Mix_all_samples_{chromosome}_{contam_rate*100}",
                        attributes={
                            "contam_rate": f"{contam_rate * 100}\%",
                            "chromosome": chromosome,
                            "job_type": "mixing_samples",
                        },
                    )
                    if chromosome == 'chr2':
                        j_mix.storage("100Gi").memory("50Gi")
                    else:
                        j_mix.storage("100Gi").memory("15Gi")
                    cram_paths = [
                        f"{MY_BUCKET}/contam_free/crams/cram_by_chrom/{sample_id}/{sample_id}_contam_free_{chromosome}.cram"
                        for sample_id in samples]
                    inputs = [b.read_input_group(**{"cram": path, "cram.crai": f"{path}.crai"}) for path in
                                         cram_paths]
                    outputs = [j_mix[f'{sample}_{chromosome}'] for sample in samples]

                    j_mix.call(
                        mixing_many_samples,
                        input_list=inputs,
                        output_list=outputs,
                        sample_list=samples,
                        input_ref_fasta=input_ref_fasta['fasta'],
                        chromosome=chromosome,
                        contam_rate=contam_rate
                    )
                    j_reheader_depend_on = j_mix
                    mix_cram_job_depend_list.append(j_mix)

                    for i in range(len(samples)):
                        mixing_samples_label = f"{samples[i]}_{contam_rate * 100}"
                        output_cram_path = output_mix_cram_paths[i]
                        b.write_output(outputs[i], output_cram_path)
                        mixed_crams[mixing_samples_label][chromosome] = outputs[i]


                for i in range(len(samples)):
                    sample = samples[i]
                    mixing_samples_label = f"{sample}_{contam_rate * 100}"
                    j_gvcf_depend_on = None
                    if not hl.hadoop_exists(f"{output_mix_cram_paths[i]}.crai"):
                        j_reheader = b.new_job(
                            name=f"Reheader_contam_free_file_{mixing_samples_label}_{chromosome}",
                            attributes={
                                "sample": sample,
                                "contam_rate": f"{contam_rate * 100}\%",
                                "chromosome": chromosome,
                                "job_type": "run_reheader",
                            },
                        )
                        j_reheader.image(SAMTOOLS_IMAGE).storage("40Gi").memory("15Gi")

                        if j_reheader_depend_on is not None:
                            j_reheader.depends_on(j_reheader_depend_on)
                        logger.info(
                            f"Reheadering contamination free cram: {sample}_{chromosome}_{contam_rate * 100}_percent_contamination..."
                        )
                        header_cram_file = b.read_input_group(
                            cram=cram_files['HGDP00021'][0], index=cram_files['HGDP00021'][1]
                        )
                        tmp_mix_cram = b.read_input(output_mix_cram_paths[i])
                        j_reheader.command(
                            f"samtools view -H {header_cram_file['cram']} > {j_reheader.header}"
                        )
                        j_reheader.command(
                            f"samtools reheader {j_reheader.header} {tmp_mix_cram} > {j_reheader.ofile1}"
                        )
                        j_reheader.command(
                            f"samtools index {j_reheader.ofile1} -o {j_reheader.ofile2}"
                        )
                        mix_cram_job_depend_list.append(j_reheader)
                        b.write_output(j_reheader.ofile1, output_mix_cram_paths[i])
                        b.write_output(j_reheader.ofile2, f"{output_mix_cram_paths[i]}.crai")
                        mixed_crams[mixing_samples_label][chromosome] = j_reheader.ofile1
                        j_gvcf_depend_on = j_reheader
                    else:
                        path = output_mix_cram_paths[i]
                        input_mixed_cram = b.read_input_group(**{"cram": path, "cram.crai": f"{path}.crai"})
                        mixed_crams[mixing_samples_label][chromosome] = input_mixed_cram['cram']

                    if not hl.hadoop_exists(output_mix_vcf_paths[i]):
                        logger.info(
                            f"Running haplotype caller: {sample}_{chromosome}_{contam_rate*100}_percent_contamination..."
                        )
                        input_chrom_cram_file = b.read_input_group(
                            **{
                                "cram": output_mix_cram_paths[i],
                                "cram.crai": f"{output_mix_cram_paths[i]}.crai",
                            }
                        )
                        tmp_gvcf, tmp_job = haplotype_caller_gatk(
                            b=b,
                            depend_on=j_gvcf_depend_on,
                            input_bam=input_chrom_cram_file,
                            ref_fasta=input_ref_fasta,
                            interval_list_file=chromosome,
                            out_dir=f"{MY_BUCKET}/mixed_samples/v2",
                            contamination=0.0,
                            bam_filename_no_ext=f"{sample}_{contam_rate * 100}_{chromosome}",
                            storage=40,
                            interval_list_name=None,
                            memory=26,
                        )
                        mix_gvcf_job_depend_list.append(tmp_job)
                        mixed_gvcfs[mixing_samples_label][chromosome] = tmp_gvcf
                    else:
                        tmp_gvcf = b.read_input(output_mix_vcf_paths[i])
                        mixed_gvcfs[mixing_samples_label][chromosome] = tmp_gvcf

                if args.test:
                    break

        for sample in samples:
            mixing_samples_label = f"{sample}_{contam_rate * 100}"
            mixed_labels.append(mixing_samples_label)
            output_mix_full_cram_path = (
                f"{MY_BUCKET}/mixed_samples/v2/crams/{mixing_samples_label}_mixed.cram"
            )
            j_cat_depend_on = None
            if not hl.hadoop_exists(output_mix_full_cram_path):
                logger.info(f"Concatenating mixed crams: {mixing_samples_label}...")
                j_cat = b.new_job(
                    name=f"Concatenate_mixed_file_{mixing_samples_label}",
                    attributes={
                        "sample": sample,
                        "contam_rate": f"{contam_rate * 100}\%",
                        "job_type": "concatenate_mixed_crams",
                    },
                )
                if len(mix_cram_job_depend_list) > 0:
                    j_cat.depends_on(*mix_cram_job_depend_list)
                j_cat.image(SAMTOOLS_IMAGE).storage("60Gi").memory("10Gi")
                tmp_cram_lst = reduce(
                    lambda x, y: x + " " + y, mixed_crams[mixing_samples_label].values()
                )
                header_cram_file = b.read_input_group(
                    cram=cram_files['HGDP00021'][0], index=cram_files['HGDP00021'][1]
                )
                j_cat.command(
                    f"samtools view -H {header_cram_file['cram']} > {j_cat.header}"
                )
                j_cat.command(
                    f"samtools cat -h {j_cat.header} -o {j_cat.ofile1} {tmp_cram_lst}"
                )
                j_cat.command(f"samtools index {j_cat.ofile1} -o {j_cat.ofile2}")
                b.write_output(j_cat.ofile1, output_mix_full_cram_path)
                b.write_output(j_cat.ofile2, f"{output_mix_full_cram_path}.crai")
                j_cat_depend_on = j_cat

            if hl.hadoop_exists(output_mix_full_cram_path) and not hl.hadoop_exists(
                f"{MY_BUCKET}/mixed_samples/v2/verifybamid/{mixing_samples_label}.selfSM"
            ):
                logger.info(f"Running VerifyBam ID: {mixing_samples_label}...")
                run_verifybamid(
                    b=b,
                    input_cram_path=output_mix_full_cram_path,
                    input_crai_path=f"{output_mix_full_cram_path}.crai",
                    cram_project_id="broad-mpg-gnomad",
                    ref_fasta_path=REF_FASTA,
                    contamination_sites_path=CONTAM_SITES,
                    output_path=f"{MY_BUCKET}/mixed_samples/v2/verifybamid/",
                    output_prefix=mixing_samples_label,
                    depend_on=j_cat_depend_on,
                    disable_sanity_check=args.disable_sanity_check,
                )
            if not hl.hadoop_exists(
                f"{MY_BUCKET}/mixed_samples/v2/merged-gvcf/{mixing_samples_label}.g.vcf.gz"
            ) and hl.hadoop_exists(f'{MY_BUCKET}/mixed_samples/v2/variant-calling/{sample}_{contam_rate * 100}/{sample}_{contam_rate * 100}_chr1.g.vcf'):
                logger.info(f"Running merge gvcfs: {mixing_samples_label}...")
                gvcfs_to_merge = hl.utils.hadoop_ls(
                    f'{MY_BUCKET}/mixed_samples/v2/variant-calling/{mixing_samples_label}/*.vcf')
                gvcfs_list = []
                gvcfs_sizes_sum = 0
                for file in gvcfs_to_merge:
                    gvcfs_list.append(file['path'])
                    gvcfs_sizes_sum += bytes_to_gb(file['path'])
                merge_disk_size = round(gvcfs_sizes_sum * 2.5) + 15
                merged_vcf, j_merge = merge_vcf(
                    b=b,
                    gvcf_list=gvcfs_list,
                    depend_on=mix_gvcf_job_depend_list,
                    storage=merge_disk_size,
                    output_vcf_name=f"{mixing_samples_label}",
                    out_dir=f"{MY_BUCKET}/mixed_samples/v2",
                    memory="50",
                )
            else:
                merged_vcf = b.read_input(
                    f"{MY_BUCKET}/mixed_samples/v2/merged-gvcf/{mixing_samples_label}.g.vcf.gz"
                )

            if not hl.hadoop_exists(
                f"{MY_BUCKET}/mixed_samples/v2/merged-gvcf/{mixing_samples_label}.g.vcf.gz.tbi"
            ) and hl.hadoop_exists(
                f"{MY_BUCKET}/mixed_samples/v2/merged-gvcf/{mixing_samples_label}.g.vcf.gz"
            ):
                logger.info(f"Indexing gvcf: {mixing_samples_label}...")
                index_gvcf(
                    b=b,
                    input_vcf=merged_vcf,
                    output_vcf_ind_name=f"{mixing_samples_label}",
                    out_dir=f"{MY_BUCKET}/mixed_samples/v2",
                    storage="15",
                    memory="15"
                )
        if args.test:
            break
    b.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-sample-info",
        help="Path to a csv of paths to cram and gvcf files with sample IDs",
        nargs="?",
    )
    parser.add_argument(
        "--selected-samples",
        help="Path to a list of sample IDs to run and mix",
        nargs="?",
    )
    parser.add_argument(
        "--test",
        help="Whether to run test",
        action="store_true",
    )
    parser.add_argument(
        "--run_mixing",
        help="Whether to run mixing chromosomes",
        action="store_true",
    )
    parser.add_argument(
        "--run-merge-gvcf-files",
        help="Whether to run freemix summary table",
        action="store_true",
    )
    args = parser.parse_args()
    print(args)
    main()