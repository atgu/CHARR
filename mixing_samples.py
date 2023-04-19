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

def edit_read_group(read, rg_name):
    rg_ind = [i for i, tuple in enumerate(read.tags) if 'RG' in tuple]
    if len(rg_ind)>0:
        temp_tag = [list(ele) for ele in read.tags]
        temp_tag[rg_ind[0]][-1] = rg_name
        read.tags = [tuple(ele) for ele in temp_tag]
    return read

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
    main_rgs = [get_read_groups(input) for input in inputs]
    next_reads = [[next(input.fetch(chromosome, until_eof=True))] for input in inputs]
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
                    next_read = edit_read_group(next(inputs[sample]), main_rgs[sample])
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
                               range(0, len(inputs) - 1)))
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

    if args.test:
        samples = samples[:3]

    for i in range(len(samples)):
        sample = samples[i]
        input_cram_file = b.read_input_group(
            cram=cram_files[sample][0], index=cram_files[sample][1]
        )
        j_header = b.new_job(
            name=f"get_header_{sample}",
            attributes={"sample_id": sample, "job_type": "get_header"},
        )
        j_header.image(SAMTOOLS_IMAGE).storage("15Gi").memory('15Gi')
        j_header.command(
            f"samtools view -H {input_cram_file['cram']} > {j_header[f'{sample}']}"
        )

    for contam_rate in CONTAM_RATES:
        for chromosome in CHROMOSOMES:
            if args.test:
                chromosome = 'chr21'
            output_mix_cram_paths = [f'{MY_BUCKET}/mixed_samples/v2/crams/cram_by_chrom/{sample}/contam_rate_{contam_rate*100}/{sample}_{chromosome}_{contam_rate*100}.cram' for sample in samples]
            if any([not hl.hadoop_exists(path) for path in output_mix_cram_paths]):
                logger.info(
                    f"Mixing contamination free crams: {chromosome}_{contam_rate*100}_percent_contamination..."
                )
                j_mix = b.new_python_job(
                    name=f"Mix_all_samples_{chromosome}",
                    attributes={
                        "contam_rate": f"{contam_rate * 100}\%",
                        "chromosome": chromosome,
                        "job_type": "mixing_samples",
                    },
                )
                j_mix.storage("50Gi").memory("10Gi")
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

                for i in range(len(samples)):
                    output_cram_path = output_mix_cram_paths[i]
                    b.write_output(outputs[i], output_cram_path)

            elif any([not hl.hadoop_exists(f"{path}.crai") for path in output_mix_cram_paths]):
                logger.info(
                    f"Reheadering contamination free cram: {chromosome}_{contam_rate*100}_percent_contamination..."
                )
                j_reheader = b.new_job(
                    name=f"Reheader_contam_free_file_{chromosome}",
                    attributes={
                        "contam_rate": f"{contam_rate * 100}\%",
                        "chromosome": chromosome,
                        "job_type": "run_reheader",
                    },
                )
                j_reheader.image(SAMTOOLS_IMAGE).storage("20Gi").memory("3.75Gi")
                j_reheader.depends_on(j_header)

                for i in range(len(samples)):
                    if not hl.hadoop_exists(f"{output_mix_cram_paths[i]}.crai"):
                        tmp_mix_cram = b.read_input(output_mix_cram_paths[i])
                        sample = samples[i]
                        j_reheader.command(
                            f"samtools reheader {j_header[f'{sample}']} {tmp_mix_cram} > {j_reheader.ofile1}"
                        )
                        j_reheader.command(
                            f"samtools index {j_reheader.ofile1} -o {j_reheader.ofile2}"
                        )
                        b.write_output(j_reheader.ofile1, output_mix_cram_paths[i])
                        b.write_output(j_reheader.ofile2, f"{output_mix_cram_paths[i]}.crai")
            if args.test:
                break
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
    args = parser.parse_args()
    print(args)
    main()
