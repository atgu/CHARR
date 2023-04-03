#!/usr/bin/env python3

__author__ = 'Wenhan Lu'


import argparse
from functools import reduce

import pandas as pd

from variant_calling.haplotype_caller import haplotype_caller_gatk
from variant_calling.merge_gvcfs import merge_vcf
from variant_calling import var_call_pipeline

from batch_verifybamid import *
from simulation_utils.cram_decontam_func import *
from simulation_utils.cram_mixing_func import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("CHARR simulation pipeline")
logger.setLevel(logging.INFO)

# Get five samples from each ancestry: EUR, AFR, AMR, EAS, SAS, MID -> 30 samples
# Apply 5 contamination rates: 0.005, 0.01, 0.02, 0.05, 0.1
# Overall: 450*5 = 2250 mixed new samples

def main():
    backend = hb.ServiceBackend(
        billing_project=BILLING_PROJECT,
        remote_tmpdir=TMP_DIR,
    )
    b = hb.Batch(
        name=f"Decontam-Mixing-Cram-Files",
        requester_pays_project="daly-ibd",
        default_python_image="us-central1-docker.pkg.dev/broad-mpg-gnomad/wlu/hail/hail-pysam-samtools:latest",
        backend=backend,
    )

    logger.info("Loading sample info...")
    sample_paths = pd.read_csv(args.input_sample_info, sep=",")
    sample_ids = pd.read_csv(args.selected_samples, sep=",")
    sample_paths = sample_paths[sample_paths["s"].isin(sample_ids["hgdp_id"])]
    samples = []
    cram_files = {}
    gvcf_files = {}
    sample_pops = {}
    for sample in sample_paths.iterrows():
        samples.append(sample[1][0])
        cram_files[sample[1][0]] = (f"{sample[1][1]}", f"{sample[1][1]}.crai")
        gvcf_files[sample[1][0]] = (f"{sample[1][2]}", f"{sample[1][2]}.tbi")
        sample_pops[sample[1][0]] = sample[1][3]

    input_ref_fasta = b.read_input_group(
        fasta=REF_FASTA_PATH, index=REF_FASTA_INDEX, dict=REF_DICT
    )

    full_contam_free_crams = []
    contam_free_crams = {}
    i = 0
    for s in samples:
        i += 1
        sample_id = s
        contam_free_crams[sample_id] = {}
        logger.info(f"-------Decontaminating {sample_id}-------")
        input_cram_file = b.read_input_group(
            cram=cram_files[s][0], index=cram_files[s][1]
        )
        if not args.skip_gvcf_dict:
            logger.info(f"Generating gvcf dict - {sample_id}...")
            gvcf_path = gvcf_files[s]
            gvcf_file_name = gvcf_path[0].split("/")[-1][:-7]
            gvcf_dict, j_gvcf = run_gvcf_dict(b, gvcf_file_name, gvcf_path)

        cram_job_depend_list = []
        for chromosome in CHROMOSOMES:
            output_cram_path = f"{MY_BUCKET}/contam_free/crams/cram_by_chrom/{sample_id}/{sample_id}_contam_free_{chromosome}.cram"
            # logger.info(f'Generating contamination free cram: {sample_id}-{chromosome}...')
            if not hl.hadoop_exists(output_cram_path):
                j_cram = b.new_python_job(
                    name=f"Run_contam_free_file_{sample_id}_{chromosome}",
                    attributes={"sample_id": sample_id, "chromosome": chromosome},
                )
                j_cram.depends_on(j_gvcf)
                j_cram.storage("500Gi")
                j_cram.call(
                    write_contam_free_cram_file,
                    gvcf_dict,
                    input_cram_file,
                    input_ref_fasta,
                    j_cram.ofile,
                    chromosome,
                )
                cram_job_depend_list.append(j_cram)
                contam_free_crams[sample_id][chromosome] = j_cram.ofile

                b.write_output(j_cram.ofile, output_cram_path)
            else:
                tmp_cram = b.read_input(output_cram_path)
                contam_free_crams[sample_id][chromosome] = tmp_cram

        cram_file_name = cram_files[s][0].split("/")[-1][:-5]
        output_full_cram_path = (
            f"{MY_BUCKET}/contam_free/crams/{cram_file_name}_contam_free.cram"
        )
        full_contam_free_crams.append(
            (output_full_cram_path, output_full_cram_path + ".crai")
        )
        j_cat_depend_on = None
        if not hl.hadoop_exists(output_full_cram_path):
            # with hl.hadoop_open(output_full_cram_path, 'w') as f:
            #     f.write(cram_file_name)
            logger.info(f"Concatenating crams: {sample_id}...")
            j_cat = b.new_job(
                name=f"Concatenate_contam_free_file_{sample_id}",
                attributes={"sample_id": sample_id},
            )
            if len(cram_job_depend_list) > 0:
                j_cat.depends_on(*cram_job_depend_list)
            j_cat.image(SAMTOOLS_IMAGE).storage("60Gi").memory("10Gi")
            tmp_cram_lst = reduce(
                lambda x, y: x + " " + y, contam_free_crams[sample_id].values()
            )
            j_cat.command(
                f"samtools view -H {input_cram_file['cram']} > {j_cat.header}"
            )
            j_cat.command(
                f"samtools cat -h {j_cat.header} -o {j_cat.ofile1} {tmp_cram_lst}"
            )
            j_cat.command(f"samtools index {j_cat.ofile1} -o {j_cat.ofile2}")
            j_cat_depend_on = j_cat

            b.write_output(j_cat.ofile1, output_full_cram_path)
            b.write_output(j_cat.ofile2, f"{output_full_cram_path}.crai")

        if hl.hadoop_exists(output_full_cram_path) and not hl.hadoop_exists(
            f"{MY_BUCKET}/contam_free/verifybamid/{sample_id}.selfSM"
        ):
            logger.info(f"Running VerifyBam ID: {sample_id}...")
            run_verifybamid(
                b=b,
                input_cram_path=output_full_cram_path,
                input_crai_path=f"{output_full_cram_path}.crai",
                cram_project_id="broad-mpg-gnomad",
                ref_fasta_path=REF_FASTA,
                contamination_sites_path=CONTAM_SITES,
                output_path=f"{MY_BUCKET}/contam_free/verifybamid/",
                output_prefix=sample_id,
                depend_on=j_cat_depend_on,
                disable_sanity_check=args.disable_sanity_check,
            )
        if i == 3:
            break  # test purpose

    if args.run_merge_gvcf_files:
        var_call_pipeline(
            input_bams=full_contam_free_crams,
            ref_fasta=REF_FASTA,
            ref_ind=REF_FASTA_INDEX,
            ref_dict=REF_DICT,
            calling_int_list=CALLING_INTERVAL_LIST,
            evaluation_int_list=EVALUATION_INTERVAL_LIST,
            dbsnp_vcf=DBSNP_VCF,
            dbsnp_vcf_ind=DBSNP_VCF_INDEX,
            contamination=0.0,
            scatter_count=50,
            out_dir=f"{MY_BUCKET}/contam_free",
            backend=backend,
        )

    if args.run_freemix_sum_table:
        contam_est = []
        for s in samples:
            contam_est.append(
                check_contam(1, f"{MY_BUCKET}/contam_free/verifybamid/", s)
            )
        sample_pairs = hl.Table.from_pandas(
            pd.DataFrame(list(zip(samples, contam_est)), columns=["sample_id", "freemix_score"])
        )
        sample_table = hl.Table.from_pandas(sample_pairs)
        sample_table.write(f"{MY_BUCKET}/hgdp_selected_sample_freemix_score.ht")

    logger.info(f"-------Mixing samples-------")
    logger.info(f"Preparing sample pairs...")
    MAIN = []
    CONTAM = []
    contam_samples_left = sample_pops
    for pop in POPs:
        main_samples = list(
            {key for (key, value) in sample_pops.items() if value == pop}
        )
        MAIN = MAIN + main_samples
        rest_pops = [p for p in POPs if p != pop]
        contam_samples = []
        for pop in rest_pops:
            sub_pop = {
                key for (key, value) in contam_samples_left.items() if value == pop
            }
            s2 = (
                list(sub_pop)[random.randint(0, len(sub_pop) - 1)]
                if len(sub_pop) > 0
                else sub_pop
            )
            contam_samples.append(s2)
            contam_samples_left = {
                key: value for (key, value) in contam_samples_left.items() if key != s2
            }
        CONTAM = CONTAM + contam_samples
    sample_pairs = hl.Table.from_pandas(
        pd.DataFrame(list(zip(MAIN, CONTAM)), columns=["original", "contaminant"])
    )
    sample_pairs.write(f"{MY_BUCKET}/mixed_samples/sample_pairs.ht", overwrite=True)

    mixed_crams = {}
    mixed_gvcfs = {}
    mixed_labels = []
    for contam_rate in CONTAM_RATES:
        for i in range(len(MAIN)):
            # s1 = MAIN[i]
            # s2 = CONTAM[i]
            s1 = "HGDP00021"
            s2 = "HGDP00105"
            OUT_BUCKET = f"{MY_BUCKET}/mixing_samples/crams/cram_by_chrom/contam_rate_{int(contam_rate*100)}/{s1}_{s2}"
            input_main_cram_file = b.read_input_group(
                cram=cram_files[s1][0], index=cram_files[s1][1]
            )
            input_contam_cram_file = b.read_input_group(
                cram=cram_files[s2][0], index=cram_files[s2][1]
            )
            mixing_samples_label = f"{s1}_{s2}_{int(contam_rate*100)}"
            mixed_labels.append(mixing_samples_label)
            mixed_crams[mixing_samples_label] = {}
            mixed_gvcfs[mixing_samples_label] = {}
            logger.info(
                f"Mixing contamination free crams: {mixing_samples_label}_percent_contamination..."
            )
            mix_cram_job_depend_list = []
            mix_gvcf_job_depend_list = []
            for chromosome in CHROMOSOMES:
                chromosome = "chr21"
                output_mix_cram_path = (
                    f"{OUT_BUCKET}/{mixing_samples_label}_{chromosome}.cram"
                )
                j_mix = b.new_python_job(
                    name=f"Mix_{mixing_samples_label}",
                    attributes={
                        "main": s1,
                        "contaminant": s2,
                        "contam_rate": f"{int(contam_rate*100)}\%",
                        "chromosome": chromosome,
                    },
                )
                j_mix.storage("500Gi").memory('50Gi')
                if len(cram_job_depend_list) > 0:
                    j_mix.depends_on(*cram_job_depend_list)
                j_mix.call(
                    mixing_two_crams,
                    s1,
                    s2,
                    input_main_cram_file["cram"],
                    input_contam_cram_file["cram"],
                    j_mix.ofile,
                    contam_rate,
                    chromosome,
                )

                mix_cram_job_depend_list.append(j_mix)
                mixed_crams[mixing_samples_label][chromosome] = j_mix.ofile
                b.write_output(j_mix.ofile, output_mix_cram_path)
                break

            output_mix_full_cram_path = (
                f"{MY_BUCKET}/mixed_samples/crams/{mixing_samples_label}_mixed.cram"
            )
            full_mix_crams = []
            j_cat_depend_on = None
            if not hl.hadoop_exists(output_mix_full_cram_path):
                logger.info(f"Concatenating mixed crams: {mixing_samples_label}...")
                j_cat = b.new_job(
                    name=f"Concatenate_mixed_file_{mixing_samples_label}",
                    attributes={
                        "main": s1,
                        "contaminant": s2,
                        "contam_rate": f"{int(contam_rate*100)}\%",
                    },
                )
                if len(mix_cram_job_depend_list) > 0:
                    j_cat.depends_on(*mix_cram_job_depend_list)
                j_cat.image(SAMTOOLS_IMAGE).storage("60Gi").memory("10Gi")
                tmp_cram_lst = reduce(
                    lambda x, y: x + " " + y, mixed_crams[mixing_samples_label].values()
                )
                j_cat.command(
                    f"samtools view -H {input_main_cram_file['cram']} > {j_cat.header}"
                )
                j_cat.command(
                    f"samtools cat -h {j_cat.header} -o {j_cat.ofile1} {tmp_cram_lst}"
                )
                j_cat.command(f"samtools index {j_cat.ofile1} -o {j_cat.ofile2}")
                b.write_output(j_cat.ofile1, output_mix_full_cram_path)
                b.write_output(j_cat.ofile2, f"{output_mix_full_cram_path}.crai")
                j_cat_depend_on = j_cat
                full_mix_crams.append(output_mix_full_cram_path)

            if hl.hadoop_exists(output_mix_full_cram_path) and not hl.hadoop_exists(
                f"{MY_BUCKET}/mixed_samples/verifybamid/{mixing_samples_label}.selfSM"
            ):
                logger.info(f"Running VerifyBam ID: {mixing_samples_label}...")
                run_verifybamid(
                    b=b,
                    input_cram_path=output_mix_full_cram_path,
                    input_crai_path=f"{output_mix_full_cram_path}.crai",
                    cram_project_id="broad-mpg-gnomad",
                    ref_fasta_path=REF_FASTA,
                    contamination_sites_path=CONTAM_SITES,
                    output_path=f"{MY_BUCKET}/mixed_samples/verifybamid/",
                    output_prefix=mixing_samples_label,
                    depend_on=j_cat_depend_on,
                    disable_sanity_check=args.disable_sanity_check,
                )
            break
        if args.run_merge_gvcf_files:
            var_call_pipeline(
                input_bams=full_mix_crams,
                ref_fasta=REF_FASTA,
                ref_ind=REF_FASTA_INDEX,
                ref_dict=REF_DICT,
                calling_int_list=CALLING_INTERVAL_LIST,
                evaluation_int_list=EVALUATION_INTERVAL_LIST,
                dbsnp_vcf=DBSNP_VCF,
                dbsnp_vcf_ind=DBSNP_VCF_INDEX,
                contamination=0.0,
                scatter_count=50,
                out_dir=f"{MY_BUCKET}/contam_free",
                backend=backend,
            )

        if args.run_freemix_sum_table:
            contam_est = []
            split_labels = pd.DataFrame([i.split("_") for i in mixed_labels], columns = ['original', 'contaminant', 'contam_rate'])
            for s in mixed_labels:
                contam_est.append(
                    check_contam(1, f"{MY_BUCKET}/mixed_samples/verifybamid/", s)
                )
            split_labels['freemix_score'] = contam_est
            mixed_table = hl.Table.from_pandas(split_labels)
            mixed_table.write(f"{MY_BUCKET}/hgdp_mixed_sample_freemix_score.ht")

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
        "--skip-gvcf-dict",
        help="Whether to skip running gvcf dict",
        action="store_true",
    )
    parser.add_argument(
        "--run-merge-gvcf-files",
        help="Whether to run freemix summary table",
        action="store_true",
    )
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
    args = parser.parse_args()
    print(args)
    main()


# python3 pipeline.py --input-sample-info gs://gnomad-wenhan/charr_simulation/hgdp_sample_path_info.csv --selected-samples gs://gnomad-wenhan/charr_simulation/hgdp_selected_sample_id.csv --skip-gvcf-dict
