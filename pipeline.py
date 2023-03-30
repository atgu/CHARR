#!/usr/bin/env python3

__author__ = 'Wenhan Lu'


import argparse
from functools import reduce

from variant_calling.haplotype_caller import haplotype_caller_gatk

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
        requester_pays_project='daly-ibd',
        default_python_image='us-central1-docker.pkg.dev/broad-mpg-gnomad/wlu/hail/hail-pysam-samtools:latest',
        backend=backend,
    )


    logger.info('Loading sample info...')
    sample_paths = pd.read_csv(args.input_sample_info, sep=',')
    sample_ids = pd.read_csv(args.selected_samples, sep=',')
    sample_paths = sample_paths[sample_paths['s'].isin(sample_ids['hgdp_id'])]
    samples = []
    cram_files = {}
    gvcf_files = {}
    sample_pops = {}
    for sample in sample_paths.iterrows():
        samples.append(sample[1][0])
        cram_files[sample[1][0]] = (f'{sample[1][1]}', f'{sample[1][1]}.crai')
        gvcf_files[sample[1][0]] = (f'{sample[1][2]}', f'{sample[1][2]}.tbi')
        sample_pops[sample[1][0]]= sample[1][3]

    input_ref_fasta = b.read_input_group(
        **{"fasta": REF_FASTA_PATH,
           "fasta.fai": REF_FASTA_INDEX,
           "dict":REF_DICT}
    )

    contam_free_gvcfs = {}
    contam_free_crams = {}
    for s in samples:
        sample_id = s
        contam_free_crams[sample_id] = {}
        contam_free_gvcfs[sample_id] = {}
        logger.info(f'-------Decontaminating {sample_id}-------')
        input_cram_file = b.read_input_group(
            **{"cram": cram_files[s][0], "cram.crai": cram_files[s][1]}
        )
        if not args.skip_gvcf_dict:
            logger.info(f'Generating gvcf dict - {sample_id}...')
            gvcf_path = gvcf_files[s]
            gvcf_file_name = gvcf_path[0].split('/')[-1][:-7]
            gvcf_dict = run_gvcf_dict(b, gvcf_file_name, gvcf_path)

        cram_job_depend_list = []
        gvcf_job_depend_list = []
        for chromosome in chromosomes:
            chromosome = 'chr21'
            output_cram_path = f"{MY_BUCKET}/contam_free/{sample_id}/cram_by_chrom/{sample_id}_contam_free_{chromosome}.cram"
            logger.info(f'Generating contamination free cram: {sample_id}-{chromosome}...')
            depend_on = None
            if not hl.hadoop_exists(output_cram_path):
                j_cram = b.new_python_job(name=f"Run_{sample_id}_{chromosome}_contam_free_file")
                j_cram._machine_type = "n1-highmem-16"
                j_cram.storage("500Gi")
                j_cram.call(
                    write_contam_free_cram_file,
                    gvcf_dict,
                    input_cram_file,
                    input_ref_fasta,
                    j_cram.ofile,
                    chromosome
                )

                logger.info(f'Indexing contamination free cram: {sample_id}-{chromosome}...')
                j_index = b.new_job(name=f"Index_{sample_id}_{chromosome}_contam_free_file")
                depend_on = j_index
                cram_job_depend_list.append(j_index)
                j_index.image('gcr.io/genomics-tools/samtools')
                j_index._machine_type = "n1-highmem-16"
                j_index.storage("50Gi")
                j_index.command(f"samtools reheader -i {input_cram_file['cram']} {j_cram.ofile} > {j_index.ofile1}")
                j_index.command(f"samtools index {j_cram.ofile} -o {j_index.ofile2}")
                contam_free_crams[sample_id][chromosome] = j_index.ofile1

                b.write_output(j_index.ofile1, output_cram_path)
                b.write_output(j_index.ofile2, f"{output_cram_path}.crai")

            elif not hl.hadoop_exists(f'{output_cram_path}.crai'):
                tmp_cram = b.read_input(output_cram_path)
                logger.info(f'Indexing contamination free cram: {sample_id}-{chromosome}...')
                j_index = b.new_job(name=f"Index_{sample_id}_{chromosome}_contam_free_file")
                depend_on = j_index
                cram_job_depend_list.append(j_index)
                j_index.image('gcr.io/genomics-tools/samtools')
                j_index._machine_type = "n1-highmem-16"
                j_index.storage("50Gi")
                # j_index.command(f"samtools reheader -i {input_cram_file['cram']} {tmp_cram} > {j_index.ofile1}")
                # j_index.command(f"samtools index {j_index.ofile1} -o {j_index.ofile2}")
                j_index.command(f"samtools index {tmp_cram} -o {j_index.ofile2}")
                contam_free_crams[sample_id][chromosome] = j_index.ofile1

                b.write_output(j_index.ofile1, output_cram_path)
                # b.write_output(j_index.ofile2, f"{output_cram_path}.crai")

            else:
                tmp_cram = b.read_input(output_cram_path)
                contam_free_crams[sample_id][chromosome] = tmp_cram


            if not hl.hadoop_exists(f'{MY_BUCKET}/contam_free/{sample_id}/variant-calling/{sample_id}_{chromosome}.g.vcf'):
                logger.info(f'Running haplotype caller: {sample_id}-{chromosome}...')
                tmp_gvcf, tmp_job = haplotype_caller_gatk(b=b,
                                                 depend_on=depend_on,
                                                 input_bam=contam_free_crams[sample_id][chromosome],
                                                 ref_fasta=input_ref_fasta,
                                                 out_dir=f'{MY_BUCKET}/contam_free/{sample_id}',
                                                 contamination=0.0,
                                                 bam_filename_no_ext=f'{sample_id}_{chromosome}',
                                                 storage=40,
                                                 interval_list_name=None)
                gvcf_job_depend_list.append(tmp_job)
                contam_free_gvcfs[sample_id][chromosome] = tmp_gvcf
            break  # test

        logger.info(f'Concatenating crams: {sample_id}...')
        j_cat = b.new_job(name=f"Concatenate_contam_free_file_{sample_id}")
        if len(cram_job_depend_list) > 0:
            j_cat.depends_on(*cram_job_depend_list)
        j_cat.image('gcr.io/genomics-tools/samtools')
        tmp_cram_lst = reduce(lambda x,y: x+' '+y, contam_free_crams[sample_id].values())
        print(tmp_cram_lst)
        j_cat.command(f'samtools cat -o {j_cat.ofile1} {tmp_cram_lst}')
        j_cat.command(f'samtools index {j_cat.ofile1} -o {j_cat.ofile2}')
        cram_file_name = cram_files[s][0].split('/')[-1][:-5]
        b.write_output(j_cat.ofile1, f"{MY_BUCKET}/contam_free/crams/{cram_file_name}_contam_free.cram")
        b.write_output(j_cat.ofile2, f"{MY_BUCKET}/contam_free/crams/{cram_file_name}_contam_free.cram.crai")

        logger.info(f'Running VerifyBam ID: {sample_id}...')
        run_verifybamid(
            b=b,
            input_cram_path=j_cat.ofile1,
            input_crai_path=j_cat.ofile2,
            ref_fasta_path=REF_FASTA,
            contamination_sites_path=CONTAM_SITES,
            output_path=f"{MY_BUCKET}/contam_free/verifybamid/",
            output_prefix=sample_id,
            disable_sanity_check=args.disable_sanity_check,
        )

        # gvcf_path = gvcf_files[s]
        # gvcf_file_name = gvcf_path[0].split('/')[-1][:-7]
        # merge_vcf(b=b,
        #           gvcf_list=list(contam_free_gvcfs[s].values()),
        #           depend_on=gvcf_job_depend_list,
        #           storage='50Gi',
        #           output_vcf_name=f'{gvcf_file_name}_contam_free',
        #           out_dir=f'{MY_BUCKET}/contam_free')
        break # test purpose

    logger.info(f'-------Mixing samples-------')
    logger.info(f'Preparing sample pairs...')
    MAIN = []
    CONTAM = []
    contam_samples_left = sample_pops
    for pop in POPs:
        main_samples = list({key for (key, value) in sample_pops.items() if value == pop})
        MAIN = MAIN + main_samples
        rest_pops = [p for p in POPs if p != pop]
        contam_samples = []
        for pop in rest_pops:
            sub_pop = {key for (key, value) in contam_samples_left.items() if value == pop}
            s2 = list(sub_pop)[random.randint(0, len(sub_pop) - 1)] if len(sub_pop) > 0 else sub_pop
            contam_samples.append(s2)
            contam_samples_left = {key: value for (key, value) in contam_samples_left.items() if key != s2}
        CONTAM = CONTAM + contam_samples

    mixed_crams = {}
    mixed_gvcfs = {}
    mixed_labels = []
    for contam_rate in CONTAM_RATES:
        for i in range(len(MAIN)):
            s1 = MAIN[i]
            s2 = CONTAM[i]
            OUT_BUCKET = f'{MY_BUCKET}/mixing_samples/contam_rate_{int(contam_rate*100)}/{s1}_{s2}'
            input_main_cram_file = b.read_input_group(
                **{"cram": cram_files[s1][0], "cram.crai": cram_files[s1][1]}
            )
            input_contam_cram_file = b.read_input_group(
                **{"cram": cram_files[s2][0], "cram.crai": cram_files[s2][1]}
            )
            mixing_samples_label = f'{s1}_{s2}_{int(contam_rate*100)}'
            mixed_labels.append(mixing_samples_label)
            mixed_crams[mixing_samples_label] = {}
            mixed_gvcfs[mixing_samples_label] = {}
            logger.info(f'Mixing contamination free crams: {mixing_samples_label}_percent_contamination...')
            mix_cram_job_depend_list = []
            mix_gvcf_job_depend_list = []
            for chromosome in chromosomes:
                chromosome = 'chr21'
                output_mix_cram_path = f"{OUT_BUCKET}/cram_by_chrom/{mixing_samples_label}_{chromosome}.cram"
                j_mix = b.new_python_job(name=f'Mix_{mixing_samples_label}')
                j_mix.depends_on(*cram_job_depend_list)
                j_mix.call(
                    mixing_two_crams,
                    s1,
                    s2,
                    input_main_cram_file['cram'],
                    input_contam_cram_file['cram'],
                    j_mix.ofile,
                    contam_rate,
                    chromosome
                )

                logger.info(f'Indexing mixed contam free cram: {mixing_samples_label}-{chromosome}...')
                j_index = b.new_job(name=f"Index_{mixing_samples_label}_{chromosome}_mixed_file")
                depend_on = j_index
                mix_cram_job_depend_list.append(j_index)
                j_index.image('gcr.io/genomics-tools/samtools')
                j_index._machine_type = "n1-highmem-16"
                j_index.storage("50Gi")
                j_index.command(f"samtools reheader -i {input_main_cram_file['cram']} {j_mix.ofile} > {j_index.ofile1}")
                j_index.command(f"samtools index {j_mix.ofile} -o {j_index.ofile2}")
                mixed_crams[mixing_samples_label][chromosome] = j_index.ofile1

                b.write_output(j_index.ofile1, output_mix_cram_path)
                b.write_output(j_index.ofile2, f"{output_mix_cram_path}.crai")

                if not hl.hadoop_exists(f'{OUT_BUCKET}/variant-calling/{mixing_samples_label}_{chromosome}.g.vcf'):
                    logger.info(f'Running haplotype caller on mixed cram: {mixing_samples_label}_{chromosome}...')
                    tmp_gvcf, tmp_job = haplotype_caller_gatk(b=b,
                                                              depend_on=depend_on,
                                                              input_bam=mixed_crams[mixing_samples_label][chromosome],
                                                              ref_fasta=input_ref_fasta,
                                                              out_dir=f'{OUT_BUCKET}',
                                                              contamination=0.0,
                                                              bam_filename_no_ext=f'{mixing_samples_label}_{chromosome}',
                                                              storage=40,
                                                              interval_list_name=None)
                    mix_gvcf_job_depend_list.append(tmp_job)
                    mixed_gvcfs[mixing_samples_label][chromosome] = tmp_gvcf
                break

            logger.info(f'Concatenating mixed crams: {mixing_samples_label}...')
            j_cat = b.new_job(name=f"Concatenate_mixed_file_{mixing_samples_label}")
            if len(mix_cram_job_depend_list) > 0:
                j_cat.depends_on(*mix_cram_job_depend_list)
            j_cat.image('gcr.io/genomics-tools/samtools')
            tmp_cram_lst = reduce(lambda x, y: x + ' ' + y, mixed_crams[mixing_samples_label].values())
            print(tmp_cram_lst)
            j_cat.command(f'samtools cat -o {j_cat.ofile1} {tmp_cram_lst}')
            j_cat.command(f'samtools index {j_cat.ofile1} -o {j_cat.ofile2}')
            b.write_output(j_cat.ofile1, f"{MY_BUCKET}/mixed_samples/crams/{mixing_samples_label}_mixed.cram")
            b.write_output(j_cat.ofile2, f"{MY_BUCKET}/mixed_samples/crams/{mixing_samples_label}_mixed.cram.crai")

            logger.info(f'Running VerifyBam ID: {mixing_samples_label}...')
            run_verifybamid(
                b=b,
                input_cram_path=j_cat.ofile1,
                input_crai_path=j_cat.ofile2,
                ref_fasta_path=REF_FASTA,
                contamination_sites_path=CONTAM_SITES,
                output_path=f"{MY_BUCKET}/mixed_samples/verifybamid/",
                output_prefix=mixing_samples_label,
                disable_sanity_check=args.disable_sanity_check,
            )

            # logger.info(f'Merging gvcfs: {mixing_samples_label}...')
            # merge_vcf(b=b,
            #           gvcf_list=list(mixed_gvcfs[mixing_samples_label].values()),
            #           depend_on=mix_gvcf_job_depend_list,
            #           storage='50Gi',
            #           output_vcf_name=f'{mixing_samples_label}_mixed',
            #           out_dir=f'{MY_BUCKET}/mixed_samples')
            break
        break

    b.run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-sample-info",
        help="Path to a csv of paths to cram and gvcf files with sample IDs",
        type=str,
        nargs="?",
    )
    parser.add_argument(
        "--selected-samples",
        help="Path to a list of sample IDs to run and mix",
        type=str,
        nargs="?",
    )
    parser.add_argument('--skip-gvcf-dict', help='Whether to skip running gvcf dict', action='store_true')
    parser.add_argument('--run-merge-gvcfs', help='Whether to run merging gvcfs', action='store_true')
    parser.add_argument('--disable-sanity-check', help='Whether to use sanity check in verifybamID', action="store_true")
    args = parser.parse_args()
    print(args)
    main()



# python3 pipeline.py --input-sample-info

