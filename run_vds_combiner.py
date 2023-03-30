import hail as hl
import argparse
from simulation_utils.generics import *
from variant_calling.get_file_size import bytes_to_gb
from variant_calling.merge_gvcfs import merge_vcf


def main():
    if args.run_merge_gvcfs:
        backend = hb.ServiceBackend(
            billing_project=BILLING_PROJECT,
            remote_tmpdir=TMP_DIR,
        )
        b = hb.Batch(
            name=f"Merging-Files",
            default_python_image='us-central1-docker.pkg.dev/broad-mpg-gnomad/wlu/hail/hail-pysam-samtools:latest',
            backend=backend,
        )

        sample_ids = pd.read_csv(args.selected_samples, sep=',')
        for i in range(len(sample_ids)):
            gvcfs_to_merge = hl.utils.hadoop_ls(f'{args.input_gvcf_bucket}/{sample_ids[i]}/variant-calling/*.vcf')
            gvcfs_list = []
            gvcfs_sizes_sum = 0
            for file in gvcfs_to_merge:
                gvcfs_list.append(file['path'])
                gvcfs_sizes_sum += bytes_to_gb(file['path'])

            merge_disk_size = round(gvcfs_sizes_sum * 2.5) + 10
            merge_vcf(b=b,
                      gvcf_list=gvcfs_list,
                      storage=merge_disk_size,
                      output_vcf_name=f'{sample_ids[i]}_contam_free',
                      out_dir=f'{args.input_gvcf_bucket}')

    gvcfs_to_combine = hl.utils.hadoop_ls(f'{args.input_gvcf_bucket}/merged-gvcfs/*.vcf.gz')
    gvcfs_list = []
    for file in gvcfs_to_combine:
        gvcfs_list.append(file['path'])

    combiner = hl.vds.new_combiner(
        output_path=f"{MY_BUCKET}/{args.output_vds_name}.vds",
        temp_path=TMP_DIR,
        gvcf_paths=gvcfs_list,
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
        type=str,
        nargs="?",
    )
    parser.add_argument(
        "--input-gvcf-bucket",
        help="Input gvcf bucket",
        type=str,
        nargs="?",
        default="gs://gnomad-wenhan/charr_simulation/contam_free/"
    )
    parser.add_argument(
        "--output-vds-name",
        help="Output VDS name",
        type=str,
        nargs="?",
    )
    args = parser.parse_args()
    print(args)
    main()
