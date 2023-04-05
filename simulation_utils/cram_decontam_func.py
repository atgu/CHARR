
from pysam import VariantFile
from typing import Tuple
import os
import logging
import pickle
import hailtop.batch as hb
import hail as hl
from .generics import *

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("CHARR simulation pipeline")
logger.setLevel(logging.INFO)


def sample_variant_site(genotype, error, ref, alt1, alt2):
    p_err = np.random.binomial(1, error)
    if genotype == (1, 1):
        allele = random.sample(ALLELES, 1) if p_err else alt1
    elif genotype == (0, 1):
        allele = random.sample(ALLELES, 1) if p_err else random.sample((alt1, ref), 1)
    elif genotype == (1, 2):
        allele = random.sample(ALLELES, 1) if p_err else random.sample((alt1, alt2), 1)
    elif genotype == (2, 2):
        allele = random.sample(ALLELES, 1) if p_err else alt2
    else:
        raise Exception("Expecting non reference genotype")
    return allele


def get_var_dict_from_gvcf(input_gvcf: str):
    # gvcf -> nested dictionary: {chromosome : { position : infomration ...}}
    # os.system(f"ls {input_gvcf}")
    print(pysam.__version__)
    var_dict = {}
    gvcf_in = VariantFile(input_gvcf)  # auto-detect input format
    for rec in gvcf_in.fetch():
        if (rec.samples[0]["GT"][0] != 0) | (rec.samples[0]["GT"][1] != 0):
            if rec.contig not in var_dict.keys():
                var_dict[rec.contig] = {}
            var_dict[rec.contig][rec.pos] = (
                rec.ref,
                rec.alts[0],
                rec.alts[1],
                rec.samples[0]["GT"],
            )
    return var_dict


def write_var_dict_from_gvcf(output: str, gvcf_dict: dict):
    with open(output, "wb") as f:
        pickle.dump(gvcf_dict, f)
    f.close()


def read_pickle_dict(dict_path: str):
    dict = pd.read_pickle(dict_path)
    return dict


def run_gvcf_dict(b: hb.batch, gvcf_file_name: str, gvcf_path: Tuple):
    gvcf_dict_path = f"{MY_BUCKET}/gvcf_dicts/{gvcf_file_name}_gvcf_dict.pkl"
    sample_id = gvcf_file_name.split(".")[0]
    j = b.new_python_job(name=f"Run_gvcf_dict_{sample_id}")
    j._machine_type = "n1-highmem-16"
    j.storage("500Gi")
    if hl.hadoop_exists(gvcf_dict_path):
        input_gvcf_dict = b.read_input(gvcf_dict_path)
        gvcf_dict = j.call(read_pickle_dict, input_gvcf_dict)
    else:
        input_gvcf = b.read_input_group(gvcf=gvcf_path[0], index=gvcf_path[1])

        gvcf_dict = j.call(get_var_dict_from_gvcf, input_gvcf.gvcf)
        j.call(write_var_dict_from_gvcf, j.ofile, gvcf_dict)
        b.write_output(j.ofile, gvcf_dict_path)
    return gvcf_dict, j


def write_contam_free_cram_file(
    gvcf_dict: dict,
    input_cram_file: hb.resource.ResourceGroup,
    input_ref_fasta: hb.resource.ResourceGroup,
    output_cram_file: hb.ResourceFile,
    # output_crai_file: hb.ResourceFile,
    chromosome: str,
):
    print(input_cram_file)
    cram_in = pysam.AlignmentFile(input_cram_file["cram"], "rc")
    pipe = pipes.Template()
    pipe.append(
        f"samtools view -C -T {input_ref_fasta['fasta']} -h -o {output_cram_file}", "-."
    )
    f = pipe.open(f"contam_free.sam", "w")
    ref_fasta = pysam.FastaFile(input_ref_fasta["fasta"])
    j = 0
    for read in cram_in.fetch(chromosome):
        j+=1
        if read.reference_id < 22:
            chrom = f"chr{read.reference_id + 1}"
        elif read.reference_id == 22:
            chrom = "chrX"
        elif read.reference_id == 23:
            chrom = "chrY"
        start_pos = read.pos
        end_pos = start_pos + 151
        error = read.qual
        tags = read.tags
        read.query_sequence = ref_fasta.fetch(chrom, start_pos, end_pos)
        no_indel = True
        for i in range(151):
            var_pos = gvcf_dict[chrom]
            current_pos = start_pos + i + 1
            current_err = 10 ** (-ord(error[i]) / 10)
            if current_pos in var_pos:
                ref, alt1, alt2, genotype = var_pos[current_pos]
                alleles = (ref, alt1, alt2)
                if (genotype[0] is not None) & (genotype[1] is not None):
                    try:
                        if (len(alleles[genotype[0]]) > 1) or (
                            len(alleles[genotype[1]]) > 1
                        ):
                            no_indel = False
                            break
                        else:
                            temp = list(read.query_sequence)
                            temp[i] = sample_variant_site(
                                genotype, current_err, ref, alt1, alt2
                            )[0]
                            read.query_sequence = "".join(temp)
                    except TypeError:
                        print(f"Read info: \n {read}")
                        print(f"Var info: \n {gvcf_dict[chrom][current_pos]}")
        # If a read overlaps with an indel sites don't write the read out
        if no_indel:
            read.qual = error
            read.tags = tags
            if j < 10:
                print(read)
                print(read.tostring())
            f.write(read.tostring() + "\n")
    cram_in.close()
    f.close()
