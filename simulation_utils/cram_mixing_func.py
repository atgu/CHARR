from .generics import *

def edit_read_group(read, rg_name):
    if "RG" in read.tags[-1]:
        temp_tag = [list(ele) for ele in read.tags]
        temp_tag[-1][-1] = rg_name
        read.tags = [tuple(ele) for ele in temp_tag]
    return read

def mixing_two_crams(
        main_cram_id: str,
        contam_cram_id: str,
        main_cram_cram: str,
        contam_cram_cram: str,
        output_cram: str,
        contam_rate: float,
        chromosome: str,
        main_rg: str = "ERR1349775",
):
    pipe = pipes.Template()
    pipe.append(
        f"samtools view -C -T {REF_FASTA_PATH} -h -o {output_cram}",
        "-.",
    )
    f = pipe.open(f"{main_cram_id}_{contam_cram_id}_{int(contam_rate * 100)}_percent.sam", "w")
    n_contam_read = 0
    with pysam.AlignmentFile(
            f"{main_cram_cram}", "rc"
    ) as file1, pysam.AlignmentFile(f"{contam_cram_cram}", "rc") as file2:
        # Tests:
        read1 = edit_read_group(next(file1.fetch(until_eof=True)), main_rg)
        read2 = edit_read_group(next(file2.fetch(until_eof=True)), main_rg)
        print(chromosome)
        next_read1 = False
        next_read2 = False
        for i in range(CHROM_LENGTHS[chromosome[3]]):
            while True:
                if (i < read1.pos or next_read1) and (i < read2.pos or next_read2):
                    break
                if i == read1.pos and (not next_read1):
                    if np.random.binomial(1, 1 - contam_rate):
                        f.write(read1.tostring() + "\n")
                    try:
                        read1 = edit_read_group(next(file1), main_rg)
                    except StopIteration:
                        read1.pos = 1e100
                elif read1.reference_id != chrom:
                    next_read1 = True
                if i == read2.pos and (not next_read2):
                    if np.random.binomial(1, contam_rate):
                        n_contam_read += 1
                        f.write(read2.tostring() + "\n")
                    try:
                        read2 = edit_read_group(next(file2), main_rg)
                    except StopIteration:
                        read2.pos == 1e100
                elif read2.reference_id != chrom:
                    next_read2 = True
        print(f"Number of reads inserted: {n_contam_read}")
        file1.close()
        file2.close()
        f.close()