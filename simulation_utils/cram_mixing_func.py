from .generics import *


def edit_read_group(read, rg_name):
    rg_ind = [i for i, tuple in enumerate(read.tags) if "RG" in tuple]
    if len(rg_ind) > 0:
        temp_tag = [list(ele) for ele in read.tags]
        temp_tag[rg_ind[0]][-1] = rg_name
        read.tags = [tuple(ele) for ele in temp_tag]
    return read


def mixing_two_crams(
    main_cram_id: str,
    contam_cram_id: str,
    main_cram: str,
    contam_cram: str,
    ref_fasta: str,
    output_cram: str,
    contam_rate: float,
    chromosome: str,
    main_rg: str = "ERR1349775",
):
    pipe = pipes.Template()
    pipe.append(
        f"samtools view -C -T {ref_fasta} -h -o {output_cram}  2> /dev/null",
        "-.",
    )
    f = pipe.open(
        f"{main_cram_id}_{contam_cram_id}_{contam_rate * 100}_percent.sam", "w"
    )
    n_contam_read = 0
    n_main_read = 0
    chrom_id = chromosome[3:]
    if chrom_id == "X":
        chrom_id = 22
    elif chrom_id == "Y":
        chrom_id = 23
    else:
        chrom_id = int(chrom_id) - 1

    with pysam.AlignmentFile(
        f"{main_cram}", "rc", reference_filename=ref_fasta
    ) as file1, pysam.AlignmentFile(
        f"{contam_cram}", "rc", reference_filename=ref_fasta
    ) as file2:
        read1 = next(file1.fetch(until_eof=True))
        rg_ind = [i for i, tuple in enumerate(read1.tags) if "RG" in tuple]
        main_rg = read1.tags[rg_ind[0]][-1]
        read2 = edit_read_group(next(file2.fetch(until_eof=True)), main_rg)
        next_read1 = False
        next_read2 = False
        for i in range(CHROM_LENGTHS[chromosome]):
            while True:
                if (i < read1.pos or next_read1) and (i < read2.pos or next_read2):
                    break
                if i == read1.pos and (not next_read1):
                    if np.random.binomial(1, 1 - contam_rate):
                        n_main_read += 1
                        f.write(read1.tostring() + "\n")
                    try:
                        read1 = edit_read_group(next(file1), main_rg)
                    except StopIteration:
                        next_read1 = True
                elif read1.reference_id != chrom_id:
                    next_read1 = True
                if i == read2.pos and (not next_read2):
                    if np.random.binomial(1, contam_rate):
                        n_contam_read += 1
                        f.write(read2.tostring() + "\n")
                    try:
                        read2 = edit_read_group(next(file2), main_rg)
                    except StopIteration:
                        next_read2 = True
                elif read2.reference_id != chrom_id:
                    next_read2 = True
        print(f"Number of reads inserted: {n_contam_read}")
        print(f"Number of original reads: {n_main_read}")
        file1.close()
        file2.close()
        f.close()
