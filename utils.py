import hail as hl
import logging
from typing import Optional, Union
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.reference_genome import get_reference_genome

gnomad_ref_af = {
    "genomes": "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht",
    "exomes": "gs://gcp-public-data--gnomad/release/2.1.1/ht/exomes/gnomad.exomes.r2.1.1.sites.ht",
}

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("CHARR")
logger.setLevel(logging.INFO)


def add_liftover_rg37_to_rg38_ht(
    t: Union[hl.MatrixTable, hl.Table],
):
    """
    Add liftover to Hail Table or MatrixTable from rg37 to rg38
    :param t: Hail Table or MatrixTable to add liftover on
    :return: Hail Table or MatrixTable
    """
    rg37 = hl.get_reference("GRCh37")
    rg38 = hl.get_reference("GRCh38")
    if not rg37.has_liftover("GRCh38"):
        rg37.add_liftover(
            "gs://hail-common/references/grch37_to_grch38.over.chain.gz", rg38
        )
    is_mt = isinstance(t, hl.MatrixTable)
    if is_mt:
        t = t.annotate_rows(new_locus=hl.liftover(t.locus, "GRCh38"))
        t = t.filter_rows(hl.is_defined(t.new_locus))
        t = t.key_rows_by(locus=t.new_locus, alleles=t.alleles)
    else:
        t = t.annotate(new_locus=hl.liftover(t.locus, "GRCh38"))
        t = t.filter(hl.is_defined(t.new_locus))
        t = t.key_by(locus=t.new_locus, alleles=t.alleles)
    return t


def import_genotype_data(format: str, path: str, num_partitions: Optional[int] = 1000):
    """
    Imports genotype data into a MatrixTable or VariantDataset.
    :param format: Format of the original data, select from ['vds', 'mt', 'vcf', 'gvcf'].
    :param path: Path to input genotype information.
    :param num_partitions: Number of partitions to use for importing vcf file
    :return: Hail MatrixTable or VariantDataset
    """
    if format not in ("vds", "mt", "vcf", "gvcf"):
        raise Exception("Select format as one of ['vds', 'mt', 'vcf', 'gvcf']")
    if format == "vds":
        mtds = hl.vds.read_vds(path)
    elif format == "mt":
        mtds = hl.read_matrix_table(path)
    elif format == "vcf":
        mtds = hl.import_vcf(
            path,
            force_bgz=True,
            reference_genome="GRCh38",
        ).repartition(num_partitions)
    elif format == "gvcf":
        mtds = hl.import_gvcfs(path, reference_genome="GRCh38")
        mtds = hl.MatrixTable(
            hl.ir.MatrixKeyRowsBy(
                mtds._mir, ["locus", "alleles"], is_sorted=True
            )  # Prevents hail from running sort on genotype MT which is already sorted by a unique locus
        )
    return mtds


def compute_charr(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    min_af: float = 0.05,
    max_af: float = 0.95,
    min_dp: int = 10,
    max_dp: int = 100,
    min_gq: int = 20,
    ref_AF_field: Optional[str] = None,
    data_type: Optional[str] = "genomes",
):
    """
    Compute CHARR, the DNA sample contamination estimator, from a Hail Sparse Matrixtable or VariantDataset.

    :param mtds: Input MatrixTable or VariantDataset.
    :param min_af: Minimum reference allele frequency to filter variants. Default: 5%.
    :param max_af: Maximum reference allele frequency to filter variants. Default: 95%.
    :param min_dp: Minimum sequencing depth to filter variants. Default: 20.
    :param max_dp: Maximum sequencing depth to filter variants. Default: 100.
    :param min_gq: Minimum genotype quality to filter variants. Default: 20.
    :param ref_AF_field: Name of the reference allele frequency (ref-AF) field. If not set, then ref-AF is either estimated from the data when the sample size is larger than 10,000 or obtained from gnomAD.
    :param data_type: Type of callset, select from ['genomes', 'exomes']. This is only used when ref_AF_field is not applicable and the sample size of the data is below 10,000.
    :return: Hail MatrixTable with CHARR annotated in the column field.
    """

    # Determine whether the input data is in the VDS format; if not, convert matrixtable to VDS and extract only the variant call information
    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if is_vds:
        mt = mtds.variant_data
    else:
        logger.info("Converting MatrixTable to VariantDataset...")
        mtds = mtds.select_entries("END", "LA", "LGT", "GQ", "DP", "LAD")
        vds = hl.vds.VariantDataset.from_merged_representation(mtds)
        mt = vds.variant_data

    # Annotate reference allele frequency when it is not defined in the original data, and name it 'ref_AF'.
    if ref_AF_field is None:
        n_samples = mt.count_cols()
        if n_samples < 10000:
            if data_type not in ("exomes", "genomes"):
                raise Exception(
                    "Select data_type as one of 'genomes' or 'exomes' to use the correct reference allele frequency information in gnomAD"
                )
            af_ht = hl.read_table(gnomad_ref_af[data_type])
            af_ht = af_ht.annotate(ref_AF=1 - af_ht.freq[1].AF)
            build_mt = get_reference_genome(mt.locus).name
            build_ht = get_reference_genome(af_ht.locus).name
            # Liftover to GRCh38 if two tables are using different builds
            if (build_mt != build_ht) & (build_mt == "GRCh37"):
                mt = add_liftover_rg37_to_rg38_ht(mt)
            elif (build_mt != build_ht) & (build_ht == "GRCh37"):
                af_ht = add_liftover_rg37_to_rg38_ht(af_ht)
            mt = mt.annotate_rows(ref_AF=af_ht[mt.row_key].ref_AF)

        else:
            n_alleles = 2 * n_samples
            mt = mt.annotate_rows(
                ref_AF=1 - hl.agg.sum(mt.LGT.n_alt_alleles()) / n_alleles
            )
        ref_AF_field = "ref_AF"

    # Filter to autosomal biallelic SNVs with reference allele frequency within the range (min_af, max_af)
    mt = filter_to_autosomes(mt)
    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2)
        & hl.is_snp(mt.alleles[0], mt.alleles[1])
        & (mt[ref_AF_field] > min_af)
        & (mt[ref_AF_field] < max_af)
    )

    # Filter to variant calls with GQ above min_gq and DP within the range (min_dp, max_dp)
    mt = mt.filter_entries(
        mt.LGT.is_hom_var() & (mt.GQ > min_gq) & (mt.DP > min_dp) & (mt.DP < max_dp)
    )

    # Compute CHARR
    mt = mt.annotate_cols(
        charr=hl.agg.mean((mt.LAD[0] / (mt.LAD[0] + mt.LAD[1])) / mt[ref_AF_field])
    )

    mt = mt.annotate_globals(
        af_min=min_af,
        af_max=max_af,
        dp_min=min_dp,
        dp_max=max_dp,
        gq_min=min_gq,
    )

    return mt


def write_col_table(
    mt: hl.MatrixTable, output_dir: str, output_name: str, extension: str = "ht", overwrite: bool = False
):
    """
    Write the column table of a Hail MatrixTable.

    :param mt: Input MatrixTable.
    :param output_dir: Directory to write the table to.
    :param output_name: Filename to use to write the table.
    :param extension: Extension of the output file, select from ['ht', 'tsv'].
    :param overwrite: If ``True``, overwrite an existing Hail Table at the destination.
    :return: None.
    """

    if extension not in ("ht", "tsv"):
        raise Exception("Select extension as one of 'ht' or 'tsv'")
    ht = mt.cols()
    if extension == "ht":
        ht.write(f"{output_dir}/{output_name}.ht", overwrite=overwrite)
    elif extension == "tsv":
        ht.export(f"{output_dir}/{output_name}.tsv")

def run_charr(
        format: str, path: str,
        output_dir: Optional[str],
        output_name: Optional[str],
        min_af: float = 0.05,
        max_af: float = 0.95,
        min_dp: int = 10,
        max_dp: int = 100,
        min_gq: int = 20,
        ref_AF_field: Optional[str] = None,
        data_type: Optional[str] = "genomes",
        num_partitions: Optional[int] = 1000,
        write_charr: bool=True,
        extension: Optional[str] = "ht",
        overwrite: bool = False
):
    mtds=import_genotype_data(format, path, num_partitions)
    mt = compute_charr(mtds, min_af, max_af, min_dp, max_dp, min_gq, ref_AF_field, data_type)
    if write_charr:
        write_col_table(mt, output_dir, output_name, extension, overwrite)
    return mt