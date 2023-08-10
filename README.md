# CHARR
Code used in the analysis for CHARR (Contamination from Homozygous Alternate Reference Reads), including early-stage explorations, simulations, freemix score recomputation, and producing figures.

## Project Overview
### Description
CHARR, Contamination from Homozygous Alternate Reference Reads, a contamination estimator which leverages the infiltration of reference reads within homozygous alternate variant calls. CHARR uses a small proportion of variant-level genotype information and thus can be computed from single-sample gVCFs or callsets in VCF or BCF formats, as well as efficiently stored variant calls in Hail VDS format. Our results demonstrate that CHARR accurately recapitulates results from existing tools with substantially reduced costs, improving the accuracy and efficiency of downstream analyses of ultra-large whole genome and exome sequencing datasets. 

### Contributors
- Wenhan Lu ([@wlu04](https://github.com/wlu04))


### Reference
- Cite: https://www.biorxiv.org/content/10.1101/2023.06.28.545801v1

## Data Overview
We used data from the Genome Aggregation Database, including:
- 59,765 release whole genome samples in gnomAD v3 joint called in (Chen et al., 2022). 58,986 samples sequenced at the Broad Institute + 779 HGDP samples
- 102,063 release whole exome samples in gnomAD v2 joint called in (Karczewski et al., 2020). 103,027 samples from gnomAD v2, excluding 10 samples with fewer than 10,000 heterozygous variants and 954 samples with an old version of freemix used
- 948 HGDP samples in gnomAD v3 joint called in (Chen et al., 2022) and described in (Koenig et al., 2023)

## Methods Overview

```math
CHARR = \frac{1}{m}\sum_j\frac{RR_i}{p_j(RR_j + AR_j)}
```
$RR_j:$ Number of reference reads called for variant j \
$AR_j:$ Number of alternate reads called for variant j \
$p_j:$ Reference allele frequency of variant j \
$m:$ total number of high-quality homozygous alternate variants <br />
Default parameter configuration: autosomal, biallelic homozygous alternate SNVs with GQ $\geq$ 20, 100 $\geq$ DP $\geq$ 20 and 0.9 $\geq$ ref_AF $\geq$ 0.1.

## Implementations
We built two implementations of CHARR, one in Hail and one using the VCF format:
- Hail (Hail Team, 2023): https://hail.is/docs/0.2/methods/genetics.html#hail.methods.compute_charr 
  * The function computes CHARR on Hail MatrixTables and Variant Datasets (VDS), both of which can be generated by importing a VCF file or set of gVCFs
- SCE-VCF: https://github.com/HTGenomeAnalysisUnit/SCE-VCF

## Analyses

### Comparison between CHARR and VerifyBamID
* Run VerifyBamID:
  * Hail Batch pipeline: https://github.com/atgu/CHARR/blob/main/batch_verifybamid.py
* Run Customized CHARR:
  * Functions: https://github.com/atgu/CHARR/blob/main/utils.py
  * Pipeline: https://github.com/atgu/CHARR/blob/main/compute_charr.py


### Simulation Framework
In order to compare the accuracies of CHARR and VerifyBamID, we designed a pipeline to simulate potential contamination scenarios and manually introduce a series of known contamination rates. We randomly selected 30 samples from the HGDP dataset, which included 5 samples from each of the 6 genetic ancestry groups, with their original contamination rates approximately distributed uniformly within each group. <br />
For each sample, we apply the 3 steps below: 
1. **Decontamination**: Decontaminating their short-read data by incorporating information from their corresponding gVCF files and the reference genome.
   * Functions: https://github.com/atgu/CHARR/blob/main/simulation_utils/cram_decontam_func.py
2. **Two-way mixing simulation**: Randomly pairing the samples and introducing short reads from a contaminating sample to a target sample at a range of contamination rates.
   * Functions: https://github.com/atgu/CHARR/blob/main/simulation_utils/cram_mixing_func.py
3. **N-way mixing simulation**: Mixing reads from all decontaminated samples at a range of contamination rates.
- Simulation code:
  - Generic settings: https://github.com/atgu/CHARR/blob/main/simulation_utils/generics.py
  - Complete simulation pipeline on Hail Batch: https://github.com/atgu/CHARR/blob/main/pipeline.py
  - Combine the gVCF files called from the simulated data into a VDS file: https://github.com/atgu/CHARR/blob/main/run_vds_combiner.py


### Figures
* Generic settings: https://github.com/atgu/CHARR/blob/main/R/constants.R
* Final figures: https://github.com/atgu/CHARR/blob/main/R/Final_CHARR_figures.Rmd
* Drafts:
  * https://github.com/atgu/CHARR/blob/main/R/final_supplement_figures.Rmd
  * https://github.com/atgu/CHARR/blob/main/R/supplement_figures.Rmd
  * https://github.com/atgu/CHARR/blob/main/R/contamination_results.Rmd 
