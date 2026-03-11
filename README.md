# GCAD VCF QC Pipeline

This pipeline performs **variant and sample-level quality control (QC)** on VCF files and produces summary statistics used for downstream filtering and dataset evaluation. The workflow processes genotype data, applies configurable QC thresholds, and outputs both **sample-level** and **variant-level** metrics.

The pipeline supports:

- Whole Genome Sequencing (WGS)
- Family-based datasets
- Multiallelic variant handling
- Mendelian inconsistency calculations

---

# Inputs

The pipeline requires the following inputs:

| Argument | Description |
|---|---|
| `--vcf` | Input VCF file containing genotype calls |
| `--fam` | FAM file containing sample metadata (must include headers) |
| `--out_dir` | Output directory for pipeline results (default: `./`) |

Optional inputs:

| Argument | Description |
|---|---|
| `--region` | Restrict processing to a genomic region (example: `chr3:100-200`) |

---

# Parameters

## Configuration Defaults

The pipeline uses the following default QC thresholds:

| Parameter | Default Value | Description |
|-----------|--------------|-------------|
| `MINDP` | `10` | Minimum genotype depth (DP) |
| `MINGQ` | `20` | Minimum genotype quality (GQ) |
| `hwe_maf` | `0.01` | Minor allele frequency threshold used for HWE filtering |
| `minTranche` | `99.7` | Minimum VQSR tranche score |
| `miss_rate` | `0.2` | Maximum allowed genotype missingness |
| `max_dp` | `500` | Maximum mean depth threshold |
| `isFam` | `1` | Indicates family-based dataset (1 = yes, 0 = no) |
| `hwe_pval` | `5e-06` | Hardy-Weinberg equilibrium p-value threshold |
| `flank_size` | `7` | Flanking region size (bp) used for WES target intervals |

## Genotype Quality Filters
The pipeline allows configurable QC thresholds.
| Parameter | Description | Default |
|---|---|---|
| `--min_dp` | Minimum genotype depth (DP) | `10` |
| `--min_gq` | Minimum genotype quality (GQ) | `20` |
| `--max_dp` | Maximum mean depth threshold | `500` |

## Dataset Settings

| Parameter | Description | Default |
|---|---|---|
| `--is_fam` | Indicates family data (1 = yes, 0 = no) | `1` |
| `--is_multiallelic` | Flag indicating multiallelic VCF | `False` |
| `--is_chrx` | Flag indicating chromosome X processing | `False` |

## Variant Filtering

| Parameter | Description | Default |
|---|---|---|
| `--min_tranche` | Minimum VQSR tranche score | `99.7` |
| `--miss_rate` | Maximum allowed genotype missingness | `0.2` |
| `--hwe_pval` | Hardy-Weinberg p-value threshold | `5e-06` |
| `--hwe_maf` | MAF threshold for HWE filtering | `0.01` |

## Heterozygosity Limits

| Parameter | Description | Default |
|---|---|---|
| `--hetz_lim1` | Heterozygosity limit parameter | `99999` |
| `--hetz_lim2` | Heterozygosity limit parameter | `99999` |

---


# Outputs

The pipeline generates TSV summary files containing QC statistics.

---

## Sample-Level Output

### `summary.snv.tsv`

This file contains **one row per sample** summarizing genotype counts, missingness, depth metrics, singleton statistics, and other QC indicators.

### Columns

| Column | Description |
|---|---|
| `SampleID` | Sample identifier |
| `SEX` | Sex value from the FAM metadata |
| `total_nRR_Bi` | Number of passing homozygous reference genotypes (biallelic sites) |
| `total_nRA_Bi` | Number of passing heterozygous genotypes (biallelic sites) |
| `total_nAA_Bi` | Number of passing homozygous alternate genotypes (biallelic sites) |
| `Missing_Bi` | Number of missing genotypes at biallelic sites |
| `Set_Missing_Bi` | Genotypes set to missing due to QC filters |
| `IndDepthSum_Bi` | Total read depth across biallelic sites |
| `IndMeanDepth_Bi` | Mean read depth across biallelic sites |
| `total_nRR_Multi` | Homozygous reference genotypes at multiallelic sites |
| `total_nRA_Multi` | Heterozygous genotypes at multiallelic sites |
| `total_nAA_Multi` | Homozygous alternate genotypes at multiallelic sites |
| `Missing_Multi` | Missing genotypes at multiallelic sites |
| `Set_Missing_Multi` | Genotypes set to missing at multiallelic sites |
| `IndDepthSum_Multi` | Total depth across multiallelic sites |
| `IndMeanDepth_Multi` | Mean depth across multiallelic sites |
| `Singleton_Bi` | Singleton variants at biallelic sites |
| `Private_Doubleton_Bi` | Private doubletons at biallelic sites |
| `Doubleton_Bi` | Doubletons at biallelic sites |
| `Singleton_Multi` | Singleton variants at multiallelic sites |
| `Private_Doubleton_Multi` | Private doubletons at multiallelic sites |
| `Doubleton_Multi` | Doubletons at multiallelic sites |
| `HetHom_Bi` | Heterozygous-to-homozygous alternate ratio (biallelic) |
| `HetHom_Multi` | Heterozygous-to-homozygous alternate ratio (multiallelic) |
| `Ti_Bi` | Transition count (biallelic) |
| `Tv_Bi` | Transversion count (biallelic) |
| `TiTvRatio_Bi` | Transition/transversion ratio (biallelic) |
| `Ti_Multi` | Transition count (multiallelic) |
| `Tv_Multi` | Transversion count (multiallelic) |
| `TiTvRatio_Multi` | Transition/transversion ratio (multiallelic) |
| `1P_MI` | Mendelian inconsistencies in one-parent families |
| `2P_MI` | Mendelian inconsistencies in two-parent families |
| `MI_pairs` | Number of informative Mendelian pairs |
| `Non_Missing_Indels_Bi` | Number of non-missing indel genotypes |

---

## Variant-Level Output

### `summary.by-sample.tsv`

This file contains **one row per variant site** summarizing genotype pass/fail counts, allele frequencies, and QC flags.

### Columns

| Column | Description |
|---|---|
| `CHR` | Chromosome |
| `POS` | Genomic position |
| `PassRR` | Passing homozygous reference calls |
| `PassRA` | Passing heterozygous calls |
| `PassAA` | Passing homozygous alternate calls |
| `FailRR` | Failing homozygous reference calls |
| `FailRA` | Failing heterozygous calls |
| `FailAA` | Failing homozygous alternate calls |
| `Missing` | Missing genotype calls |
| `GT_Failed` | Genotypes failing genotype-level QC |
| `rsID` | Variant rsID |
| `RefAllele` | Reference allele |
| `AltAllele` | Alternate allele(s) |
| `VTYPE` | Variant type (SNV or indel) |
| `MultiAllele` | Indicator for multiallelic sites |
| `QUAL` | Variant quality score |
| `GLNexusPass` | Indicator of GLNexus filter pass |
| `Mono` | Indicator of monomorphic site |
| `CallRate` | Fraction of non-missing genotypes |
| `CallBad` | Call rate below missingness threshold |
| `AltAF` | Alternate allele frequency |
| `MeanDepth` | Mean depth across samples |
| `HiDepth` | Mean depth exceeds threshold |
| `ABHet` | Allele balance for heterozygous calls |
| `RecommendedDrop` | Site recommended for removal |
| `Mend_Incon` | Mendelian inconsistencies |
| `Mend_pairs` | Informative Mendelian comparisons |
| `propMI` | Proportion of Mendelian inconsistencies |

---

# Example Usage

Run the pipeline:
```
python3 /path/to/gcad_vcf_qc.py --fam /path/to/famfile --vcf /path/to/vcffile
```

Run the pipeline on ChrX:
```
python3 /path/to/gcad_vcf_qc.py --fam /path/to/famfile --vcf /path/to/vcffile --is_chrX
```

Run pipeline on specific region:
```
python3 /path/to/gcad_vcf_qc.py --fam /path/to/famfile --vcf /path/to/vcffile --region chr10:1000000-2000000
```
