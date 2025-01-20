# SensSpec

A comprehensive R script for calculating and visualizing sensitivity and specificity metrics for variant calling pipelines.

## Overview

SensSpec is a tool designed to assess the performance of variant calling pipelines by calculating sensitivity, specificity, and other important metrics. It compares pipeline output against a truth set and generates detailed performance metrics and visualizations.

## Features

- Calculates key performance metrics:
  - Sensitivity (True Positive Rate)
  - Specificity
  - False Positive Rate
  - Positive Predictive Value (PPV)
  - F-measure
- Analyzes both SNPs and INDELs separately
- Evaluates performance across different:
  - Variant Allele Frequencies (VAF)
  - Read Depths
- Generates ROC curves for visual performance assessment
- Supports command-line operation for easy integration into pipelines

## Requirements

Required R packages:
- ggplot2
- vcfR
- dbplyr
- sqldf
- reshape2
- optparse
- gridExtra
- grid
- lattice
- tidyverse

## Installation

```R
# Install required packages
required_packages <- c("ggplot2", "vcfR", "dbplyr", "sqldf", "reshape2", 
                      "optparse", "gridExtra", "grid", "lattice", "tidyverse")
install.packages(required_packages)
```

## Usage

```bash
Rscript SensSpec.R -f1 <truth_set_file> -f2 <test_file> -b <total_base_pairs>
```

### Parameters

- `-f1, --file1`: Truth set file (reference dataset)
- `-f2, --file2`: Test file (pipeline output to evaluate)
- `-b, --bp`: Total number of base pairs in the analyzed region

### Input File Requirements

#### Truth Set File (file1)
Must contain columns:
- Chromosome (in chr1 format)
- Start_position
- Reference_Allele
- Tumor_Seq_Allele2
- Variant_Type

#### Test File (file2)
Must contain columns:
- Samplename
- Chrom
- Pos
- Ref
- Alt
- Variant.Type.SnpEff
- FDP (Flow Depth)
- FRD (Flow Reference Depth)
- FAD (Flow Alternate Depth)
- FAF (Flow Allele Frequency)

## Output

The script generates:
1. ROC curves comparing SNP and INDEL detection performance
2. Performance metrics at various VAF thresholds
3. Performance metrics at various depth thresholds

## Example

```bash
Rscript SensSpec.R -f1 truth_set.maf -f2 pipeline_output.tsv -b 493868
```

## Authors

- Chase Rushton
- A. Bigdeli

## License

This project is open source. Please contact the repository owner for specific licensing information.

---
Last Updated: January 20, 2025
