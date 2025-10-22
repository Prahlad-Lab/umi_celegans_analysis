# UMI C. elegans Analysis Pipeline

This repository contains scripts and tools for analyzing UMI (Unique Molecular Identifier) data from *C. elegans* samples.

## Prerequisites

Before running the analysis pipeline, you need to install the required bioinformatics tools using conda.

### Required Tools

The pipeline requires the following tools from the bioconda channel:

- **fastqc** = 0.12.1 - Quality control for sequencing data
- **bcftools** = 1.22 - Utilities for variant calling and manipulating VCFs/BCFs
- **umi_tools** = 1.1.6 - Tools for handling Unique Molecular Identifiers
- **gatk4** = 4.6.1.0 - Genome Analysis Toolkit for variant discovery
- **fgbio** = 2.5.0 - Tools for working with genomic and sequencing data
- **snpeff** = 5.2 - Genetic variant annotation and effect prediction
- **samtools** = 1.22.1 - Tools for manipulating SAM/BAM files
- **STAR** = 2.7.11b - RNA-seq aligner

## Installation

### Step 1: Install Conda

If you don't have conda installed, download and install Miniconda or Anaconda:

- **Miniconda**: https://docs.conda.io/en/latest/miniconda.html
- **Anaconda**: https://www.anaconda.com/products/distribution

### Step 2: Set Up Conda Environments

We provide two options for setting up the required tools:

#### Option 1: Combined Environment (Recommended for Simple Setup)

Create a single conda environment with all tools:

```bash
./setup_conda_environments.sh all
```

This creates an environment called `umi_celegans_analysis` with all tools. Activate it with:

```bash
conda activate umi_celegans_analysis
```

#### Option 2: Separate Environments (Matches Pipeline Structure)

Create separate environments for each tool (as used in the pipeline scripts):

```bash
./setup_conda_environments.sh separate
```

This creates individual environments:
- `fastqc`
- `bcftools`
- `umi_tools`
- `gatk4`
- `fgbio`
- `snpeff`
- `samtools.v1.22`
- `STAR`

Activate each as needed:

```bash
conda activate fastqc
conda activate gatk4
# etc.
```

### Manual Installation

Alternatively, you can manually create the environment using the provided YAML file:

```bash
conda env create -f environment.yml
conda activate umi_celegans_analysis
```

Or install tools individually:

```bash
conda create -n fastqc -c bioconda -c conda-forge fastqc=0.12.1 -y
conda create -n bcftools -c bioconda -c conda-forge bcftools=1.22 -y
conda create -n umi_tools -c bioconda -c conda-forge umi_tools=1.1.6 -y
conda create -n gatk4 -c bioconda -c conda-forge gatk4=4.6.1.0 -y
conda create -n fgbio -c bioconda -c conda-forge fgbio=2.5.0 -y
conda create -n snpeff -c bioconda -c conda-forge snpeff=5.2 -y
conda create -n samtools.v1.22 -c bioconda -c conda-forge samtools=1.22.1 -y
conda create -n STAR -c bioconda -c conda-forge star=2.7.11b -y
```

## Usage

### Running the Complete Analysis Pipeline

The complete analysis pipeline is available in:
- `Complete_Analysis/UMI_analysis_pipeline_2.sh`

### Running Individual Steps

Individual pipeline steps are available in the `Separate_Scripts/` directory:
- `task.FastqToUbam.sh` - Convert FastQ to unmapped BAM
- `task.ExtractUMIs.sh` - Extract UMIs from BAM files
- `task.STARNoClip.sh` - STAR alignment
- `task.GroupByUMIsGenomeneClip.sh` - Group reads by UMI
- And more...

### Analysis Scripts

R scripts for downstream analysis are available in the `R_Scripts/` directory.

## Pipeline Overview

The pipeline performs the following steps:

1. **FastQC** - Quality control of raw sequencing data
2. **FastqToUbam** - Convert FastQ files to unmapped BAM format
3. **ExtractUMIs** - Extract UMI sequences from reads
4. **STAR Alignment** - Align reads to the reference genome
5. **UMI Grouping** - Group reads by their UMI sequences
6. **BAM Filtering** - Filter BAM files by family size (singletons, families)
7. **GATK BQSR Pipeline** - Base quality score recalibration
8. **Variant Calling** - Call variants using GATK HaplotypeCaller
9. **Variant Annotation** - Annotate variants using SnpEff
10. **Allele Proportion Calculation** - Calculate allele frequencies

## Directory Structure

```
.
├── Complete_Analysis/           # Complete pipeline scripts
├── Separate_Scripts/            # Individual task scripts
├── R_Scripts/                   # R analysis scripts
├── input_broad/                 # Reference files
├── Family_size1_allele_proportions_depth10/  # Output data
├── Family_size1_annotation/     # Variant annotations
├── environment.yml              # Conda environment specification
├── setup_conda_environments.sh  # Setup script
└── README.md                    # This file
```

## Support

For questions or issues, please open an issue on the GitHub repository.

## License

Please refer to the repository license file for usage terms.
