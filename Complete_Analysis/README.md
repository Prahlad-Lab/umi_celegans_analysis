# UMI Analysis Pipeline 2 - Complete Analysis

## Overview

The `UMI_analysis_pipeline_2.sh` script is a comprehensive bioinformatics pipeline for analyzing RNA-seq data with Unique Molecular Identifiers (UMIs) in *C. elegans*. This pipeline performs quality control, alignment, UMI-aware duplicate marking, variant calling, and allele proportion analysis.

## Pipeline Workflow

The pipeline consists of the following major steps:

1. **FastQC** - Quality control of raw sequencing data
2. **FastqToUbam** - Convert FASTQ files to unmapped BAM format
3. **ExtractUMIs** - Extract UMI sequences from reads
4. **Sort by QueryName** - Sort UMI-extracted BAM files
5. **STAR Alignment** - Align reads to the reference genome
6. **MergeBamAlignment** - Merge aligned and unmapped BAM files
7. **GroupByUMIs** - Group reads by UMI to identify read families
8. **Count Family Sizes** - Generate family size statistics
9. **Filter by Family Size** - Create BAM files for singletons (family size = 1)
10. **GATK BQSR Pipeline** - Split N-cigars, base quality score recalibration
11. **Variant Calling** (Path A) - HaplotypeCaller and SnpEff annotation
12. **Allele Proportion Analysis** (Path B) - bcftools mpileup and allele proportion calculation

## Prerequisites

### Required Software and Tools

The pipeline uses the following conda environments and tools:

- **fastqc** - For quality control
- **gatk4** - GATK toolkit for genomic analysis
- **fgbio** - Tools for working with genomic data
- **STAR** - RNA-seq aligner
- **umi_tools** - UMI handling tools
- **samtools.v1.22** - SAMtools for BAM manipulation
- **bcftools** - For variant calling
- **snpeff** - For variant annotation

### Required System Modules

- `gcc` - GNU Compiler Collection
- `python` - Python 3

### Input Data Requirements

1. **FASTQ Files**: Paired-end sequencing data
   - Located in: `../data/seqs/`
   - Format: `{sample}_R1.fastq.gz` and `{sample}_R2.fastq.gz`

2. **Reference Genome Files**:
   - Reference genome FASTA: `Caenorhabditis_elegans.WBcel235.dna.toplevel.fa`
   - Reference dictionary: `Caenorhabditis_elegans.WBcel235.dna.toplevel.dict`
   - STAR index directory

3. **Known Sites VCF**:
   - CB4856 variants: `CB4856.hard-filter.vcf.gz`

4. **Python Script**:
   - Allele proportion calculator: `calculate_allele_proportions_depth10.py`

## Configuration

### Customizing the Pipeline

Before running the pipeline, you need to configure the following variables in the script:

#### 1. Sample Names

Edit the sample array to include your sample identifiers:

```bash
sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3")
```

#### 2. Input Directories

Update the paths to match your data locations:

```bash
DIR_SEQS="../data/seqs"                    # Location of FASTQ files
DIR_REF_INPUT="<path_to_reference_files>"  # Reference genome and dictionary
DIR_VCF_INPUT="<path_to_vcf_files>"        # Known sites VCF file
```

#### 3. Reference Files

Ensure the following reference files are available:

```bash
REF_GENOME="${DIR_REF_INPUT}/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
REF_DICT="${DIR_REF_INPUT}/Caenorhabditis_elegans.WBcel235.dna.toplevel.dict"
KNOWN_SITES_VCF="${DIR_VCF_INPUT}/CB4856.hard-filter.vcf.gz"
STAR_INDEX="../data/output/star_index"
```

#### 4. Python Script Path

Update the path to the allele proportion calculation script:

```bash
PYTHON_SCRIPT_PATH="<path_to>/calculate_allele_proportions_depth10.py"
```

#### 5. Output Directory

Configure the base output directory:

```bash
BASE_OUTPUT_DIR="/vscratch/grp-vprahlad/umi_analysis/data/output_broad"
```

#### 6. SLURM Configuration

Adjust SLURM parameters based on your cluster configuration:

```bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500G
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
```

Update the email notification:

```bash
#SBATCH --mail-user=your.email@institution.org
```

## Directory Structure

The pipeline will create the following output directory structure:

```
BASE_OUTPUT_DIR/
├── fastqc/                              # FastQC quality reports
├── FastqToUbam/                         # Unmapped BAM files
├── ExtractUmisFromBam/                  # UMI-extracted BAM files
├── STARNoClip/                          # STAR alignment outputs
├── UMIAwareDuplicateMarkingGenomeNoClip/# UMI-grouped and merged BAMs
├── FilterBambyFamilySizeGenome/         # Filtered BAM files (singletons)
├── Family_size1_haplotype_caller/       # Variant calling outputs
├── Family_size1_annotation/             # SnpEff annotation outputs
├── Family_size1_bcftools_mpileup/       # bcftools mpileup outputs
└── Family_size1_allele_proportions_depth10/ # Allele proportion results
```

## Running the Pipeline

### On a SLURM Cluster

Submit the job using `sbatch`:

```bash
sbatch UMI_analysis_pipeline_2.sh
```

### Monitoring Progress

Check the job status:

```bash
squeue -u $USER
```

View the output log:

```bash
tail -f stdout/Combined_UMI_Pipeline.out
```

### Running Specific Steps

If you need to run specific steps individually, you can comment out unwanted sections in the script using `#` or extract specific step commands.

## Output Files

### Key Output Files by Step:

1. **FastQC Reports** (`DIR_FASTQC`):
   - `*.html` - HTML quality reports
   - `*.zip` - Detailed quality metrics

2. **UMI Family Analysis** (`DIR_UMI_MARKING`):
   - `{sample}.grouped.tsv` - UMI group information
   - `{sample}.family_size_summary.tsv` - Family size distribution

3. **Filtered BAMs** (`DIR_FILTER_BAM`):
   - `{sample}.sort_coood_singletons.bam` - Singleton reads only
   - `{sample}.singletons_id.txt` - Singleton UMI identifiers

4. **Variant Calls** (`DIR_HAPLOTYPE_CALLER`):
   - `{sample}.singletons.first.vcf.gz` - Raw variant calls
   - `{sample}.singletons.BSQR.bam` - Recalibrated BAM files

5. **Variant Annotations** (`DIR_ANNOTATION`):
   - `{sample}.singletons.ann.vcf` - Annotated variants
   - `{sample}_singleton_summary.html` - SnpEff summary report
   - `{sample}.singletons.ann.bed` - BED format annotations

6. **Allele Proportions** (`DIR_ALLELE_PROPORTIONS`):
   - `{sample}.singletons.per_position_results.tsv` - Per-position allele data
   - `{sample}.singletons.summary_totals.txt` - Summary statistics

## UMI Read Structure

The pipeline uses the following UMI read structure:

- **Read 1**: `+T` (full read)
- **Read 2**: `8M6S+T` (8 bp UMI, 6 bp skip, then template)

This extracts an 8-base UMI from Read 2 with a 6-base skip region.

## Resource Requirements

### Computational Resources:

- **CPUs**: 32 cores (configurable)
- **Memory**: 500 GB RAM (may need adjustment based on data size)
- **Storage**: Significant disk space required for intermediate and output files
- **Runtime**: Several hours to days depending on:
  - Number of samples
  - Sequencing depth
  - Cluster load

### Disk Space Estimates:

- FASTQ files: Variable (input)
- Intermediate BAM files: ~5-10x input size
- Final outputs: ~2-3x input size

## Troubleshooting

### Common Issues:

1. **Conda environment not found**:
   - Ensure all required conda environments are installed
   - Activate environments manually to test: `conda activate <env_name>`

2. **File not found errors**:
   - Verify all input paths are correct
   - Check that reference files exist at specified locations
   - Ensure FASTQ files follow the naming convention

3. **Out of memory errors**:
   - Increase `--mem` in SLURM header
   - Reduce `MAX_RECORDS_IN_RAM` in GATK SortSam steps
   - Process fewer samples in parallel

4. **STAR alignment fails**:
   - Verify STAR index was built correctly
   - Check that reference genome matches the index
   - Ensure sufficient disk space in temp directory

5. **Permission denied**:
   - Verify write permissions to output directories
   - Check that `mkdir -p` commands can create directories

### Log Files:

- **Main log**: `stdout/Combined_UMI_Pipeline.out`
- **UMI grouping logs**: `{DIR_UMI_MARKING}/{sample}.group.log`
- **STAR logs**: `{DIR_STAR_ALIGN}/{sample}.Log.final.out`

## Notes

- The pipeline processes multiple samples in parallel using GNU `parallel`
- Some steps run sequentially (for loops) while others use parallel execution
- The pipeline forks into two paths after BQSR:
  - **Path A**: GATK variant calling and SnpEff annotation
  - **Path B**: bcftools mpileup and allele proportion analysis
- Both paths run independently and can be monitored separately

## Additional Resources

For more information on individual tools:

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [GATK](https://gatk.broadinstitute.org/)
- [fgbio](https://fulcrumgenomics.github.io/fgbio/)
- [STAR](https://github.com/alexdobin/STAR)
- [UMI-tools](https://umi-tools.readthedocs.io/)
- [bcftools](http://samtools.github.io/bcftools/)
- [SnpEff](http://pcingola.github.io/SnpEff/)

## Citation

If you use this pipeline, please cite the appropriate tools and methods used in your analysis.
