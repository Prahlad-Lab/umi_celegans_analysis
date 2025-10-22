# Quick Start Guide - Conda Environment Setup

This guide provides quick commands to set up the required conda environments for the UMI C. elegans analysis pipeline.

## TL;DR - Quick Setup

```bash
# Make the script executable (if needed)
chmod +x setup_conda_environments.sh

# Create all environments in one combined environment
./setup_conda_environments.sh all

# OR create separate environments for each tool
./setup_conda_environments.sh separate
```

## Required Tools and Versions

All tools are installed from the **bioconda** channel:

| Tool       | Version   | Purpose                                    |
|------------|-----------|-------------------------------------------|
| fastqc     | 0.12.1    | Quality control for sequencing data       |
| bcftools   | 1.22      | VCF/BCF manipulation                      |
| umi_tools  | 1.1.6     | UMI handling                              |
| gatk4      | 4.6.1.0   | Variant discovery                         |
| fgbio      | 2.5.0     | Genomic/sequencing data tools             |
| snpeff     | 5.2       | Variant annotation                        |
| samtools   | 1.22.1    | SAM/BAM manipulation                      |
| STAR       | 2.7.11b   | RNA-seq alignment                         |

## Setup Options

### Option 1: Combined Environment (Easiest)

One environment with all tools:

```bash
./setup_conda_environments.sh all
conda activate umi_celegans_analysis
```

**Pros:** Simple, all tools in one place  
**Cons:** Potential for dependency conflicts

### Option 2: Separate Environments (Recommended)

Individual environments for each tool (matches pipeline structure):

```bash
./setup_conda_environments.sh separate
```

Activate environments as needed:
```bash
conda activate fastqc
conda activate gatk4
conda activate umi_tools
# etc.
```

**Pros:** Avoids dependency conflicts, matches pipeline  
**Cons:** More environments to manage

## Manual Installation

### Using the YAML file:
```bash
conda env create -f environment.yml
```

### Create individual environments:
```bash
conda create -n fastqc -c bioconda fastqc=0.12.1 -y
conda create -n bcftools -c bioconda bcftools=1.22 -y
conda create -n umi_tools -c bioconda umi_tools=1.1.6 -y
conda create -n gatk4 -c bioconda gatk4=4.6.1.0 -y
conda create -n fgbio -c bioconda fgbio=2.5.0 -y
conda create -n snpeff -c bioconda snpeff=5.2 -y
conda create -n samtools.v1.22 -c bioconda samtools=1.22.1 -y
conda create -n STAR -c bioconda star=2.7.11b -y
```

## Verifying Installation

Check that environments are created:
```bash
conda env list
```

Test a tool (after activating its environment):
```bash
conda activate fastqc
fastqc --version
```

## Troubleshooting

### Conda not found
```bash
# Add conda to PATH (adjust path as needed)
export PATH="$HOME/miniconda3/bin:$PATH"
source ~/.bashrc
```

### Environment already exists
The setup script will prompt you to remove and recreate existing environments.

Or manually remove:
```bash
conda env remove -n umi_celegans_analysis
```

### Installation fails
Try updating conda first:
```bash
conda update -n base conda
```

### Slow installation
Add mamba for faster package resolution:
```bash
conda install -n base -c conda-forge mamba
mamba env create -f environment.yml
```

## Next Steps

After installation, see the main [README.md](README.md) for:
- Pipeline usage instructions
- Analysis workflow details
- Directory structure
- Script descriptions

## Getting Help

Run the setup script with help flag:
```bash
./setup_conda_environments.sh help
```
