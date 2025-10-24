#!/bin/bash

# ---
# Conda Environment Setup for R Variant Analysis
#
# This script creates a new conda environment named 'r_variant_analysis'
# with all the required R packages (from CRAN and Bioconductor)
# to run both R scripts.
# ---

# Define the name of the environment
ENV_NAME="r_variant_analysis"

echo "--- Creating Conda environment: $ENV_NAME ---"

# Create the environment with all packages
# We specify the channels in order of priority:
# 1. conda-forge: For R base and most CRAN packages
# 2. bioconda: For Bioconductor packages
# 3. defaults: For any remaining dependencies
conda create -n $ENV_NAME \
  -c conda-forge \
  -c bioconda \
  -c defaults \
  'r-base>=4.2' \
  'r-tidyverse' \
  'r-scales' \
  'r-ggpubr' \
  'r-cowplot' \
  'r-extrafont' \
  'r-plotly' \
  'r-ggprism' \
  'r-vcfr' \
  'r-rstatix' \
  'r-knitr' \
  'r-kableextra' \
  'r-officer' \
  'r-flextable' \
  'bioconductor-rtracklayer' \
  'bioconductor-plyranges' \
  -y

# Check if creation was successful
if [ $? -eq 0 ]; then
  echo ""
  echo "--- Environment '$ENV_NAME' created successfully. ---"
  echo ""
  echo "To activate this environment, run:"
  echo "  conda activate $ENV_NAME"
  echo ""
  echo "To deactivate it, run:"
  echo "  conda deactivate"
  echo ""
else
  echo ""
  echo "--- ERROR: Failed to create Conda environment '$ENV_NAME'. ---"
  echo "Please check your Conda installation and for any package conflicts."
  echo ""
fi
