#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name="extractUMIs"
#SBATCH --output=stdout/task.ExtractUmisFromBam.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc


conda activate fgbio


# Define input and output directories
INPUT_DIR=$"/vscratch/grp-vprahlad/umi_analysis/data/output_broad/FastqToUbam" # Replace with your actual input directory
OUTPUT_DIR=$"/vscratch/grp-vprahlad/umi_analysis/data/output_broad/ExtractUmisFromBam" # Replace with your actual output directory

# Prepare output directory
if [ ! -d $OUTPUT_DIR ]; then      # if doesn't exit
    mkdir ${OUTPUT_DIR}            # create it
fi

# Create a bash array with sample names
sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3") # Add your sample names here

parallel --xapply fgbio ExtractUmisFromBam \
      --input ${INPUT_DIR}/{}.unmapped.bam \
      --read-structure +T \
      --read-structure 8M6S+T  \
      --molecular-index-tags RX \
      --output  ${OUTPUT_DIR}/{}.extractUMIs.out.bam ::: "${sample[@]}"