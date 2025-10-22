#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name="FastqToUbam"
#SBATCH --output=stdout/task.FastqToUbam.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc


conda activate gatk4

# Define input and output directories
INPUT_DIR=$"../data/seqs" # Replace with your actual input directory
OUTPUT_DIR=$"/vscratch/grp-vprahlad/umi_analysis/data/output_broad/FastqToUbam" # Replace with your actual output directory

# Prepare output directory
if [ ! -d $OUTPUT_DIR ]; then      # if doesn't exit
    mkdir ${OUTPUT_DIR}            # create it
fi


# Create a bash array with sample names
sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3") # Add your sample names here

parallel --xapply FastqToSam \
      SORT_ORDER=unsorted \
      F1=${INPUT_DIR}/{}_R1.fastq.gz \
      F2=${INPUT_DIR}/{}_R2.fastq.gz \
      SM={} \
      LB="L1" \
      PL="ILLUMINA" \
      PU="barcode1" \
      RG="RG1" \
      CN="BI" \
      O=${OUTPUT_DIR}/{}.unmapped.bam \
      ALLOW_EMPTY_FASTQ=True ::: "${sample[@]}"