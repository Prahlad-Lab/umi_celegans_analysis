#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name="calculate_allele_proportions_singletons"
#SBATCH --output=stdout/task.calculate_allele_proportions_singletons.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc


module load gcc python

# --- Define Directories ---
# Input directory is the output from the bcftools mpileup script
pileup_dir="/vscratch/grp-vprahlad/umi_analysis/data/output_broad/Family_size1_bcftools_mpileup"
# Output directory for this allele proportion calculation step
proportions_dir="/vscratch/grp-vprahlad/umi_analysis/data/output_broad/Family_size1_allele_proportions"

# --- NOTE ---
# The original script assumed the python script 'calculate_allele_proportions.py'
# was located in the same directory as the input pileup files.
# This script maintains that assumption. If your script is located elsewhere,
# you may need to adjust the path in the python command below.
SCRIPT_PATH="/projects/rpci/vprahlad/umi_celegans/analysis_broad/calculate_allele_proportions.py"


# --- Prepare output directory ---
if [ ! -d $proportions_dir ]; then      # if directory doesn't exist
    mkdir -p ${proportions_dir}         # create it
fi

# --- Sample List ---
sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3")

echo "--- Starting Allele Proportion Calculation ---"
for i in ${sample[@]}
do
    echo "Processing sample: ${i}"
    python ${SCRIPT_PATH} \
        ${pileup_dir}/"${i}".singletons.bcftools.pileup \
        ${proportions_dir}/"${i}".singletons.per_position_results.tsv \
        ${proportions_dir}/"${i}".singletons.summary_totals.txt
done

echo "--- Allele proportion calculation complete for all samples. ---"
