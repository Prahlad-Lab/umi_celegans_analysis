#!/bin/bash -l
#
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name="bcftools.mpileup.singletons"
#SBATCH --output=stdout/task.bcftoolsMpileupSingletons.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
#SBATCH --mem=500G

module load gcc python

# --- Define Directories ---
# Input directory is the output from the GATK BQSR pipeline
gatk_output_dir="/vscratch/grp-vprahlad/umi_analysis/data/output_broad/Family_size1_haplotype_caller"
# Output directory for this bcftools step
bcftools_output_dir="/vscratch/grp-vprahlad/umi_analysis/data/output_broad/Family_size1_bcftools_mpileup"
# Directory for reference genome
ref_dir="/projects/rpci/vprahlad/umi_celegans/data/input"


# --- Prepare output directory ---
if [ ! -d $bcftools_output_dir ]; then      # if directory doesn't exist
    mkdir -p ${bcftools_output_dir}         # create it
fi

# --- Sample List ---
sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3")

# --- Activate Conda Environment ---
conda activate bcftools

echo "--- Starting bcftools mpileup for all samples ---"
# --- Process each sample ---
# Need to adjust --max-depth(-d) to keep all reads at each locus, the default is 250
for i in ${sample[@]}
do
    echo "Processing sample: ${i}"
    bcftools mpileup \
        -d 1000000 \
        -f ${ref_dir}/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
        ${gatk_output_dir}/"${i}".singletons.BSQR.bam > ${bcftools_output_dir}/"${i}".singletons.bcftools.pileup 2>&1 &
done

wait # Wait for all background processes to finish
echo "All bcftools mpileup tasks completed."
