#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name="GATK_BQSR_Pipeline_by_Step"
#SBATCH --output=stdout/task.GATK_BQSR_Pipeline.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
#SBATCH --mem=100GB

# --- Define Directories ---
vcf_dir="/vscratch/grp-vprahlad/umi_analysis/data/input"
input_dir="/projects/rpci/vprahlad/umi_celegans/data/input"
bam_dir="/vscratch/grp-vprahlad/umi_analysis/data/output_broad/FilterBambyFamilySizeGenome"
output_dir="/vscratch/grp-vprahlad/umi_analysis/data/output_broad/Family_size1_haplotype_caller"

# --- Reference Files ---
REFERENCE_FASTA="${input_dir}/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
REFERENCE_DICT="${input_dir}/Caenorhabditis_elegans.WBcel235.dna.toplevel.dict"
KNOWN_SITES_VCF="${vcf_dir}/CB4856.hard-filter.vcf.gz"

# --- Prepare output directory ---
if [ ! -d $output_dir ]; then      # if directory doesn't exist
    mkdir -p ${output_dir}         # create it
fi

# --- Sample List ---
sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3")

# --- Activate Conda Environment ---
conda activate gatk4

# ---
# STEP 1: Run SplitNCigarReads for ALL samples
# ---
echo "--- Starting Step 1: SplitNCigarReads for all samples ---"
for i in ${sample[@]}
do
    echo "Step 1 - Processing sample: ${i}"
    # Prepares RNA-seq data for HaplotypeCaller by splitting reads with N-CIGAR operations.
    gatk SplitNCigarReads \
        -R ${REFERENCE_FASTA} \
        -I ${bam_dir}/"${i}".sort_coood_singletons.bam \
        -O ${output_dir}/"${i}".singletons.ready.haplotype.caller.bam
done
echo "--- Completed Step 1 for all samples. ---"


# ---
# STEP 2: Run BaseRecalibrator for ALL samples
# ---
echo "--- Starting Step 2: BaseRecalibrator for all samples ---"
for i in ${sample[@]}
do
    echo "Step 2 - Processing sample: ${i}"
    # Generates a recalibration table to detect and correct systematic errors in base quality scores.
    gatk BaseRecalibrator \
        -R ${REFERENCE_FASTA} \
        -I ${output_dir}/"${i}".singletons.ready.haplotype.caller.bam \
        --use-original-qualities \
        -O ${output_dir}/"${i}".singletons.recal.table \
        --known-sites ${KNOWN_SITES_VCF}
done
echo "--- Completed Step 2 for all samples. ---"


# ---
# STEP 3: Run ApplyBQSR for ALL samples
# ---
echo "--- Starting Step 3: ApplyBQSR for all samples ---"
for i in ${sample[@]}
do
    echo "Step 3 - Processing sample: ${i}"
    # Applies the recalibration table to the BAM file, adjusting the base quality scores.
    gatk ApplyBQSR \
        -R ${REFERENCE_FASTA} \
        -I ${output_dir}/"${i}".singletons.ready.haplotype.caller.bam \
        -O ${output_dir}/"${i}".singletons.BSQR.bam \
        --bqsr-recal-file ${output_dir}/"${i}".singletons.recal.table \
        -OBI \
        --sequence-dictionary ${REFERENCE_DICT}
done
echo "--- Completed Step 3 for all samples. ---"


echo "All pipeline steps have been successfully processed."

