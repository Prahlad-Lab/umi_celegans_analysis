#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name="SortNIndexFilterBambyFamilySize"
#SBATCH --output=stdout/task.SortNIndexFilterBambyFamilySize.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc


# No modules needed for this script

INPUT_DIR="/vscratch/grp-vprahlad/umi_analysis/data/output_broad/UMIAwareDuplicateMarkingGenomeNoClip"
OUTPUT_DIR="/vscratch/grp-vprahlad/umi_analysis/data/output_broad/FilterBambyFamilySizeGenome"

sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3")

# Prepare output directory
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi


conda activate samtools.v1.22

for i in "${sample[@]}"
do
    echo "Processing sample: $i"
	samtools sort -@ 8   -o "${OUTPUT_DIR}/${i}.sort_coood_families_2-10.bam" "${OUTPUT_DIR}/${i}.families_2-10.bam";
	samtools index -@ 8 "${OUTPUT_DIR}/${i}.sort_coood_families_2-10.bam"	
done;
echo "All families_2-10 BAM created."

for i in "${sample[@]}"
do
    echo "Processing sample: $i"
	samtools sort -@ 8   -o "${OUTPUT_DIR}/${i}.sort_coood_singletons.bam" "${OUTPUT_DIR}/${i}.singletons.bam";
	samtools index -@ 8 "${OUTPUT_DIR}/${i}.sort_coood_singletons.bam" 
done;


echo "All singletons BAM created."
