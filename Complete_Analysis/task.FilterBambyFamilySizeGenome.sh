#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name="FilterBambyFamilySizeGenome"
#SBATCH --output=stdout/task.FilterBambyFamilySizeGenome.out
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

for i in "${sample[@]}"
do
    echo "Processing sample: $i"
	awk 'BEGIN{FS="\t"} NR > 1 && $8 == 1 {print $9}'  "${INPUT_DIR}/${i}.grouped.tsv" | sort -u  > "${OUTPUT_DIR}/${i}.singletons_id.txt"
done;
echo "All singletons files created."

conda activate samtools.v1.22

for i in "${sample[@]}"
do
    echo "Processing sample: $i"
	
	samtools view -@ 8 -h --tag-file UG:"${OUTPUT_DIR}/${i}.singletons_id.txt" -o "${OUTPUT_DIR}/${i}.singletons.bam" "${INPUT_DIR}/${i}.grouped_by_UMI.bam"
done;

echo "All singletons BAM created."



for i in "${sample[@]}"
do
    echo "Processing sample: $i"
	samtools sort -@ 8   -o "${OUTPUT_DIR}/${i}.sort_coood_singletons.bam" "${OUTPUT_DIR}/${i}.singletons.bam";
	samtools index -@ 8 "${OUTPUT_DIR}/${i}.sort_coood_singletons.bam" 
done;

