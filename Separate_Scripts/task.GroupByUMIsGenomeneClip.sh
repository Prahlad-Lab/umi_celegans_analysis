#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name="GroupByUMIs"
#SBATCH --output=stdout/task.GroupByUMIsGenomeNoClip.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
#SBATCH --mem=200G

	
module load gcc python	


INPUT_DIR=$"/vscratch/grp-vprahlad/umi_analysis/data/output_broad/ExtractUmisFromBam"
OUTPUT_DIR=$"/vscratch/grp-vprahlad/umi_analysis/data/output_broad/UMIAwareDuplicateMarkingGenomeNoClip"

# Prepare output directory
if [ ! -d $OUTPUT_DIR ]; then      # if doesn't exit
    mkdir ${OUTPUT_DIR}            # create it
fi ;

sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3") 	
 

conda activate umi_tools
export PYTHONHASHSEED=0
 
parallel --xapply umi_tools group -I ${OUTPUT_DIR}/{}.MergeBamAlignment.coordinate_sorted.bam  --paired --no-sort-output --output-bam --stdout ${OUTPUT_DIR}/{}.grouped_by_UMI.bam --umi-tag-delimiter "-" \
    --log ${OUTPUT_DIR}/{}.group.log  --group-out ${OUTPUT_DIR}/{}.grouped.tsv  --extract-umi-method tag --umi-tag RX --unmapped-reads use ::: "${sample[@]}" ;
