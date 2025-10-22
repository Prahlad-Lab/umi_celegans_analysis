#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name="MergeBamAlignment"
#SBATCH --output=stdout/task.MergeBamAlignmentGenomeNoClip.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
	
conda activate gatk4	

INPUT_DIR=$"/vscratch/grp-vprahlad/umi_analysis/data/output_broad/ExtractUmisFromBam"
OUTPUT_DIR=$"/vscratch/grp-vprahlad/umi_analysis/data/output_broad/UMIAwareDuplicateMarkingGenomeNoClip"

# Prepare output directory
if [ ! -d $OUTPUT_DIR ]; then      # if doesn't exit
    mkdir ${OUTPUT_DIR}            # create it
fi ;

sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3") 	

parallel --xapply gatk     MergeBamAlignment \
            --REFERENCE_SEQUENCE /projects/rpci/vprahlad/umi_celegans/analysis_broad/../data/input/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa \
            --UNMAPPED_BAM ${INPUT_DIR}/{}.querysort.extractUMIs.out.bam \
            --ALIGNED_BAM ${OUTPUT_DIR}/{}.querysort.Aligned.out.bam  \
            --OUTPUT ${OUTPUT_DIR}/{}.MergeBamAlignment.bam \
            --INCLUDE_SECONDARY_ALIGNMENTS false \
            --VALIDATION_STRINGENCY SILENT ::: "${sample[@]}" ;