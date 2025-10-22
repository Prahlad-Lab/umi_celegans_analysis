#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name="TransferReadTags"
#SBATCH --output=stdout/task.TransferReadTagsGenomeNoClip.out
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
	
# Done already
# parallel --xapply  gatk SortSam \
      # INPUT=${INPUT_DIR}/{}.extractUMIs.out.bam  \
      # OUTPUT=${INPUT_DIR}/{}.querysort.extractUMIs.out.bam  \
      # SORT_ORDER="queryname" \
      # CREATE_MD5_FILE=true \
      # MAX_RECORDS_IN_RAM=300000  ::: "${sample[@]}" ;
  
parallel --xapply gatk TransferReadTags \
    -I ${OUTPUT_DIR}/{}.querysort.Aligned.out.bam \
    --unmapped-sam ${INPUT_DIR}/{}.querysort.extractUMIs.out.bam  \
    -O ${OUTPUT_DIR}/{}.ReadTags.querysort.Aligned.out.bam \
    --read-tags RX ::: "${sample[@]}" ;
