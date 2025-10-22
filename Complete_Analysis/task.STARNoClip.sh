#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name="STARNoClip"
#SBATCH --output=stdout/task.STARNoClip.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
#SBATCH --mem=100G

conda activate fgbio


# Define input and output directories
INPUT_DIR=$"/vscratch/grp-vprahlad/umi_analysis/data/output_broad/ExtractUmisFromBam" # Replace with your actual input directory
OUTPUT_DIR=$"/vscratch/grp-vprahlad/umi_analysis/data/output_broad/STARNoClip" # Replace with your actual output directory
# Prepare output directory
if [ ! -d $OUTPUT_DIR ]; then      # if doesn't exit
    mkdir ${OUTPUT_DIR}            # create it
fi


# Create a bash array with sample names
sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3") 
# Add your sample names here



parallel --xapply STAR --runThreadN 8 \
         --genomeDir ../data/output/star_index \
         --readFilesIn  ${INPUT_DIR}/{}.extractUMIs.out.bam \
         --outSAMtype BAM Unsorted \
         --outFileNamePrefix ${OUTPUT_DIR}/{}. \
		 --readFilesType SAM PE \
		 --readFilesCommand samtools view -h \
		 --runRNGseed 12345 \
		 --outSAMunmapped Within \
		 --outFilterType BySJout \
		 --outFilterMultimapNmax 20 \
		 --outFilterScoreMinOverLread 0.33 \
		 --outFilterMatchNminOverLread 0.33 \
		 --outFilterMismatchNmax 999 \
		 --outFilterMismatchNoverLmax 0.1 \
		 --alignIntronMin 20 \
		 --alignIntronMax 1000000 \
		 --alignMatesGapMax 1000000 \
		 --alignSJoverhangMin 8 \
		 --alignSJDBoverhangMin 1 \
		 --alignSoftClipAtReferenceEnds Yes \
		 --chimSegmentMin 15 \
		 --chimMainSegmentMultNmax 1 \
		 --chimOutType WithinBAM SoftClip \
		 --chimOutJunctionFormat 0 \
		 --twopassMode Basic \
		 --quantMode TranscriptomeSAM \
		 --alignEndsProtrude 20 ConcordantPair ::: "${sample[@]}"