#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --job-name="Combined_UMI_Pipeline"
#SBATCH --output=../Combined_UMI_Pipeline.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc
#SBATCH --mem=100G

#==============================================================================
# START TIMER
#==============================================================================
START_TIME=$SECONDS
echo "--- Pipeline started at $(date) ---"

#==============================================================================
# DEFINE VARIABLES AND DIRECTORIES
#==============================================================================

# --- Input/Reference ---
DIR_SEQS="../Sample_Fastq"
DIR_REF_INPUT="../input"
REF_GENOME="${DIR_REF_INPUT}/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa"
REF_DICT="${DIR_REF_INPUT}/Caenorhabditis_elegans.WBcel235.dna.toplevel.dict"
REF_GTF_FILE="${DIR_REF_INPUT}/Caenorhabditis_elegans.WBcel235.114.gtf "
KNOWN_SITES_VCF="${DIR_REF_INPUT}/random_subset.vcf.gz"
STAR_INDEX="../star_index" # This directory will be created
PYTHON_SCRIPT_PATH="../Separate_Scripts/calculate_allele_proportions_depth.py"

# --- Output Directories ---
BASE_OUTPUT_DIR="../output_subset"
DIR_FASTQC="${BASE_OUTPUT_DIR}/fastqc"
DIR_FASTQ_TO_UBAM="${BASE_OUTPUT_DIR}/FastqToUbam"
DIR_EXTRACT_UMIS="${BASE_OUTPUT_DIR}/ExtractUmisFromBam"
DIR_STAR_ALIGN="${BASE_OUTPUT_DIR}/STARNoClip"
DIR_UMI_MARKING="${BASE_OUTPUT_DIR}/UMIAwareDuplicateMarkingGenomeNoClip"
DIR_FILTER_BAM="${BASE_OUTPUT_DIR}/FilterBambyFamilySizeGenome"
DIR_HAPLOTYPE_CALLER="${BASE_OUTPUT_DIR}/Family_size1_haplotype_caller"
DIR_ANNOTATION="${BASE_OUTPUT_DIR}/Family_size1_annotation"
DIR_BCFTOOLS_PILEUP="${BASE_OUTPUT_DIR}/Family_size1_bcftools_mpileup"
DIR_ALLELE_PROPORTIONS="${BASE_OUTPUT_DIR}/Family_size1_allele_proportions_depth"


# --- Sample Array ---
sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3")

# --- Create all output directories ---
mkdir -p ${DIR_FASTQC}
mkdir -p ${DIR_FASTQ_TO_UBAM}
mkdir -p ${DIR_EXTRACT_UMIS}
mkdir -p ${DIR_STAR_ALIGN}
mkdir -p ${DIR_UMI_MARKING}
mkdir -p ${DIR_FILTER_BAM}
mkdir -p ${DIR_HAPLOTYPE_CALLER}
mkdir -p ${DIR_ANNOTATION}
mkdir -p ${DIR_BCFTOOLS_PILEUP}
mkdir -p ${DIR_ALLELE_PROPORTIONS}
mkdir -p stdout # For SBATCH logs

#==============================================================================
# ACTIVATE MASTER CONDA ENVIRONMENT
#==============================================================================
echo "--- Activating Conda Environment: umi_celegans_analysis ---"
conda activate umi_celegans_analysis
echo "--- Conda environment activated. ---"

#==============================================================================
# STEP 1: FASTQC (Quality Control)
#==============================================================================
echo "--- Starting Step 1: FastQC ---"

# Store the full fastq file paths in an environmental variable
raw_fastq_files=$(ls ${DIR_SEQS}/*fastq.gz)

# Run fastqc with multi-threading
fastqc -t 15 -q -o ${DIR_FASTQC} ${raw_fastq_files}

echo "--- Finished Step 1: FastQC ---"

#==============================================================================
# STEP 2: FastqToUbam (Convert FastQ to unmapped BAM)
#==============================================================================
echo "--- Starting Step 2: FastqToUbam ---"

parallel --xapply gatk FastqToSam \
      SORT_ORDER=unsorted \
      F1=${DIR_SEQS}/{}_R1.subset.fastq.gz  \
      F2=${DIR_SEQS}/{}_R2.subset.fastq.gz \
      SM={} \
      LB="L1" \
      PL="ILLUMINA" \
      PU="barcode1" \
      RG="RG1" \
      CN="BI" \
      O=${DIR_FASTQ_TO_UBAM}/{}.unmapped.bam \
      ALLOW_EMPTY_FASTQ=True ::: "${sample[@]}"

echo "--- Finished Step 2: FastqToUbam ---"

#==============================================================================
# STEP 3: ExtractUMIs (Extract UMIs from uBAM)
#==============================================================================
echo "--- Starting Step 3: ExtractUMIs ---"

parallel --xapply fgbio ExtractUmisFromBam \
      --input ${DIR_FASTQ_TO_UBAM}/{}.unmapped.bam \
      --read-structure +T \
      --read-structure 8M6S+T \
      --molecular-index-tags RX \
      --output ${DIR_EXTRACT_UMIS}/{}.extractUMIs.out.bam ::: "${sample[@]}"

echo "--- Finished Step 3: ExtractUMIs ---"

#==============================================================================
# STEP 4: Sort UMI-extracted BAM by QueryName
#==============================================================================
echo "--- Starting Step 4: Sort UMI-extracted BAM by QueryName ---"

parallel --xapply gatk SortSam \
      INPUT=${DIR_EXTRACT_UMIS}/{}.extractUMIs.out.bam \
      OUTPUT=${DIR_EXTRACT_UMIS}/{}.querysort.extractUMIs.out.bam \
      SORT_ORDER="queryname" \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000 ::: "${sample[@]}"

echo "--- Finished Step 4: Sort UMI-extracted BAM by QueryName ---"

#==============================================================================
# STEP 4.5: Generate STAR Index 
#==============================================================================
echo "--- Starting Step 4.5: Generate STAR Index ---"

# Ensure the STAR_INDEX directory exists and is empty, or create it
if [ -d "${STAR_INDEX}" ]; then
    echo "STAR index directory ${STAR_INDEX} already exists. Removing for fresh build."
    rm -rf "${STAR_INDEX}"
fi
mkdir -p ${STAR_INDEX}

echo "Generating STAR index in ${STAR_INDEX} using ${REF_GENOME}"
STAR --runThreadN ${SLURM_CPUS_PER_TASK:-32} \
     --runMode genomeGenerate \
     --genomeDir ${STAR_INDEX} \
     --genomeFastaFiles ${REF_GENOME}\
     --sjdbGTFfile ${REF_GTF_FILE} \
     --sjdbOverhang 100 # This is often (ReadLength - 1)

echo "--- Finished Step 4.5: Generate STAR Index ---"


#==============================================================================
# STEP 5: STAR Alignment
#==============================================================================
echo "--- Starting Step 5: STAR Alignment ---"

parallel --xapply STAR --runThreadN 8 \
         --genomeDir ${STAR_INDEX} \
         --readFilesIn ${DIR_EXTRACT_UMIS}/{}.extractUMIs.out.bam \
         --outSAMtype BAM Unsorted \
         --outFileNamePrefix ${DIR_STAR_ALIGN}/{}. \
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

echo "--- Finished Step 5: STAR Alignment ---"

#==============================================================================
# STEP 6: Sort Aligned BAM by QueryName
#==============================================================================
echo "--- Starting Step 6: Sort Aligned BAM by QueryName ---"

parallel --xapply gatk SortSam \
      INPUT=${DIR_STAR_ALIGN}/{}.Aligned.out.bam \
      OUTPUT=${DIR_UMI_MARKING}/{}.querysort.Aligned.out.bam \
      SORT_ORDER="queryname" \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000 ::: "${sample[@]}"

echo "--- Finished Step 6: Sort Aligned BAM by QueryName ---"

#==============================================================================
# STEP 7: MergeBamAlignment (Merge aligned BAM with unmapped BAM)
#==============================================================================
echo "--- Starting Step 7: MergeBamAlignment ---"

parallel --xapply gatk MergeBamAlignment \
            --REFERENCE_SEQUENCE ${REF_GENOME} \
            --UNMAPPED_BAM ${DIR_EXTRACT_UMIS}/{}.querysort.extractUMIs.out.bam \
            --ALIGNED_BAM ${DIR_UMI_MARKING}/{}.querysort.Aligned.out.bam \
            --OUTPUT ${DIR_UMI_MARKING}/{}.MergeBamAlignment.bam \
            --INCLUDE_SECONDARY_ALIGNMENTS false \
            --VALIDATION_STRINGENCY SILENT ::: "${sample[@]}"

echo "--- Finished Step 7: MergeBamAlignment ---"

#==============================================================================
# STEP 8: Sort Merged BAM by Coordinate
#==============================================================================
echo "--- Starting Step 8: Sort Merged BAM by Coordinate ---"

parallel --xapply gatk SortSam \
      INPUT=${DIR_UMI_MARKING}/{}.MergeBamAlignment.bam \
      OUTPUT=${DIR_UMI_MARKING}/{}.MergeBamAlignment.coordinate_sorted.bam \
      SORT_ORDER="coordinate" \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000 ::: "${sample[@]}"

echo "--- Finished Step 8: Sort Merged BAM by Coordinate ---"

#==============================================================================
# STEP 9: GroupByUMIs (Group reads by UMI)
#==============================================================================
echo "--- Starting Step 9: GroupByUMIs ---"

export PYTHONHASHSEED=0

parallel --xapply umi_tools group \
    -I ${DIR_UMI_MARKING}/{}.MergeBamAlignment.coordinate_sorted.bam \
    --paired \
    --no-sort-output \
    --output-bam \
    --stdout ${DIR_UMI_MARKING}/{}.grouped_by_UMI.bam \
    --umi-tag-delimiter "-" \
    --log ${DIR_UMI_MARKING}/{}.group.log \
    --group-out ${DIR_UMI_MARKING}/{}.grouped.tsv \
    --extract-umi-method tag \
    --umi-tag RX \
    --unmapped-reads use ::: "${sample[@]}"

echo "--- Finished Step 9: GroupByUMIs ---"

#==============================================================================
# STEP 10: awkCountFamilysize (Generate family size report)
#==============================================================================
echo "--- Starting Step 10: Count Family Sizes ---"
# No conda env needed for awk

for i in "${sample[@]}"
do
    echo "Processing family size summary for sample: $i"

    awk '
        BEGIN { FS="\t"; OFS="\t" }
        NR > 1 {
            family_size = $8;
            unique_id = $9;
            if ( !(unique_id in seen) ) {
                seen[unique_id] = 1;
                total_reads += family_size;
                reads_per_size[family_size] += family_size;
            }
        }
        END {
            print "family_size", "num_reads", "proportion";
            PROCINFO["sorted_in"] = "@ind_num_asc";
            for (size in reads_per_size) {
                if (total_reads > 0) {
                    proportion = reads_per_size[size] / total_reads;
                } else {
                    proportion = 0;
                }
                print size, reads_per_size[size], proportion;
            }
        }
    ' "${DIR_UMI_MARKING}/${i}.grouped.tsv" > "${DIR_UMI_MARKING}/${i}.family_size_summary.tsv"
done

echo "--- Finished Step 10: Count Family Sizes ---"

#==============================================================================
# STEP 11: FilterBambyFamilySize (Create BAMs for singletons)
#==============================================================================
echo "--- Starting Step 11: Filter BAMs by Family Size ---"

# --- Part A: Process Singletons (Family Size == 1) ---
echo "--- Step 11a: Extracting Singleton IDs ---"
for i in "${sample[@]}"
do
    echo "Processing singletons for sample: $i"
	awk 'BEGIN{FS="\t"} NR > 1 && $8 == 1 {print $9}' "${DIR_UMI_MARKING}/${i}.grouped.tsv" | sort -u > "${DIR_FILTER_BAM}/${i}.singletons_id.txt"
done;


echo "--- Step 11b: Creating Singleton BAMs ---"
for i in "${sample[@]}"
do
    echo "Filtering BAM for singletons: $i"
	samtools view -@ 8 -h --tag-file UG:"${DIR_FILTER_BAM}/${i}.singletons_id.txt" -o "${DIR_FILTER_BAM}/${i}.singletons.bam" "${DIR_UMI_MARKING}/${i}.grouped_by_UMI.bam"
done;

echo "--- Step 11c: Sorting and Indexing Singleton BAMs ---"
for i in "${sample[@]}"
do
    echo "Sorting/Indexing singletons: $i"
	samtools sort -@ 8 -o "${DIR_FILTER_BAM}/${i}.sort_coood_singletons.bam" "${DIR_FILTER_BAM}/${i}.singletons.bam";
	samtools index -@ 8 "${DIR_FILTER_BAM}/${i}.sort_coood_singletons.bam" 
done;



echo "--- Finished Step 11: Filter BAMs by Family Size ---"

#==============================================================================
# STEP 12: GATK SplitNCigarReads for RNA seq + BQSR for recalibration
#==============================================================================
echo "--- Starting Step 12: GATK BQSR Pipeline (SplitN, BQSR, ApplyBQSR) ---"

# --- Step 12a: SplitNCigarReads ---
echo "--- Step 12a: SplitNCigarReads for all samples ---"
for i in ${sample[@]}
do
    echo "Step 12a - Processing sample: ${i}"
    gatk SplitNCigarReads \
        -R ${REF_GENOME} \
        -I ${DIR_FILTER_BAM}/"${i}".sort_coood_singletons.bam \
        -O ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.ready.haplotype.caller.bam
done

# --- Step 12b: BaseRecalibrator ---
echo "--- Step 12b: BaseRecalibrator for all samples ---"
for i in ${sample[@]}
do
    echo "Step 12b - Processing sample: ${i}"
    gatk BaseRecalibrator \
        -R ${REF_GENOME} \
        -I ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.ready.haplotype.caller.bam \
        --use-original-qualities \
        -O ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.recal.table \
        --known-sites ${KNOWN_SITES_VCF}
done

# --- Step 12c: ApplyBQSR ---
echo "--- Step 12c: ApplyBQSR for all samples ---"
for i in ${sample[@]}
do
    echo "Step 12c - Processing sample: ${i}"
    gatk ApplyBQSR \
        -R ${REF_GENOME} \
        -I ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.ready.haplotype.caller.bam \
        -O ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.BSQR.bam \
        --bqsr-recal-file ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.recal.table \
        -OBI \
        --sequence-dictionary ${REF_DICT}
done
echo "--- Finished Step 12: GATK BQSR Pipeline ---"


#==============================================================================
# --- PIPELINE FORK ---
# The output from Step 12 ("${i}".singletons.BSQR.bam) is now used for
# two separate, parallel analyses.
#
# PATH A: GATK Variant Calling and SnpEff Annotation (Step 13)
# PATH B: bcftools mpileup and Allele Proportion Calculation (Steps 14 & 15)
#==============================================================================


#==============================================================================
# PATH A - STEP 13: Variant Calling (HaplotypeCaller) & Variant Annotation (SnpEff)
#==============================================================================
echo "--- Starting PATH A, Step 13: Variant Calling and Annotation ---"

# --- Step 13a: HaplotypeCaller ---
echo "--- Step 13a: HaplotypeCaller for all samples ---"
for i in ${sample[@]}
do
    echo "Step 13a - Processing sample: ${i}"
    gatk HaplotypeCaller \
        -R ${REF_GENOME} \
        -I ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.BSQR.bam \
        --dont-use-soft-clipped-bases true \
        -O ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.first.vcf.gz \
        --standard-min-confidence-threshold-for-calling 20 \
        --dbsnp ${KNOWN_SITES_VCF}
done

# --- Step 13b: SnpEff Annotation ---
echo "--- Step 13b: SnpEff Annotation for all samples ---"
for i in ${sample[@]}
do
    echo "Step 13b - Annotating sample: ${i}"
    # Annotate the VCF file
    snpEff WBcel235.105 -stats ${DIR_ANNOTATION}/"${i}"_singleton_summary.html ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.first.vcf.gz > ${DIR_ANNOTATION}/"${i}".singletons.ann.vcf

    # Generate additional BED annotation and CSV stats
    snpEff ann WBcel235.105 -csvStats ${DIR_ANNOTATION}/"${i}"_singletons_summary.csv -o bedANN ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.first.vcf.gz > ${DIR_ANNOTATION}/"${i}".singletons.ann.bed
done
echo "--- Finished PATH A, Step 13 ---"


#==============================================================================
# PATH B - STEP 14: bcftools mpileup 
#==============================================================================
echo "--- Starting PATH B, Step 14: bcftools mpileup ---"

for i in ${sample[@]}
do
    echo "Step 14 - Processing sample: ${i}"
    bcftools mpileup \
        -d 1000000 \
        -f ${REF_GENOME} \
        ${DIR_HAPLOTYPE_CALLER}/"${i}".singletons.BSQR.bam > ${DIR_BCFTOOLS_PILEUP}/"${i}".singletons.bcftools.pileup 2>&1 &
done

wait # Wait for all background mpileup processes to finish
echo "--- Finished PATH B, Step 14 ---"


#==============================================================================
# PATH B - STEP 15: Calculate Allele Proportions
#==============================================================================
echo "--- Starting PATH B, Step 15: Calculate Allele Proportions ---"
 # To use python 3


for i in ${sample[@]}
do
    echo "Step 15 - Processing sample: ${i}"
    python ${PYTHON_SCRIPT_PATH} \
        ${DIR_BCFTOOLS_PILEUP}/"${i}".singletons.bcftools.pileup \
        ${DIR_ALLELE_PROPORTIONS}/"${i}".singletons.per_position_results.tsv \
        ${DIR_ALLELE_PROPORTIONS}/"${i}".singletons.summary_totals.txt --min-depth 2
done
echo "--- Finished PATH B, Step 15 ---"


echo "--- All pipeline steps complete. ---"

#==============================================================================
# CALCULATE AND SAVE TOTAL RUNTIME
#==============================================================================
END_TIME=$SECONDS
DURATION=$(( END_TIME - START_TIME ))
RUNTIME_FILE="${BASE_OUTPUT_DIR}/pipeline_runtime.txt"

echo "--- Pipeline execution finished at $(date) ---"
echo "Total runtime: ${DURATION} seconds"
echo "Total runtime (HH:MM:SS): $(printf '%02d:%02d:%02d\n' $((DURATION/3600)) $((DURATION%3600/60)) $((DURATION%60)))" | tee ${RUNTIME_FILE}

echo "Runtime information saved to ${RUNTIME_FILE}"
