#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name="bctools"
#SBATCH --output=bctools.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc



conda activate bcftools

# --- Configuration ---
# The input VCF file
INPUT_VCF="CB4856.hard-filter.vcf.gz"

# The desired number of random variants
NUM_VARIANTS=1000

# The name for the output subset file (will be compressed)
OUTPUT_VCF_UNCOMPRESSED="random_subset.vcf"
OUTPUT_VCF_COMPRESSED="${OUTPUT_VCF_UNCOMPRESSED}.gz"
# ---------------------

echo "Starting VCF subsetting process..."

# Check if required commands are available
command -v bcftools >/dev/null 2>&1 || { echo >&2 "Error: 'bcftools' is not installed or not in PATH."; exit 1; }
command -v shuf >/dev/null 2>&1 || { echo >&2 "Error: 'shuf' (from coreutils) is not installed or not in PATH."; exit 1; }
command -v bgzip >/dev/null 2>&1 || { echo >&2 "Error: 'bgzip' (from htslib) is not installed or not in PATH."; exit 1; }

# Check if input file exists
if [ ! -f "$INPUT_VCF" ]; then
    echo "Error: Input file '$INPUT_VCF' not found."
    exit 1
fi

# 1. Extract the VCF header and write it to the new file
#    (This overwrites the output file if it already exists)
echo "Step 1: Extracting header from $INPUT_VCF..."
bcftools view -h "$INPUT_VCF" > "$OUTPUT_VCF_UNCOMPRESSED"
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract header."
    exit 1
fi

# 2. Extract only the variant lines (no header),
#    pipe them to 'shuf' to select N random lines,
#    then sort them by chromosome (col 1, version sort) and position (col 2, numeric sort)
#    and append them to the file containing the header.
echo "Step 2: Sampling $NUM_VARIANTS random variants and sorting them..."
bcftools view -H "$INPUT_VCF" | shuf -n "$NUM_VARIANTS" | sort -k1,1V -k2,2n >> "$OUTPUT_VCF_UNCOMPRESSED"
if [ $? -ne 0 ]; then
    echo "Error: Failed to sample variants. This could be a pipe error or 'shuf' error."
    # Clean up the partial file
    rm "$OUTPUT_VCF_UNCOMPRESSED"
    exit 1
fi

echo "Successfully created uncompressed subset: $OUTPUT_VCF_UNCOMPRESSED"

# 3. Compress the resulting VCF file with bgzip
#    -f forces overwrite if the .gz file already exists
echo "Step 3: Compressing file with bgzip..."
bgzip -f "$OUTPUT_VCF_UNCOMPRESSED"
if [ $? -ne 0 ]; then
    echo "Error: bgzip compression failed."
    exit 1
fi

echo "Successfully compressed file: $OUTPUT_VCF_COMPRESSED"

# 4. Index the compressed VCF file with bcftools
#    -t for tabix index (.tbi)
echo "Step 4: Indexing file with bcftools..."
bcftools index -t "$OUTPUT_VCF_COMPRESSED"
if [ $? -ne 0 ]; then
    echo "Error: bcftools indexing failed."
    exit 1
fi

echo "Successfully indexed file: ${OUTPUT_VCF_COMPRESSED}.tbi"
echo "---"
echo "Process complete! Your file is ready at $OUTPUT_VCF_COMPRESSED"