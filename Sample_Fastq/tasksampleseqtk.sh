#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --job-name="seqtk_sample"
#SBATCH --output=stdout/task.sample.seqtk.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc


conda activate seqtk

# --- Define Input and Output Directories ---
IN_DIR="/projects/rpci/vprahlad/umi_celegans/data/seqs"
OUT_DIR="/projects/rpci/vprahlad/umi_celegans/data/seqs_subset"

# --- Create output directory if it doesn't exist ---
echo "Creating output directory: $OUT_DIR"
mkdir -p "$OUT_DIR"

# --- Change to the input directory to find files ---
echo "Changing to input directory: $IN_DIR"
cd "$IN_DIR"


# Set the number of reads you want
NUM_READS=100000

# Set a random seed (any integer works, '100' is just an example)
SEED=100

echo "Starting subsampling of $NUM_READS reads with seed $SEED..."

# Loop over all R1 files in the current directory (which is now $IN_DIR)
for r1_file in *_R1.fastq.gz; do
    
    # 1. Get the sample prefix by removing "_R1.fastq.gz"
    prefix="${r1_file%_R1.fastq.gz}"
    
    # 2. Define the R2 filename
    r2_file="${prefix}_R2.fastq.gz"
    
   # 4. Check if the R2 file actually exists
    if [ -f "$r2_file" ]; then
        echo "Processing pair: $r1_file and $r2_file"

        # --- Define full output paths ---
        out_r1="$OUT_DIR/${prefix}_R1.subset.fastq.gz"
        out_r2="$OUT_DIR/${prefix}_R2.subset.fastq.gz"
        
        # --- Pipe seqtk output to gzip and write to OUT_DIR ---
        seqtk sample -s$SEED "$r1_file" $NUM_READS | gzip > "$out_r1"
        
        # Run seqtk on R2 with the *identical seed*
        seqtk sample -s$SEED "$r2_file" $NUM_READS | gzip > "$out_r2"
        
    else
        echo "WARNING: Could not find partner file $r2_file for $r1_file"
    fi
done

echo "Subsampling complete. Output files are in $OUT_DIR"
