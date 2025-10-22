#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --job-name="fastqc"
#SBATCH --output=stdout/fastqc.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

conda activate fastqc

# Store input and output directory paths in environmental variables
dir_in=$"../data/seqs"
dir_out=$"/vscratch/grp-vprahlad/umi_analysis/data/output_broad/fastqc"

# Prepare output directory
if [ ! -d $dir_out ]; then      # if doesn't exit
    mkdir ${dir_out}            # create it
fi

# Store the full fastq file paths in an environmental variable
raw_fastq_files=$(ls ${dir_in}/*fastq.gz)

# CORRECT: Run fastqc with multi-threading on 15 fastq files simoultaneously
fastqc "$i" -t 15 -q -o ${dir_out} ${raw_fastq_files}
