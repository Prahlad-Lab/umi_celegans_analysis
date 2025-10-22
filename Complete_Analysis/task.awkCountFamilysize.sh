#!/bin/bash -l
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name="awkCountFamilysize"
#SBATCH --output=stdout/task.awkCountFamilysize.out
#SBATCH --mail-user=Johnny.CruzCorchado@RoswellPark.org
#SBATCH --mail-type=begin,end
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

# No modules needed for this script

INPUT_DIR="/vscratch/grp-vprahlad/umi_analysis/data/output_broad/UMIAwareDuplicateMarkingGenomeNoClip"
OUTPUT_DIR="/vscratch/grp-vprahlad/umi_analysis/data/output_broad/UMIAwareDuplicateMarkingGenomeNoClip"

sample=("N2.30min.HS.1" "N2.30min.HS.2" "N2.30min.HS.3" "PRDE1.30min.HS.1" "PRDE1.30min.HS.2" "PRDE1.30min.HS.3")

# Prepare output directory
if [ ! -d "$OUTPUT_DIR" ]; then
    mkdir -p "$OUTPUT_DIR"
fi

# Loop through each sample
for i in "${sample[@]}"
do
    echo "Processing sample: $i"

    awk '
        # Set input and output field separators to a tab
        BEGIN { FS="\t"; OFS="\t" }

        # For every line EXCEPT the header (line 1)...
        NR > 1 {
            # Assign columns to variables for clarity
            family_size = $8;  # final_umi_count
            unique_id = $9;    # unique_id
        
            # STEP 1: Process each unique family only once.
            if ( !(unique_id in seen) ) {
            
                # Mark this unique_id as seen
                seen[unique_id] = 1;

                # STEP 2: Build the totals using this family's size.
                total_reads += family_size;
                reads_per_size[family_size] += family_size;
            }
        }
        # After reading the entire file...
        END {
            # Print a header for our output
            print "family_size", "num_reads", "proportion";
        
            # Use a feature of GNU awk to sort the output numerically
            PROCINFO["sorted_in"] = "@ind_num_asc";
        
            # STEP 3: Calculate proportions and print the final report.
            for (size in reads_per_size) {
                # Prevent division by zero if file is empty
                if (total_reads > 0) {
                    proportion = reads_per_size[size] / total_reads;
                } else {
                    proportion = 0;
                }
                print size, reads_per_size[size], proportion;
            }
        }
    ' "${INPUT_DIR}/${i}.grouped.tsv" > "${OUTPUT_DIR}/${i}.family_size_summary.tsv"
done

echo "All samples processed."