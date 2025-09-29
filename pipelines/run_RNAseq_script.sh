#!/bin/bash

#SBATCH --job-name=Run_RNAseq_alignment
#SBATCH --account=indikar1
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --time=1:00:00
#SBATCH --mail-user=jrcwycy@umich.edu
#SBATCH --mail-type=FAIL,END


# Directory containing .fastq files
FASTQ_DIR="/scratch/indikar_root/indikar1/jrcwycy/RECODE/script_test/data"

# Check if FASTQ_DIR exists and is a directory
if [ ! -d "$FASTQ_DIR" ]; then
    echo "Error: $FASTQ_DIR is not a directory or does not exist."
    exit 1
fi

# Loop through each .fastq file in the directory
for fastq_file in "$FASTQ_DIR"/*.fastq.gz; do
    # Check if the file is a regular file
    if [ -f "$fastq_file" ]; then
        # Run RNAseq_script.sh on the file
        sbatch RNAseq_script.sh "$fastq_file"
    fi
done