#!/bin/bash

#SBATCH --job-name=RNAseq_alignment
#SBATCH --account=indikar1
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --mail-user=jrcwycy@umich.edu


### Script to align RNAseq data and obtain gene counts with HTseq ###


#Define directories
SCRIPT_DIR="script_test"
DATA_DIR="$SCRIPT_DIR/data"
REF_DIR="/nfs/turbo/umms-indikar/shared/projects/RECODE/references"
#REF_DIR="/scratch/indikar_root/indikar1/jrcwycy/RECODE/reference/ensembl"
FASTQC_DIR="$SCRIPT_DIR/fastqc"
ALNS_DIR="$SCRIPT_DIR/alignments"
STATS_DIR="$SCRIPT_DIR/bamstats"
COUNTS_DIR="$SCRIPT_DIR/counts"

#Load tools
module load Bioinformatics
module load fastqc
module load minimap2
module load samtools
module load bamtools
module load subread

#Make directories for output if they don't already exist
mkdir -p "$ALNS_DIR"
mkdir -p "$FASTQC_DIR"
mkdir -p "$COUNTS_DIR"
mkdir -p "$STATS_DIR"



#Perform FastQC on each merged .fastq.gz file (merged by barcode)
#for fastq in "$DATA_DIR"/*.fastq.gz; do

#	fastqc "$fastq" -t 8 --quiet --outdir="$FASTQC_DIR"

#done

#Align each .fastq.gz to the reference genome with minimap2
for fastq in "$DATA_DIR"/*.fastq.gz; do

	#Extract file name without extension
	filename=$(basename -- "$fastq")
	filename_no_ext="${filename%.fastq.gz}"

	#Align to ref 
	minimap2 -ax splice --secondary=no --MD -t 8 "$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" "$fastq" | samtools sort -@ 8 -O bam -o "$ALNS_DIR/${filename_no_ext}.bam"

	#Index sorted bam files
	samtools index "$ALNS_DIR/${filename_no_ext}.bam" "$ALNS_DIR/${filename_no_ext}.bam.bai"

	#Get alignment stats
	bamtools stats -in "$ALNS_DIR/${filename_no_ext}.bam" > "$STATS_DIR/${filename_no_ext}.bamstats.txt"

	#Use htseq-count to get gene counts
	htseq-count -f bam -s no "$ALNS_DIR/${filename_no_ext}.bam" "$REF_DIR/Homo_sapiens.GRCh38.111.gtf.gz" > "$COUNTS_DIR/${filename_no_ext}.counts.txt"

	#FeatureCounts
	featureCounts -L -T 8 -a "$REF_DIR/Homo_sapiens.GRCh38.111.gtf.gz" -o "$COUNTS_DIR/${filename_no_ext}.featcounts.txt" "$ALNS_DIR/${filename_no_ext}.bam"

done







