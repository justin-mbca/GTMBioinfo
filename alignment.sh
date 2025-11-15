#!/bin/bash

# Reference genome
reference="/Volumes/T7/GTMBioinfo/reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa"

# Paths
fastq="/Volumes/T7/GTMBioinfo/trimmed"
alignment="/Volumes/T7/GTMBioinfo/aligned"
mkdir -p "$alignment"

# Sample-to-SRR mapping
samples=("C_MB014-2" "C_MB014-3" "FFPE_C_MB014" "PBMC_C_MB014" "PBMC_C_MB014-2")
ids=("SRR21014913" "SRR21014914" "SRR21015062" "SRR21015019" "SRR21015008")

# Loop through samples
for i in "${!samples[@]}"; do
  s="${samples[$i]}"
  id="${ids[$i]}"
  fq1="$fastq/${id}_1.trimmed.fastq"
  fq2="$fastq/${id}_2.trimmed.fastq"
  bam="$alignment/${s}.chr12.sorted.bam"

  # Time the alignment
  start=$(date +%s)

  # Alignment pipeline
  bwa mem -t 8 -R "@RG\tID:$s\tSM:$s\tPL:ILLUMINA" "$reference" "$fq1" "$fq2" | \
    samtools view -bS - | \
    samtools sort -o "$bam"

  # Index the BAM
  samtools index "$bam"
done
