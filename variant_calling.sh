#!/bin/bash

# Set paths
reference="/Volumes/T7/GTMBioinfo/reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa"
aligned="/Volumes/T7/GTMBioinfo/aligned"
variants="/Volumes/T7/GTMBioinfo/variants"
mkdir -p "$variants"

# Sample list (Patient 014)
samples=("C_MB014-2" "C_MB014-3" "FFPE_C_MB014" "PBMC_C_MB014" "PBMC_C_MB014-2")

for sample in "${samples[@]}"; do

  # Input/output files
  input_bam="${aligned}/${sample}.chr12.sorted.bam"
  dedup_bam="${aligned}/${sample}.chr12.dedup.bam"
  metrics="${aligned}/${sample}.metrics.txt"
  unfiltered_vcf="${variants}/${sample}.unfiltered.vcf.gz"
  filtered_vcf="${variants}/${sample}.filtered.vcf.gz"
  kras_vcf="${variants}/${sample}.KRAS.vcf.gz"

  # Mark duplicates
  
  gatk MarkDuplicates \
    -I "$input_bam" \
    -O "$dedup_bam" \
    -M "$metrics" \
    --VALIDATION_STRINGENCY SILENT \
    --REMOVE_DUPLICATES false

  samtools index "$dedup_bam"

  # call Mutect2 (somatic variant caller)
  gatk Mutect2 \
    -R "$reference" \
    -I "$dedup_bam" \
    -tumor "$sample" \
    -O "$unfiltered_vcf"

  # Filter Mutect2 calls
  gatk FilterMutectCalls \
    -R "$reference" \
    -V "$unfiltered_vcf" \
    -O "$filtered_vcf"

  # Select variant of KRAS gene region

  kras_region="12:25245350-25250930"

  gatk SelectVariants \
    -R "$reference" \
    -V "$filtered_vcf" \
    -L "$kras_region" \
    -O "$kras_vcf"

done
