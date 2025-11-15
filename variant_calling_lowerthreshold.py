#!/bin/bash

# Set paths
reference="/Volumes/T7/GTMBioinfo/reference/Homo_sapiens.GRCh38.dna.chromosome.12.fa"
aligned="/Volumes/T7/GTMBioinfo/aligned"
variants="/Volumes/T7/GTMBioinfo/variants_low_threshold"
mkdir -p "$variants"

# Sample list (Patient 014)
samples=("C_MB014-2" "C_MB014-3" "FFPE_C_MB014" "PBMC_C_MB014" "PBMC_C_MB014-2")

for sample in "${samples[@]}"; do
  echo "----------------------------------------"
  echo "🔬 Processing $sample with relaxed thresholds..."

  # Input/output files
  input_bam="${aligned}/${sample}.chr12.sorted.bam"
  dedup_bam="${variants}/${sample}.chr12.dedup.relaxed.bam"
  metrics="${variants}/${sample}.metrics.relaxed.txt"
  unfiltered_vcf="${variants}/${sample}.unfiltered.relaxed.vcf.gz"
  filtered_vcf="${variants}/${sample}.filtered.relaxed.vcf.gz"
  kras_vcf="${variants}/${sample}.KRAS.relaxed.vcf.gz"
  kras_region="12:25245350-25250930"

  # Check input BAM
  if [[ ! -f "$input_bam" ]]; then
    echo "Missing input BAM for $sample: $input_bam"
    continue
  fi

  # Time the run
  start=$(date +%s)

  # Step 1: Mark duplicates
  gatk MarkDuplicates \
    -I "$input_bam" \
    -O "$dedup_bam" \
    -M "$metrics" \
    --VALIDATION_STRINGENCY SILENT \
    --REMOVE_DUPLICATES false 

  samtools index "$dedup_bam" 

  # Step 2: Mutect2 with relaxed thresholds
  gatk Mutect2 \
    -R "$reference" \
    -I "$dedup_bam" \
    -tumor "$sample" \
    --initial-tumor-lod 2.0 \
    --tumor-lod-to-emit 2.0 \
    --max-reads-per-alignment-start 100 \
    --min-base-quality-score 20 \
    -O "$unfiltered_vcf" 

  # Step 3: Filter Mutect2 calls
  gatk FilterMutectCalls \
    -R "$reference" \
    -V "$unfiltered_vcf" \
    --min-median-mapping-quality 20 \
    --min-median-base-quality 20 \
    -O "$filtered_vcf" 

  # Step 4: Select KRAS region
  gatk SelectVariants \
    -R "$reference" \
    -V "$filtered_vcf" \
    -L "$kras_region" \
    -O "$kras_vcf" 
  
done

