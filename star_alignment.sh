#!/bin/bash

WORKING_DIR=$1
INPUT_DIR=$2
GENOME_URL=$3
GTF_URL=$4

OUTPUT_DIR="./$WORKING_DIR/2_star_output"

# Local file names
GENOME_FASTA=$(basename "$GENOME_URL" .gz)
GTF_FILE=$(basename "$GTF_URL" .gz)

# Download and unzip the reference genome if it doesn't exist
if [ ! -f "$GENOME_FASTA" ]; then
  wget $GENOME_URL
  gunzip "${GENOME_FASTA}.gz"
fi

# Download and unzip the GTF file if it doesn't exist
if [ ! -f "$GTF_FILE" ]; then
  wget $GTF_URL
  gunzip "${GTF_FILE}.gz"
fi

# Define variables for common paths and parameters
STAR_PATH="STAR"  # Path to the STAR executable, adjust if needed
INDEX_DIR="${GENOME_FASTA}_star.index/"

if [ ! -d "$INDEX_DIR" ]; then
  echo -e "\nInitiating STAR reference indexing! \n\n"
  # Create the index directory if it doesn't exist
  mkdir -p $INDEX_DIR

  # Run STAR for genome indexing
  $STAR_PATH \
  --runMode genomeGenerate \
  --genomeDir $INDEX_DIR \
  --genomeFastaFiles $GENOME_FASTA \
  --sjdbGTFfile $GTF_FILE \
  --runThreadN 24
  echo -e "\nSTAR reference indexing is complete! \n\n"
fi

# Create output directories if they don't exist
mkdir -p "${OUTPUT_DIR}/bams"

# Read sample list into an array
readarray -t sample_list < ./$WORKING_DIR/filenames.txt

# Loop through sample list for alignment
for i in "${sample_list[@]}"; do
  STAR --runMode alignReads \
  --runThreadN 20 \
  --genomeDir ${INDEX_DIR} \
  --outSAMtype BAM SortedByCoordinate \
  --readFilesIn "${INPUT_DIR}/${i}_R1.fq.gz" "${INPUT_DIR}/${i}_R2.fq.gz" \
  --outFileNamePrefix "${OUTPUT_DIR}/${i}_" \
  --readFilesCommand gunzip -c
  
  # Rename the output .bam file to match ${i}.bam
  mv "${OUTPUT_DIR}/${i}_Aligned.sortedByCoord.out.bam" "${OUTPUT_DIR}/bams/${i}.bam"
done

echo -e "\nSTAR alignment complete! \n\n"
