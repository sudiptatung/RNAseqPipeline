#!/bin/bash

WORKING_DIR=$1
INPUT_DIR=$2
GENOME_URL=$3
GTF_URL=$4

OUTPUT_DIR="./$WORKING_DIR/2_hisat2_output"

# Define variables for common paths and parameters
THREADS=24

GENOME_FASTA=$(basename "$GENOME_URL" .gz)

INDEX_DIR="${GENOME_FASTA}_hisat2.index"

# Check for existence of the reference genome file
if [ ! -f "${GENOME_FASTA}" ]; then
  wget ${GENOME_URL}
  gunzip "${GENOME_FASTA}.gz"
fi


if [ ! -d "${INDEX_DIR}" ]; then
  echo -e "\nInitiating HISAT2 reference indexing...\n\n"
  # Create directory for HISAT2 index if not exists
  mkdir -p ${INDEX_DIR}

  # Indexing
  hisat2-build ${GENOME_FASTA} ${INDEX_DIR}/index -p ${THREADS}
  echo -e "\nHISAT2 reference indexing is complete! \n\n"
fi


# Read sample list into an array
readarray -t sample_list < ./$WORKING_DIR/filenames.txt

# Create output directories if they don't exist
mkdir -p ${OUTPUT_DIR}/bams

# Loop through sample list for alignment
for i in "${sample_list[@]}"; do
  hisat2 -p ${THREADS} -x ${INDEX_DIR}/index \
  -1 ${INPUT_DIR}/${i}_R1.fq.gz \
  -2 ${INPUT_DIR}/${i}_R2.fq.gz \
  --dta --summary-file ${OUTPUT_DIR}/${i}_summary.txt | \
  samtools sort -@ ${THREADS} -o ${OUTPUT_DIR}/bams/${i}.bam -
done


echo -e "\nHISAT2 alignment is complete! \n\n"