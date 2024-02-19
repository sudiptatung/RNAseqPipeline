#!/bin/bash

WORKING_DIR=$1
INPUT_DIR=$2
TRANSCRIPTOME_URL=$3

OUTPUT_DIR="./$WORKING_DIR/3_salmon_output"

# Define variables for common paths and parameters
SALMON_PATH="salmon"  # Path to the Salmon executable
REF_FA=$(basename "$TRANSCRIPTOME_URL" .gz)

INDEX_DIR="${REF_FA}_salman.index"

# Download reference FASTA file if it doesn't exist
if [ ! -f "${REF_FA}.gz" ]; then
  wget $TRANSCRIPTOME_URL
fi

# Index reference if it doesn't exist

if [ ! -d $INDEX_DIR ]; then

  echo -e "\nInitiating Salman reference indexing.\n\n"

  # Create output directory if it doesn't exist
  mkdir -p $INDEX_DIR

  # Index the reference genome
  $SALMON_PATH index \
  -t ${REF_FA}.gz \
  -i $INDEX_DIR \
  -p 24
  
  echo -e "\nSalman reference indexing run complete! \n\n"

fi


# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Read sample list from a file
readarray -t sample_list < ./$WORKING_DIR/filenames.txt

# Loop through sample list for alignment
for i in "${sample_list[@]}"; do
  $SALMON_PATH quant \
  -i $INDEX_DIR \
  -l A \
  -o "${OUTPUT_DIR}/$i" \
  -1 "${INPUT_DIR}/${i}_R1.fq.gz" \
  -2 "${INPUT_DIR}/${i}_R2.fq.gz" \
  -p 24 \
  --validateMappings \
  &> "${OUTPUT_DIR}/${i}.log"

  # Process the output files
  awk 'NR==1 {print "target_id\tlength\teff_length\test_counts\ttpm"} NR>1 {print $1 "\t" $2 "\t" $3 "\t" $5 "\t" $4}' "${OUTPUT_DIR}/${i}/quant.sf" > "${OUTPUT_DIR}/${i}/abundance.tsv"
done


echo -e "\nSalman alignment run complete! \n\n"