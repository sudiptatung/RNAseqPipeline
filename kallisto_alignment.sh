#!/bin/bash

WORKING_DIR=$1
INPUT_DIR=$2
TRANSCRIPTOME_URL=$3

OUTPUT_DIR="./$WORKING_DIR/3_kallisto_output"

REF_FA=$(basename "$TRANSCRIPTOME_URL" .gz)

INDEX_FILE="${REF_FA}_kallisto.index"

# Download reference FASTA file if it doesn't exist
if [ ! -f "${REF_FA}.gz" ]; then
  wget $TRANSCRIPTOME_URL
fi

# Index reference if it doesn't exist

if [ ! -f ${INDEX_FILE} ]; then

  echo -e "\nInitiating Kallisto indexing.\n\n"
  # Build the index
  kallisto index -i $INDEX_FILE ${REF_FA}.gz
  
  echo -e "\nKallisto indexing is complete! \n\n"

fi


# Define variables for common paths and parameters

THREADS=24

# Read sample list into an array
readarray -t sample_list < ./$WORKING_DIR/filenames.txt

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Loop through sample list for quantification
for i in "${sample_list[@]}"; do
  kallisto quant -i $INDEX_FILE -o ${OUTPUT_DIR}/${i} \
  ${INPUT_DIR}/${i}_R1.fq.gz ${INPUT_DIR}/${i}_R2.fq.gz \
  -t $THREADS &> ${OUTPUT_DIR}/${i}.log
done


echo -e "\nKallisto alignment complete! \n\n"