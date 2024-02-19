#!/bin/bash

WORKING_DIR=$1
RAW_DATA_PATH=$2

# Define variables for common paths and parameters
TRIMGALORE_PATH="/home/igel/Downloads/TrimGalore-0.6.10/trim_galore"

OUTPUT_PATH="./$WORKING_DIR/1_trimgalore_output"

PAIR_DIR="${OUTPUT_PATH}/Paired"
UNPAIR_DIR="${OUTPUT_PATH}/Unpaired"
LOG_DIR="${OUTPUT_PATH}/TrimLogfiles"

# Create output directories if they don't exist
mkdir -p $PAIR_DIR $UNPAIR_DIR $LOG_DIR

# Loop through all .gz files in the directory to get unique sample names
# Remove the file extension and the last three characters from each file name i.e to remove _R1 or _R2 from each file name
# Sort the names and save them to filenames.txt

# Read the sample names into an array
readarray -t sample_list < ./$WORKING_DIR/filenames.txt


for i in "${sample_list[@]}"; do
  $TRIMGALORE_PATH \
  "${RAW_DATA_PATH}/${i}_R1.fq.gz" \
  "${RAW_DATA_PATH}/${i}_R2.fq.gz" \
  --illumina \
  --cores 8 \
  --paired \
  --length 100 \
  --output_dir $PAIR_DIR \
  --retain_unpaired \
  --gzip
  
  # Rename the output files
  mv ${PAIR_DIR}/${i}_R1_val_1.fq.gz ${PAIR_DIR}/${i}_R1.fq.gz
  mv ${PAIR_DIR}/${i}_R2_val_2.fq.gz ${PAIR_DIR}/${i}_R2.fq.gz
  mv ${PAIR_DIR}/${i}_R1_unpaired_1.fq.gz ${UNPAIR_DIR}/${i}_R1.fq.gz
  mv ${PAIR_DIR}/${i}_R2_unpaired_2.fq.gz ${UNPAIR_DIR}/${i}_R2.fq.gz
  mv ${PAIR_DIR}/${i}_R1.fq.gz_trimming_report.txt ${LOG_DIR}/${i}_R1.log
  mv ${PAIR_DIR}/${i}_R2.fq.gz_trimming_report.txt ${LOG_DIR}/${i}_R2.log
done


echo -e "\nTrimGalore run complete! \n\n"