#!/bin/bash

WORKING_DIR=$1
RAW_DATA_PATH=$2

# Define variables for common paths and parameters
TRIMMOMATIC_PATH="/home/igel/Downloads/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS="/home/igel/Downloads/Trimmomatic-0.39/adapters/TruSeq3-PE.fa"

OUTPUT_PATH="./$WORKING_DIR/1_trimmomatic_output"

# Create output directories if they don't exist
mkdir -p ${OUTPUT_PATH}/{TrimLogfiles,Paired,Unpaired}

# Loop through all .gz files in the directory to get unique sample names
# Remove the file extension and the last three characters from each file name

# Read the sample names into an array
readarray -t sample_list < ./$WORKING_DIR/filenames.txt

for i in "${sample_list[@]}"; do
  java -jar $TRIMMOMATIC_PATH PE -threads 24 \
  -trimlog "${OUTPUT_PATH}/TrimLogfiles/${i}.log" \
  "${RAW_DATA_PATH}/${i}_R1.fq.gz" \
  "${RAW_DATA_PATH}/${i}_R2.fq.gz" \
  "${OUTPUT_PATH}/Paired/${i}_R1.fq.gz" \
  "${OUTPUT_PATH}/Unpaired/${i}_R1.fq.gz" \
  "${OUTPUT_PATH}/Paired/${i}_R2.fq.gz" \
  "${OUTPUT_PATH}/Unpaired/${i}_R2.fq.gz" \
  ILLUMINACLIP:${ADAPTERS}:2:30:10:2:True \
  LEADING:3 TRAILING:3 MINLEN:100
done


echo -e "\nTrimmomatic trimming complete! \n\n"