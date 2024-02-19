#!/bin/bash

WORKING_DIR=$1
RAW_DATA_PATH=$2


# Define SortMeRNA variables
SORTMERNA_PATH="/home/igel/Downloads/sortmerna-4.3.6-bin/sortmerna" # Replace with the actual path to SortMeRNA
REFDB_PATH="./rRNA_databases_v4"  # Replace with the path to SortMeRNA pre-indexed databases
SORTMERNA_OUTPUT_PATH="./$WORKING_DIR/0_sortmerna_output"
SORTMERNA_FILTERED="${SORTMERNA_OUTPUT_PATH}/rRNA"
SORTMERNA_UNFILTERED="${SORTMERNA_OUTPUT_PATH}/cleanData"
SORTMERNA_LOG_DIR="${SORTMERNA_OUTPUT_PATH}/SortMeRNA_Logfiles"

# Download and unzip the SortMeRNA databases if they don't exist
if [ ! -d "$REFDB_PATH" ]; then
  echo "\nDownloading SortMeRNA databases.\n\n"
  wget https://github.com/biocore/sortmerna/releases/download/v4.3.4/database.tar.gz
  mkdir -p "$REFDB_PATH"
  tar -xvf database.tar.gz -C "$REFDB_PATH"
  rm database.tar.gz
fi

# Create output directories for SortMeRNA if they don't exist
mkdir -p "$SORTMERNA_FILTERED" "$SORTMERNA_UNFILTERED" "$SORTMERNA_LOG_DIR"

for file in ${RAW_DATA_PATH}/*.gz; do
    basename "$file" .fastq.gz | sed 's/...$//'
done | sort -u > ./$WORKING_DIR/filenames.txt

# Read the sample names into an array
readarray -t sample_list < ./$WORKING_DIR/filenames.txt

# Run SortMeRNA for each sample
for i in "${sample_list[@]}"; do
  # Define input paired files
  fwd_file="${RAW_DATA_PATH}/${i}_R1.fastq.gz"
  rev_file="${RAW_DATA_PATH}/${i}_R2.fastq.gz"

  # Run SortMeRNA
  $SORTMERNA_PATH \
    -ref ./rRNA_databases_v4/smr_v4.3_default_db.fasta \
    -ref ./rRNA_databases_v4/smr_v4.3_fast_db.fasta \
    -ref ./rRNA_databases_v4/smr_v4.3_sensitive_db.fasta \
    -ref ./rRNA_databases_v4/smr_v4.3_sensitive_db_rfam_seeds.fasta \
    -reads $fwd_file \
    -reads $rev_file \
    -fastx \
    -aligned $SORTMERNA_FILTERED/${i} \
    -other $SORTMERNA_UNFILTERED/${i} \
    -paired_in Ture \
    --zip-out yes \
    --threads 30
   # Move log files (uncomment and adjust if necessary)
   mv ${SORTMERNA_FILTERED}/*log "${SORTMERNA_LOG_DIR}/${i}_sortmerna.log"
   
   FORWARD_OUTPUT="${SORTMERNA_UNFILTERED}/${i}_R1.fq.gz"
   REVERSE_OUTPUT="${SORTMERNA_UNFILTERED}/${i}_R2.fq.gz"
    
   /home/igel/Downloads/bbmap/reformat.sh in="${SORTMERNA_UNFILTERED}/${i}.fq.gz" out=${FORWARD_OUTPUT} out2=${REVERSE_OUTPUT}
   
   rm "${SORTMERNA_UNFILTERED}/${i}.fq.gz"
   rm ~/sortmerna/run/kvdb/*
done

echo -e "\nSortMeRNA filtering complete! \n\n"