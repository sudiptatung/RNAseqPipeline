#!/bin/bash

WORKING_DIR=$1
ALIGNER=$2
GTF_URL=$3


OUTPUT_DIR="./$WORKING_DIR/3_featurecounts_output"

BAM_DIR="./$WORKING_DIR/2_${ALIGNER}_output/bams"
GTF_FILE=$(basename "$GTF_URL" .gz)

# Download and unzip the GTF file if it doesn't exist
if [ ! -f "$GTF_FILE" ]; then
  wget $GTF_URL
  gunzip "${GTF_FILE}.gz"
fi


# Create output directory for featureCounts if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Read sample list into an array
readarray -t sample_list < ./$WORKING_DIR/filenames.txt

# Run featureCounts for each sample
for i in "${sample_list[@]}"; do
  bam_file="${BAM_DIR}/${i}.bam"
  
  # Create a directory for each sample
  mkdir -p "${OUTPUT_DIR}/${i}"

  # Define log file path
  log_file="${OUTPUT_DIR}/${i}/featureCounts.log"

  # Run featureCounts to generate abundance counts for each sample
  /home/igel/Downloads/subread-2.0.6-Linux-x86_64/bin/featureCounts -p --countReadPairs -T 30 -a "${GTF_FILE}" -o "${OUTPUT_DIR}/${i}/abundance.tsv" -g transcript_id -s 0 "${bam_file}" > "${log_file}" 2>&1

  # Modify the header of the abundance.tsv file
  # Remove existing header
  sed -i '/^#/d' "${OUTPUT_DIR}/${i}/abundance.tsv"
  # Add new header with 'target_id' and 'est_counts'
  sed -i '1s/Geneid/target_id/' "${OUTPUT_DIR}/${i}/abundance.tsv"
  sed -i '1s|'"${bam_file}"'|est_counts|' "${OUTPUT_DIR}/${i}/abundance.tsv"
done

echo -e "\nfeatureCounts count generation complete! \n\n"