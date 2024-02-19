#!/bin/bash

WORKING_DIR=$1
ALIGNER=$2
GTF_URL=$3

OUTPUT_DIR="./$WORKING_DIR/3_htseq_output"

BAM_DIR="./$WORKING_DIR/2_${ALIGNER}_output/bams"
GTF_FILE=$(basename "$GTF_URL" .gz)


# Download and unzip the GTF file if it doesn't exist
if [ ! -f "$GTF_FILE" ]; then
  wget $GTF_URL
  gunzip "${GTF_FILE}.gz"
fi


# Read sample list into an array
readarray -t sample_list < ./$WORKING_DIR/filenames.txt

# Loop through sample list for generating abundance counts with HTSeq
for i in "${sample_list[@]}"; do
  # Create a directory for each sample in the OUTPUT_DIR
  mkdir -p "${OUTPUT_DIR}/${i}"

  # Extract the base name of the file (without the path and extension)
  bam_file="${BAM_DIR}/${i}.bam"

  # Check if BAM index file exists, if not create it
  if [ ! -f "${bam_file}.bai" ]; then
    #echo "Indexing ${bam_file}"
    samtools index "${bam_file}"
  fi
  # Create abundance.tsv with the header line
  echo -e "target_id\test_counts" > "${OUTPUT_DIR}/${i}/abundance.tsv"

  # Run HTSeq-count to generate abundance counts
  htseq-count -f bam -r pos -s no -i transcript_id -n 30 "${bam_file}" "${GTF_FILE}" >> "${OUTPUT_DIR}/${i}/abundance.tsv"
  # Remove lines starting with '__'
  sed -i '/^__/d' "${OUTPUT_DIR}/${i}/abundance.tsv"
done

echo -e "\nHTSeq count generation complete! \n\n"
