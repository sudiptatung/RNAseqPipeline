#!/bin/bash

# Assign the arguments to variables
TRIMMER=$1
ALIGNER=$2
COUNTER=$3


# Print the argument
echo -e "\nThe provided arguments are: $TRIMMER $ALIGNER $COUNTER\n\n"

WORKING_DIR=0.$TRIMMER-$ALIGNER-$COUNTER
mkdir -p $WORKING_DIR

RAW_DATA_PATH="/home/igel/Documents/users/Sudipta/Transcriptomics/Chandrakanth_analysis/Jan2024/IGELpipeline/data"

GENOME_URL="http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.55_FB2023_06/fasta/dmel-all-chromosome-r6.55.fasta.gz"
GTF_URL="http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.55_FB2023_06/gtf/dmel-all-r6.55.gtf.gz"
TRANSCRIPTOME_URL="http://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r6.55_FB2023_06/fasta/dmel-all-transcript-r6.55.fasta.gz"

mkdir -p $WORKING_DIR/QC_reports/fastqc_RAWDATA

##QC_RAWDATA
fastqc -t 30 -o $WORKING_DIR/QC_reports/fastqc_RAWDATA $RAW_DATA_PATH/*gz
multiqc -o $WORKING_DIR/QC_reports -n 1.multiqc_RAWDATA.html $WORKING_DIR/QC_reports/fastqc_RAWDATA


#---------- rRNA cleanup ------------

echo -e "\nStarting SortMeRNA filtering.\n\n"
bash sortmerna.sh $WORKING_DIR $RAW_DATA_PATH


multiqc -o $WORKING_DIR/QC_reports -n 2.multiqc_rRNA_status.html $WORKING_DIR/0_sortmerna_output/SortMeRNA_Logfiles


#---------- TRIMMING ----------------

CLEAN_DATA_PATH="./$WORKING_DIR/0_sortmerna_output/cleanData"

# Check if an argument was provided
if [ "$TRIMMER" == "trimgalore" ] || [ "$TRIMMER" == "trimmomatic" ]; then
    # Argument is not null, run trimming script
    bash "${TRIMMER}.sh" $WORKING_DIR $CLEAN_DATA_PATH
    ALIGNMENT_INPUT="./$WORKING_DIR/1_${TRIMMER}_output/Paired"
    
else
    echo -e "\nTrimming is skipped.\n\n"
    ALIGNMENT_INPUT=$CLEAN_DATA_PATH
    
fi

mkdir -p $WORKING_DIR/QC_reports/fastqc_preAlignment
fastqc -t 30 -o $WORKING_DIR/QC_reports/fastqc_preAlignment $ALIGNMENT_INPUT/*gz
multiqc -o $WORKING_DIR/QC_reports -n 3.multiqc_postTrimming.html $WORKING_DIR/QC_reports/fastqc_preAlignment

#---------- ALIGNMENT ---------------

if [ "$ALIGNER" == "star" ] || [ "$ALIGNER" == "hisat2" ]; then
    
    bash ${ALIGNER}_alignment.sh $WORKING_DIR ${ALIGNMENT_INPUT} ${GENOME_URL} ${GTF_URL}
    multiqc -o $WORKING_DIR/QC_reports -n 4.multiqc_AlignmentStatus_${ALIGNER}.html "$WORKING_DIR/2_${ALIGNER}_output"

elif [ "$ALIGNER" == "kallisto" ] || [ "$ALIGNER" == "salmon" ]; then
    bash ${ALIGNER}_alignment.sh ${WORKING_DIR} ${ALIGNMENT_INPUT} ${TRANSCRIPTOME_URL}
    multiqc -o $WORKING_DIR/QC_reports -n 4.multiqc_AlignmentStatus_${ALIGNER}.html "$WORKING_DIR/3_${ALIGNER}_output"

fi


#---------- COUNTING ---------------

# Check if an argument was provided
if [ -n "${COUNTER}" ]; then
    # Argument is not null, run trimming script
    bash "${COUNTER}_count.sh" ${WORKING_DIR} ${ALIGNER} ${GTF_URL}

fi



