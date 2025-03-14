#!/bin/bash

##############################################################

# Creates samplesheets for bwa_align_stats
# Usage: ./create_sample_sheet.sh

##############################################################

# Modify I/O as needed
FASTQ_DIR="/project/60006/mdprieto/raw_data/giardia/repositories/fastq"
HEADERS="sample,fastq_1,fastq_2,organism,reference_genome"
OUT_FILEPATH="samplesheet.csv"

#############################################################

# Create full samplesheet

echo $HEADERS  > "${OUT_FILEPATH}"

for fastq1 in $(ls ${FASTQ_DIR}/*_1.fastq*);
    do
        # get the basename of file and remove suffix
    sample_id=$(echo $fastq1 | xargs -n 1 basename | sed -E 's/.fq.*|.fastq.*//' | sed -E 's/_1$//')
    fastq2=${fastq1/_1/_2}
        # write in a new line for each sample
    row_data="$sample_id,$fastq1,$fastq2,Giardia_duodenalis,GiardiaDB-68_GintestinalisAssemblageEP15"
    echo $row_data >> "${OUT_FILEPATH}"
    done

