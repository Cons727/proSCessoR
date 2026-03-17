#!/bin/bash
set -e

# ml MultiQC/1.9-foss-2019b-Python-3.7.4
# ml parallel
# ml FastQC

concat_fastqc(){

    # Make the output directory if it doesn't already exist
    mkdir -p FastQC_Reports

    echo "Merging R1 from ${2}${1}_L00*_R1_001.fastq.gz"
zcat ${2}${1}_L00*_R1_001.fastq.gz | fastqc stdin:${1}_R1 --outdir=FastQC_Reports &
    echo "$!" >> PID_LIST

    echo "Merging R2 from ${2}${1}_L00*_R2_001.fastq.gz"
zcat ${2}${1}_L00*_R2_001.fastq.gz | fastqc stdin:${1}_R2 --outdir=FastQC_Reports &
    echo "$!" >> PID_LIST
}

export -f concat_fastqc

find $1 -name "*.fastq.gz" | \
    while read F; do basename $F | \
    rev | \
    cut -c 22- | \
    rev; done | \
    sort | \
    uniq | \
    while read FASTQ_BASE; do

    echo "Running concat_fastqc for $FASTQ_BASE"
    concat_fastqc "${FASTQ_BASE}" "${1}"

done

# Wait for all processes to finish
cat PID_LIST | while read PID; do

    # While the PID is running, wait
    while ps -p $PID > /dev/null; do
        echo "Waiting for process $PID to finish"
        sleep 10
    done

done

multiqc FastQC_Reports


# Arguments:    1 = absolute path to parent directory of fastq files
# Example run: sbatch -c 16 --mail-user=cmarinim@fredhutch.org /fh/fast/mcelrath_j/McElrath_NGS/C_Marini_output/SOP/proSCessoR/concat_fastqc.sh /shared/ngs/illumina/mcelrath_j_SR_illumina/230412_A00613_0534_BHJ7JFDMXY/Unaligned/Project_mcelrath_j_SR_illumina/
# Reports will be saved in ./FastQC_Reports
