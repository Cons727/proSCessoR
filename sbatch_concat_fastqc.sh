#!/bin/bash

#SBATCH --job-name fastqc	# Job name
#SBATCH --out fastqc.out	# name of the stdout
#SBATCH --mail-type=END,FAIL

ml MultiQC/1.9-foss-2019b-Python-3.7.4
ml parallel
ml FastQC

mkdir -p FastQC_Reports

concat_fastqc(){
    echo "Merging R1"
zcat $2"$1"_L00*_R1_001.fastq.gz | fastqc stdin:"$1"_R1 --outdir=FastQC_Reports

echo "Merging R2"
zcat $2"$1"_L00*_R2_001.fastq.gz | fastqc stdin:"$1"_R2 --outdir=FastQC_Reports
}

export -f concat_fastqc

find $1 -type f -name "*.fastq.gz" | while read F; do basename $F | rev | cut -c 22- | rev; done | sort | uniq | parallel -j8 -a - concat_fastqc ::: $1

multiqc FastQC_Reports --outdir FastQC_Reports


# Arguments:    1 = absolut path to parent directory of fastq files
# Example run: sbatch -c 16 --mail-user=cmarinim@fredhutch.org /fh/fast/mcelrath_j/McElrath_NGS/C_Marini_output/SOP/proSCessoR/concat_fastqc.sh /shared/ngs/illumina/mcelrath_j_SR_illumina/230412_A00613_0534_BHJ7JFDMXY/Unaligned/Project_mcelrath_j_SR_illumina/
# Reports will be saved in ./FastQC_Reports
