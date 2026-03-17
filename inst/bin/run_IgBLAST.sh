#!/bin/bash

# Depricated 2023/May/29

set -e

# Get command line arguments
while getopts f:r:o flag
do
    case "${flag}" in
        f) files=${OPTARG};;
        r) REF=${OPTARG};;
        o) OUT=${OPTARG};;
    esac
done

for file in ${files}
do
    filterbyname.sh \
    in=${file}/outs/multi/vdj_b/all_contig.fasta \
    out=Processing_QC/singlets_contigs/${file}_contigs.fasta \
    names=Processing_QC/${file}_singlets_barcodes.txt \
    substring=t \
    include=t \
    ow=t
    # Add $file (Pool) info to fasta
    sed -i "s/^>/>${file}_/g" Processing_QC/singlets_contigs/${file}_contigs.fasta
done

# If previous file exists, delete to overwrite with new data
if [ -e Processing_QC/singlets_contigs/${OUT}_singlets_contigs.fasta ]; then
    rm Processing_QC/singlets_contigs/${OUT}_singlets_contigs.fasta
fi

cat Processing_QC/singlets_contigs/*.fasta > Processing_QC/singlets_contigs/${OUT}_singlets_contigs.fasta

# Gene Assignment with IgBLAST-IMGT ref
AssignGenes.py igblast \
-s Processing_QC/singlets_contigs/${OUT}_singlets_contigs.fasta \
-b /fh/fast/mcelrath_j/McElrath_NGS/C_Marini_output/IMGT_ref/${REF}/igblast \
--organism human \
--loci ig \
--format airr
