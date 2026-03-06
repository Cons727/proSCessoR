#!/bin/bash

#SBATCH --mail-type=FAIL,END

# REF (full path to reference folder) and DB_NAME (imgt or airr_c) from sbatch export
# REF=//fh/fast/_VIDD/HVTN/Vaccine/B_cell_Sequencing/OGRDB_ref/ref-2024.01.11
# DB=airr_c
# sbatch run_IgSADIE.sh --mail-user=cmarinim@fredhutch.org

set -e

# FASTA and OUT from sbatch export
# FASTA name of singlets
# OUT Output directory name: project.name_IgBLAST_SADIE_Runs

GERMLINE_DB=$REF
CORES=$(nproc)

mkdir -p ${OUT}

module load IgBLAST/1.22.0-x64-linux
export IGDATA=$EBROOTIGBLAST

# Are these variables exported from sbatch_pair_clone.sh? Remove this chunk if yes
# Set variables for OGRDB
if [[ "${DB_NAME}" == "ogrdb" ]]; then
    DB_V=airr_c_human_ig.V
    DB_J=airr_c_human_ig.J
    DB_D=airr_c_human_igh.D
    DB_C=imgt_human_ig_c

    SEQTYPE=Ig
fi

# Set variables for IMGT
if [[ "${DB_NAME}" == "imgt" ]]; then
    if [[ "${VDJ}" == "vdj_b" ]]; then
        DB_V=igblast/database/imgt_human_ig_v
        DB_J=igblast/database/imgt_human_ig_j
        DB_D=igblast/database/imgt_human_ig_d
        DB_C=igblast/database/imgt_human_ig_c

        SEQTYPE=Ig
    fi
    if [[ "${VDJ}" == "vdj_t" ]]; then
        DB_V=igblast/database/imgt_human_tr_v
        DB_J=igblast/database/imgt_human_tr_j
        DB_D=igblast/database/imgt_human_tr_d
        DB_C=igblast/database/imgt_human_tr_c

        SEQTYPE=TCR
    fi
fi


echo "Annotating $FASTA with IgBLAST"
date

# SADIE penalties
# https://github.com/jwillis0720/sadie/blob/967cd1ddfbd6dfc5c67220839c52067211f596ae/src/sadie/airr/airr.py#L555
for v_gene_penalty in 1 2 3
do
  for j_gene_penalty in 1 2
  do
    echo "v_gene_penalty=-${v_gene_penalty}, d_gene_penalty=-1, j_gene_penalty=-${j_gene_penalty}"
    time igblastn \
      -query ${FASTA} \
      -germline_db_V ${GERMLINE_DB}/${DB_V} \
      -germline_db_D ${GERMLINE_DB}/${DB_D} \
      -germline_db_J ${GERMLINE_DB}/${DB_J} \
      -c_region_db ${GERMLINE_DB}/${DB_C} \
      -auxiliary_data ${IGDATA}/optional_file/human_gl.aux \
      -V_penalty -"${v_gene_penalty}" \
      -D_penalty -1 \
      -J_penalty -"${j_gene_penalty}" \
      -extend_align5end \
      -extend_align3end \
      -organism human \
      -ig_seqtype $SEQTYPE \
      -domain_system imgt \
      -show_translation \
      -out ${OUT}/SADIE_V${v_gene_penalty}_D1_J${j_gene_penalty}.tsv \
      -num_threads ${CORES} \
      -evalue 0.001 \
      -outfmt 19
  done
done

echo "v_gene_penalty=-1, d_gene_penalty=-4, j_gene_penalty=-3, allow_vdj_overlap"
time igblastn \
  -query ${FASTA} \
  -germline_db_V ${GERMLINE_DB}/${DB_V} \
  -germline_db_D ${GERMLINE_DB}/${DB_D} \
  -germline_db_J ${GERMLINE_DB}/${DB_J} \
  -c_region_db ${GERMLINE_DB}/${DB_C} \
  -auxiliary_data ${IGDATA}/optional_file/human_gl.aux \
  -V_penalty -1 \
  -D_penalty -4 \
  -J_penalty -3 \
  -extend_align5end \
  -extend_align3end \
  -allow_vdj_overlap \
  -organism human \
  -ig_seqtype $SEQTYPE \
  -domain_system imgt \
  -show_translation \
  -out ${OUT}/SADIE_V1_D4_J3_overlap.tsv \
  -num_threads ${CORES} \
  -evalue 0.001 \
  -outfmt 19

