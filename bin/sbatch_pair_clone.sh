#!/bin/bash

# Example run:
# sbatch -c 16 --mail-user=cmarinim@fredhutch.org -o prSCR-40.3.log -J prSCR-40.3 --mail-type=END,FAIL \
# --export=OUTNAME=HVTN301_BCS2040,REF=/fh/fast/_VIDD/HVTN/Vaccine/B_cell_Sequencing/OGRDB_ref/ref-2024.01.11/airr_c_human,\
# DB_NAME=ogrdb,VDJ=vdj_b,FASTA=Processing_QC/HVTN301_BCS2040_singlets_prod_contig.fasta sbatch_pair_clone.sh

set -e


# You can have any version of python listed in your PATH,
# just make sure Change-O is installed in the same python version
# that it is loaded when you load R. For example,
# ml fhR/4.1.1-foss-2020b
# loads python3.8, while
# ml fhR/4.3.1-foss-2022b
# loads python3.10
export PATH=$PATH:/home/$(whoami)/.local/lib/python3.8/site-packages:/home/$(whoami)/.local/bin:/home/$(whoami)/proSCessoR/bin
export FILES=`ls | grep "^P[0-9]"`

echo Files: $FILES
echo Reference: $REF
echo Out name: $OUTNAME
echo VDJ type: $VDJ

# You can load any version of python,
# just make sure is the same version that you exported to your PATH
ml BBMap Python/3.8.2-GCCcore-9.3.0 IgBLAST/1.22.0-x64-linux 

for file in ${FILES}
do
    # Added '|| true' to prevent script from exiting due to 'set -e' when directories are missing
    in_file=$(find ${file}/outs/multi/${VDJ} ${file}/multi/${VDJ} ${file}/outs/${VDJ} -name "all_contig.fasta" -print -quit 2>/dev/null || true)

    # Ensure we actually found a file before running the tool
    if [[ -z "$in_file" ]]; then
        echo "WARNING: Could not find all_contig.fasta for ${file}. Skipping."
        continue
    fi

    filterbyname.sh \
    in=${in_file} \
    out=Processing_QC/singlets_contigs/${file}_contigs.fasta \
    names=Processing_QC/${file}_singlets_barcodes.txt \
    substring=t \
    include=t \
    ow=t
    # Add $file (Pool) info to fasta
    sed -i "s/^>/>${file}_/g" Processing_QC/singlets_contigs/${file}_contigs.fasta
done

# If previous file exists, delete to overwrite with new data
if [ -e Processing_QC/singlets_contigs/${OUTNAME}_singlets_contigs.fasta ]; then
    rm Processing_QC/singlets_contigs/${OUTNAME}_singlets_contigs.fasta
fi

cat Processing_QC/singlets_contigs/*.fasta > Processing_QC/singlets_contigs/${OUTNAME}_singlets_contigs.fasta

echo "VDJ Annotation with IgBLAST"

# Gene Assignment with IgBLAST
CORES="$(nproc)"
export IGDATA=$EBROOTIGBLAST

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

igblastn \
      -query Processing_QC/singlets_contigs/${OUTNAME}_singlets_contigs.fasta \
      -germline_db_V ${REF}/${DB_V} \
      -germline_db_D ${REF}/${DB_D} \
      -germline_db_J ${REF}/${DB_J} \
      -c_region_db ${REF}/${DB_C} \
      -auxiliary_data ${IGDATA}/optional_file/human_gl.aux \
      -organism human \
      -ig_seqtype $SEQTYPE \
      -domain_system imgt \
      -out Processing_QC/singlets_contigs/${OUTNAME}_singlets_contigs_igblast.tsv \
      -num_threads ${CORES} \
      -outfmt 19

# You can load any version of R (proSCessoR and its dependencies should be installed in this version),
# just make sure Change-O is installed in the same python version
# that it is loaded when you load R. For example,
# ml fhR/4.1.1-foss-2020b
# loads python3.8, while
# ml fhR/4.3.1-foss-2022b
# loads python3.10
# Suggestion: use the same python version as in your PATH above
ml fhR/4.1.1-foss-2020b

echo "Contig Selection"
sbatch_contig_selection.R $OUTNAME $VDJ $IGM_IGD

echo "Running SADIE"
export OUT=Processing_QC/${OUTNAME}_IgBLAST_SADIE_Runs
run_IgSADIE.sh 

echo "Annotation and pairing"
sbatch_pairing.R $OUTNAME $VDJ

changeo-10x.sh -a Processing_QC/${OUTNAME}_ig.sadie_paired_airr.tsv -o Processing_QC/${OUTNAME}_changeo -n $OUTNAME -x auto -e gmm -z

merge_Changeo_SADIE.R $OUTNAME