#!/bin/bash

#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=IMGT_ref
#SBATCH --output=IMGT_ref.out

module load IgBLAST/1.22.0-x64-linux
module load Python/3.8.2-GCCcore-9.3.0

PATH=$PATH:/home/cmarinim/.local/bin:/home/cmarinim/proSCessoR/Get_IMGTref

DATE=$(date +"%Y.%m.%d")
IMGTDIR=/fh/fast/mcelrath_j/McElrath_NGS/C_Marini_output/IMGT_ref/ref-${DATE}/germlines/imgt
IGBLASTDIR=/fh/fast/mcelrath_j/McElrath_NGS/C_Marini_output/IMGT_ref/ref-${DATE}/igblast

# Build IgBLAST database from IMGT reference sequences
fetch_imgtdb.sh -o $IMGTDIR
imgt2igblast.sh -i $IMGTDIR -o $IGBLASTDIR


# Example
# sbatch -c 8 --mail-user=$(whoami)@fredhutch.org config_IMGT_ref.sh