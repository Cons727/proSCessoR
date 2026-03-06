#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH --job-name CellRanger	# Job name
#SBATCH --out CellRanger.out	# name of the stdout
#SBATCH --mail-type=END
#SBATCH --mail-user=$1@fredhutch.org # send-to address

module load CellRanger/7.0.1

for file in ${2}/*.csv
do
    BATCH=`basename ${file} .csv | sed 's/CellRanger_Submission_//'`
    echo cellranger multi \
    --id=$BATCH \
    --csv=${file}
done

# Arguments:    1 = username
#               2 = absolut path to directory containing all CellRanger submission files
# Example run: sbatch -c 1 /fh/fast/mcelrath_j/McElrath_NGS/C_Marini_output/SOP/run_cellranger.sh cmarinim /fh/fast/mcelrath_j/McElrath_NGS/C_Marini_output/Submission_files

EOT
