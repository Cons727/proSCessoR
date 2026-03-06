#!/bin/bash
sbatch <<EOT
#!/bin/bash

#SBATCH --job-name changeo	# Job name
#SBATCH --out ChangeO.out	# name of the stdout
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user="$1"@fredhutch.org # send-to address

module load Singularity


singularity exec -B "$2":/data \
/fh/fast/mcelrath_j/McElrath_NGS/C_Marini_output/immcantation_suite-4.4.0.sif changeo-10x \
-s /data/post_e-Vq_filtered_contig.fasta \
-a /data/post_e-Vq_filtered_contig_annotations.csv \
-o /data/"$3"/ \
-n $3 \
-g human \
-t ig \
-x auto \
-e gmm

# Arguments:    1 = username
#               2 = absolut path to parent directory of contig.fasta and contig_annotations.csv
#               3 = name of file
# Example run: sbatch -c 1 /fh/fast/mcelrath_j/McElrath_NGS/C_Marini_output/SOP/run_changeo.sh cmarinim /fh/fast/mcelrath_j/McElrath_NGS/C_Marini_output/HVTN137-scCITE-seq/Downsampling/ Testing_sh

EOT