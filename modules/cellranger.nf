process concat_fastqc {
    container "${params.container__fastqc}"
    publishDir params.output_directory, mode: 'copy', overwrite: true

    input:
    path input_directory

    output:
    path "FastQC_Reports/*"

    """#!/bin/bash
concat_fastqc.sh "${input_directory}/"
    """
}

process cellranger_multi {
    container "${params.container__cellranger}"
    publishDir "${params.output_directory}/cellranger_multi/", mode: 'copy', overwrite: true

    input:
    path cr_submission_dir

    output:
    path "*"

    """#!/bin/bash
    set -e
    
    FILES=\$(find "$cr_submission_dir" -name "*.csv")

    for file in $FILES
    do
        # Keep Pool from *_CellRanger_Submission.csv file name
        POOL=\$(echo $file | sed 's/_CellRanger_Submission.*//')

        echo "CellRanger Configuration File: \$file"
        cat "$file"

        echo "Processing batch \$POOL"
        cellranger multi \
            --id=\$POOL \
            --csv="$file" \
            --localcores=${task.cpus} \
            --localmem=${task.memory.toGiga()}
    done
    """

}

workflow cellranger {

    take:
        cr_submission_dir

    main:
    
        // Run the concat_fastqc.sh script on the input files
        concat_fastqc(input_directory)

        // Run CellRanger
        cellranger_multi(
            cr_submission_dir
        )

    emit:
        cellranger_multi.out

}