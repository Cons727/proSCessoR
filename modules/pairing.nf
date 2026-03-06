
process split {
    container "${params.container__bbmap}"

    publishDir params.output_directory, mode: 'copy', overwrite: true

    input:
    // The working directory needs the files in Processing_QC
    // and CellRanger VDJ files so the input directory must be
    // the same as merge_pools
    path "*" // input_directory
    path "*" // Processing_QC from merge_pools

    output:
    path "*" // singlets_contigs

    script:
    """
    #!/bin/bash
    set -x
    
    # Set variables
    OUTNAME="${params.project_name}"
    
    BATCHES="\$(ls | grep "^P[0-9]")"

    echo Dirs: "\$BATCHES"
    echo Out name: "${params.project_name}"

    for batch in \${BATCHES}
    do
        echo "\$batch"
        if [ -d "\${batch}/outs/multi/${params.vdj_type}" ]; then
                in_file="\${batch}/outs/multi/${params.vdj_type}/all_contig.fasta"
        fi
        if [ -d "\${batch}/multi/${params.vdj_type}" ]; then
                in_file="\${batch}/multi/${params.vdj_type}/all_contig.fasta"
        fi

        filterbyname.sh \
        in="\${in_file}" \
        out=singlets_contigs/"\${batch}"_contigs.fasta \
        names=Processing_QC/"\${batch}"_singlets_barcodes.txt \
        substring=t \
        include=t \
        ow=t
        # Add file (Pool) info to fasta
        sed -i "s/^>/>"\${batch}"_/g" singlets_contigs/"\${batch}"_contigs.fasta
    done

    # If previous file exists, delete to overwrite with new data
    if [ -e singlets_contigs/"${params.project_name}"_singlets_contigs.fasta ]; then
        rm singlets_contigs/"${params.project_name}"_singlets_contigs.fasta
    fi

    cat singlets_contigs/*.fasta > singlets_contigs/"${params.project_name}"_singlets_contigs.fasta

    """
}

process annotate {
    container "${params.container__igblast}"

    publishDir params.output_directory, mode: 'copy', overwrite: true

    input:
    path "*" // singlets_contigs
    path vdj_path // VDJ reference

    output:
    path "*_singlets_contigs_igblast.tsv"

    script:
    def ref_vdj = vdj_path.getName()
    """
    #!/bin/bash
    set -e

    # Make sure to use the most recent IMGT or OGRDB reference version, or update to a new one using config_IMGT_ref.sh script.

    # Gene Assignment with IgBLAST
    # CORES=task.cpus
    export IGDATA="/usr/local/share/igblast/"

    # Set variables for OGRDB
    if [[ "${params.ref_vdj_type}" == "ogrdb" ]]; then
        DB_V=airr_c_human_ig.V
        DB_J=airr_c_human_ig.J
        DB_D=airr_c_human_igh.D
        DB_C=imgt_human_ig_c

        SEQTYPE=Ig
    fi

    # Set variables for IMGT
    if [[ "${params.ref_vdj_type}" == "imgt" ]]; then
        if [[ "${params.vdj_type}" == "vdj_b" ]]; then
            DB_V=igblast/database/imgt_human_ig_v
            DB_J=igblast/database/imgt_human_ig_j
            DB_D=igblast/database/imgt_human_ig_d
            DB_C=igblast/database/imgt_human_ig_c

            SEQTYPE=Ig
        fi
        if [[ "${params.vdj_type}" == "vdj_t" ]]; then
            DB_V=igblast/database/imgt_human_tr_v
            DB_J=igblast/database/imgt_human_tr_j
            DB_D=igblast/database/imgt_human_tr_d
            DB_C=igblast/database/imgt_human_tr_c

            SEQTYPE=TCR
        fi
    fi


    igblastn \
        -query singlets_contigs/"${params.project_name}"_singlets_contigs.fasta \
        -germline_db_V "${ref_vdj}"/"\${DB_V}" \
        -germline_db_D "${ref_vdj}"/"\${DB_D}" \
        -germline_db_J "${ref_vdj}"/"\${DB_J}" \
        -c_region_db "${ref_vdj}"/"\${DB_C}" \
        -auxiliary_data "\${IGDATA}"/optional_file/human_gl.aux \
        -organism human \
        -ig_seqtype "\${SEQTYPE}" \
        -domain_system imgt \
        -out "${params.project_name}"_singlets_contigs_igblast.tsv \
        -num_threads 8 \
        -outfmt 19

    """
}

process select {
    container "${params.container__R}"

    publishDir params.output_directory, mode: 'copy', overwrite: true

    input:
    path "*" // singlets_contigs_igblast.tsv
    path "*" // Processing_QC, merge_pools.out
    path "*" // input_directory
    path code_repository

    output:
    path "Processing_QC/*_singlets_prod_contig.*" // fasta and tsv

    script:
    """
    #!/bin/bash
    set -e

    # Move singlets_contigs_igblast.tsv to Processing_QC/ since R functions work there
    mkdir -p Processing_QC/singlets_contigs
    mv *_singlets_contigs_igblast.tsv Processing_QC/singlets_contigs/"${params.project_name}"_singlets_contigs_igblast.tsv
    
    # Add proSCessor to PATH and make scripts executable
    export PATH=\$PATH:${code_repository}
    chmod +x ${code_repository}/bin/*

    sbatch_contig_selection.R "${params.project_name}" "${params.vdj_type}"
    """
}

process grid {
    container "${params.container__igblast}"

    publishDir params.output_directory, mode: 'copy', overwrite: true

    input:
    path "*" // singlets_prod_contig fasta and tsv
    path vdj_path // VDJ reference

    output:
    path "*" // IgBLAST_SADIE_Runs

    script:
    def ref_vdj = vdj_path.getName()
    """
    #!/bin/bash
    set -e

    # FASTA file is always saved in Processing_QC/`project.name`_singlets_prod_contig.fasta in contig_selection.R
    export FASTA="${params.project_name}"_singlets_prod_contig.fasta
    export OUT="${params.project_name}"_IgBLAST_SADIE_Runs
    
    mkdir -p "\${OUT}"
    
    export IGDATA="/usr/local/share/igblast/"
    CORES=8

    # Set variables for OGRDB
    if [[ "${params.ref_vdj_type}" == "ogrdb" ]]; then
        DB_V=airr_c_human_ig.V
        DB_J=airr_c_human_ig.J
        DB_D=airr_c_human_igh.D
        DB_C=imgt_human_ig_c

        SEQTYPE=Ig
    fi

    # Set variables for IMGT
    if [[ "${params.ref_vdj_type}" == "imgt" ]]; then
        if [[ "${params.vdj_type}" == "vdj_b" ]]; then
            DB_V=igblast/database/imgt_human_ig_v
            DB_J=igblast/database/imgt_human_ig_j
            DB_D=igblast/database/imgt_human_ig_d
            DB_C=igblast/database/imgt_human_ig_c

            SEQTYPE=Ig
        fi
        if [[ "${params.vdj_type}" == "vdj_t" ]]; then
            DB_V=igblast/database/imgt_human_tr_v
            DB_J=igblast/database/imgt_human_tr_j
            DB_D=igblast/database/imgt_human_tr_d
            DB_C=igblast/database/imgt_human_tr_c

            SEQTYPE=TCR
        fi
    fi


    echo "Annotating "\${FASTA}" with IgBLAST"
    date
    for V_penalty in 1 2 3
    do
    for J_penalty in 1 2 3
    do
        echo V_penalty=-"\${V_penalty}", D_penalty=-2, J_penalty=-"\${J_penalty}"
        time igblastn \
        -query "\${FASTA}" \
        -germline_db_V "${ref_vdj}"/"\${DB_V}" \
        -germline_db_D "${ref_vdj}"/"\${DB_D}" \
        -germline_db_J "${ref_vdj}"/"\${DB_J}" \
        -c_region_db "${ref_vdj}"/"\${DB_C}" \
        -auxiliary_data "\${IGDATA}"/optional_file/human_gl.aux \
        -V_penalty -"\${V_penalty}" \
        -D_penalty -2 \
        -J_penalty -"\${J_penalty}" \
        -extend_align5end \
        -extend_align3end \
        -organism human \
        -ig_seqtype "\${SEQTYPE}" \
        -domain_system imgt \
        -show_translation \
        -out "\${OUT}"/Run_V"\${V_penalty}"_D2_J"\${J_penalty}".tsv \
        -num_threads "\${CORES}" \
        -evalue 0.001 \
        -outfmt 19
    done
    done

    echo "V_penalty=-1, D_penalty=-4, J_penalty=-3, allow_vdj_overlap"
    time igblastn \
    -query "\${FASTA}" \
    -germline_db_V "${ref_vdj}"/"\${DB_V}" \
    -germline_db_D "${ref_vdj}"/"\${DB_D}" \
    -germline_db_J "${ref_vdj}"/"\${DB_J}" \
    -c_region_db "${ref_vdj}"/"\${DB_C}" \
    -auxiliary_data "\${IGDATA}"/optional_file/human_gl.aux \
    -V_penalty -1 \
    -D_penalty -4 \
    -J_penalty -3 \
    -extend_align5end \
    -extend_align3end \
    -allow_vdj_overlap \
    -organism human \
    -ig_seqtype "\${SEQTYPE}" \
    -domain_system imgt \
    -show_translation \
    -out "\${OUT}"/Run_V1_D4_J3_overlap.tsv \
    -num_threads "\${CORES}" \
    -evalue 0.001 \
    -outfmt 19        
    """
}

process pair {
    container "${params.container__R}"

    publishDir params.output_directory, mode: 'copy', overwrite: true

    input:
    // The working directory needs the files in Processing_QC
    // and CellRanger VDJ files so the input directory must be
    // the same as merge_pools
    path "*" // IgBLAST_SADIE_Runs
    path "*" // Processing_QC, merge_pools.out
    path "*" // singlets_prod_contig fasta and tsv, select.out
    path "*" // input_directory
    path code_repository

    output:
    path "*" // ig.sadie_paired_airr.tsv and ig.sadie_wide.tsv

    script:
    """
    #!/bin/bash
    set -e

    # Move singlets_contigs_igblast.tsv to Processing_QC/ since R functions work there
    mkdir -p Processing_QC/${params.project_name}_IgBLAST_SADIE_Runs
    cp ${params.project_name}_IgBLAST_SADIE_Runs/* Processing_QC/${params.project_name}_IgBLAST_SADIE_Runs/
    mv *_singlets_prod_contig.tsv Processing_QC/"${params.project_name}"_singlets_prod_contig.tsv

    # Add proSCessor to PATH and make scripts executable
    export PATH=\$PATH:${code_repository}
    chmod +x ${code_repository}/bin/*

    sbatch_pairing.R "${params.project_name}" "${params.vdj_type}"

    # Move output files outside Processing_QC to channel them to clone process
    cp Processing_QC/${params.project_name}_ig.sadie_paired_airr.tsv ${params.project_name}_ig.sadie_paired_airr.tsv
    cp Processing_QC/${params.project_name}_sadie_wide.tsv ${params.project_name}_sadie_wide.tsv
    """
}

process clone {
    container "${params.container__R}"

    publishDir params.output_directory, mode: 'copy', overwrite: true

    input:
    // The working directory needs the files in Processing_QC
    // and CellRanger VDJ files so the input directory must be
    // the same as merge_pools
    path "*" // pair.out, ig.sadie_paired_airr.tsv and ig.sadie_wide.tsv
    path "*" // Processing_QC, merge_pools.out
    path code_repository

    output:
    path "*" // Final airr output

    script:
    """
    #!/bin/bash
    set -e

    # Move files to Processing_QC/ since R functions work there
    cp ${params.project_name}_sadie_wide.tsv Processing_QC/${params.project_name}_sadie_wide.tsv

    # Add proSCessor to PATH and make scripts executable
    export PATH=\$PATH:${code_repository}
    chmod +x ${code_repository}/bin/*

    if [ "${params.clone}" == "true" ]; then
        changeo-10x.sh \
        -a "${params.project_name}"_ig.sadie_paired_airr.tsv \
        -o Processing_QC/"${params.project_name}"_changeo \
        -n "${params.project_name}" \
        -t ig \
        -x auto \
        -e gmm \
        -z

        merge_Changeo_SADIE.R "${params.project_name}"
    else
        merge_md_SADIE.R "${params.project_name}"
    fi
    """
}

workflow pairing {
    take:
        // Channel containing the output files produced by merge_pools
        merge_pools_out
        input_directory
        code_repository
        ref_vdj

    main:

        // Run 
       split(
             input_directory,
             merge_pools_out
        )
        annotate(
            split.out,
            ref_vdj
        )
        select(
            annotate.out,
            merge_pools_out,
            input_directory,
            code_repository
        )
        grid(
            select.out,
            ref_vdj
        )
        pair(
            grid.out,
            merge_pools_out,
            input_directory,
            select.out,
            code_repository
        )
        clone(
            pair.out,
            merge_pools_out,
            code_repository
        )

    emit:
        clone.out

}