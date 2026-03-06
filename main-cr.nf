#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

include { cellranger } from './modules/cellranger'
include { merge_pools } from './modules/merge_pools'
include { pairing } from './modules/pairing'

workflow {

    // Get the input files provided by the user-defined parameters

    // CODE REPOSITORY
    Channel
        .fromPath(
            "$baseDir",
            checkIfExists: true,
            type: 'dir'
        )
        .set { code_repository }

    // INPUT DIRECTORY
    Channel
        .fromPath(
            params.input_directory,
            checkIfExists: true,
            type: "dir"
        )
        .set { input_directory }

    // CR SUBMISSION DIRECTORY
    Channel
        .fromPath(
            params.cr_submission_dir,
            checkIfExists: true,
            type: "dir"
        )
        .set { cr_submission_dir }

    // METADATA FILE
    Channel
        .fromPath(
            params.metadata,
            checkIfExists: true,
            type: "dir"
        )
        .set { metadata }
    
    // VDJ REF
    Channel
        .fromPath(
            params.ref_vdj,
            checkIfExists: true,
            type: "dir"
        )
        .set { ref_vdj }

    // Run QC on input reads and align using CellRanger
    cellranger(
        cr_submission_dir
    )

    // Merge pools
    merge_pools(
         cellranger.out,
         metadata,
         code_repository,
    )

    // Background correction and demultiplexing
    pairing(
        merge_pools.out
        metadata
        ref_vdj
        code_repository
    )
}
