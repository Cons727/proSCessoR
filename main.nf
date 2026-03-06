#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

include { merge_pools } from './modules/merge_pools'
include { pairing } from './modules/pairing'

workflow {

    // Get the input files provided by the user-defined parameters

    // CODE REPOSITORY - Not defined by user
    Channel
        .fromPath(
            "$baseDir",
            checkIfExists: true,
            type: 'dir'
        )
        .set { code_repository }

    // INPUT DIRECTORY - Cell Ranger output (parent dir)
    Channel
        .fromPath(
            params.input_directory,
            checkIfExists: true,
            type: "dir"
        )
        .toSortedList()
        .set { input_directory }

    // METADATA FILE
    Channel
        .fromPath(
            params.metadata,
            checkIfExists: true,
            type: "file"
        )
        .set { metadata }
    
    // VDJ REF - Included in package
    Channel
        .fromPath(
            "${projectDir}/data/${params.ref_vdj}",
            checkIfExists: true,
            type: "dir"
        )
        .set { ref_vdj }

    // Merge pools, background correction and demultiplexing
    merge_pools(
         input_directory,
         metadata,
         code_repository
    )

    // merge_pools.out.view { "Processed files: $it" }

    // BCR QC and pairing
    pairing(
        merge_pools.out,
        input_directory,
        code_repository,
        ref_vdj
    )
}
