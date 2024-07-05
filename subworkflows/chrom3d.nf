// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT


workflow CHROM3D {

    take:
        beads
        args
        number_of_models

    main:

        GENERATE_SEEDS(
            beads,
            number_of_models
        )

        GENERATE_SEEDS.out.txt
            .splitText()
            .map {
                toks = it[1].trim().split("\t")
                tuple(it[0], toks[0], toks[1])
            }
            .set { seeds }

        beads.combine(seeds, by: 0)
            .set { chrom3d_tasks }

        SIMULATE(
            chrom3d_tasks,
            args
        )

        ARCHIVE(
            SIMULATE.out.cmm.groupTuple()
        )

    emit:
        tar = ARCHIVE.out.tar
}


process GENERATE_SEEDS {
    label 'process_very_short'

    tag "${sample}"

    input:
        tuple val(sample),
              path(files)
        val num_seeds

    output:
        tuple val(sample),
              path(outname), emit: txt

    shell:
        outname="${sample}.seeds.txt"
        files_str=files.join(" ")
        '''
        generate_seed_sequence.py !{files_str} --number-of-seeds='!{num_seeds}' > '!{outname}'
        '''
}

process SIMULATE {
    tag "${sample}_${id}"

    label 'error_retry'
    label 'duration_long'

    input:
        tuple val(sample),
              path(beads),
              val(id),
              val(seed)
        val args

    output:
        tuple val(sample),
              path("*.cmm"),
        emit: cmm

    shell:
        outname="${beads.simpleName}_${id}.cmm"
        args_str=args.join(" ")
        '''
        Chrom3D -o '!{outname}' \\
                -s '!{seed}' \\
                --nucleus '!{beads}' \\
                !{args}
        '''
}

process ARCHIVE {
    publishDir "${params.publish_dir}/models",
        enabled: !!params.publish_dir,
        mode: params.publish_dir_mode

    tag "${sample}"

    input:
        tuple val(sample),
              path(models)

    output:
        tuple val(sample),
              path(outname),
        emit: tar

    shell:
        outname="${sample}.tar.gz"
        '''
        tar --transform 's,^,!{sample}/,' -czf '!{outname}' *.cmm
        '''
}
