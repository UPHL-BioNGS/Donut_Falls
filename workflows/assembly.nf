include { bandage }                     from '../modules/bandage'       addParams(params)
include { circlator }                   from '../modules/circlator'     addParams(params)
include { dragonflye }                  from '../modules/dragonflye'    addParams(params)
include { flye }                        from '../modules/flye'          addParams(params)
include { gfastats }                    from '../modules/gfastats'      addParams(params)
include { medaka }                      from '../modules/medaka'        addParams(params)
include { miniasm }                     from '../modules/miniasm'       addParams(params)
include { raven }                       from '../modules/raven'         addParams(params)
include { unicycler_long as unicycler } from '../modules/unicycler'     addParams(params)

workflow assembly {
    take:
    ch_fastq

    main:
    ch_gfa       = Channel.empty()
    ch_summary   = Channel.empty()

    if (params.assembler == 'raven' ) {
        raven(ch_fastq)
        ch_gfa = ch_gfa.mix(raven.out.gfa)
    } else if (params.assembler == 'flye' ) {
        flye(ch_fastq)
        ch_gfa = ch_gfa.mix(flye.out.gfa)

        flye.out.summary
            .collectFile(
                storeDir: "${params.outdir}/flye/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "flye_summary.tsv")
            .set { flye_summary }

        ch_summary=ch_summary.mix(flye_summary)
    } else if (params.assembler == 'miniasm' ) {
        miniasm(ch_fastq)
        ch_gfa = ch_gfa.mix(miniasm.out.gfa)
    } else if (params.assembler == 'lr_unicycler' ) {
        unicycler(ch_fastq)
        ch_gfa = ch_gfa.mix(unicycler.out.gfa)
    } else if (params.assembler == 'dragonflye' ) {
        dragonflye(ch_fastq)
        ch_gfa = dragonflye.out.gfa

        dragonflye.out.summary
            .collectFile(
                storeDir: "${params.outdir}/dragonflye/",
                keepHeader: true,
                sort: { file -> file.text },
                name: "dragonflye_summary.tsv")
            .set { dragonflye_summary }

        ch_summary=ch_summary.mix(dragonflye_summary)
    }

    bandage(ch_gfa)
    gfastats(ch_gfa)
    circlator(gfastats.out.fasta)

    gfastats.out.summary
        .collectFile(
            storeDir: "${params.outdir}/gfastats/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "gfastats_summary.csv")
        .set { gfastats_summary }

    circlator.out.summary
        .collectFile(
            storeDir: "${params.outdir}/circlator/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "circlator_summary.csv")
        .set { circlator_summary }

    emit:
    fasta       = circlator.out.fasta
    assembly    = gfastats.out.assembly
    summary     = ch_summary.mix(bandage.out.summary).mix(gfastats_summary).mix(circlator_summary)
}