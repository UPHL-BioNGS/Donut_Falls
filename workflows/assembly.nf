include { circlator }                   from '../modules/circlator' addParams(params)
include { flye }                        from '../modules/flye'      addParams(params)
include { gfastats }                    from '../modules/gfastats'  addParams(params)
include { medaka }                      from '../modules/medaka'    addParams(params)
include { miniasm }                     from '../modules/miniasm'   addParams(params)
include { raven }                       from '../modules/raven'     addParams(params)
include { unicycler_long as unicycler } from '../modules/unicycler' addParams(params)

workflow assembly {
    take:
    ch_fastq

    main:
    ch_gfa   = Channel.empty()
    ch_fasta = Channel.empty()

    if (params.assembler == 'raven' ) {
        raven(ch_fastq)
        ch_gfa = ch_gfa.mix(raven.out.gfa)
    } else if (params.assembler == 'flye' ) {
        flye(ch_fastq)
        ch_gfa = ch_gfa.mix(flye.out.gfa)
    } else if (params.assembler == 'miniasm' ) {
        miniasm(ch_fastq)
        ch_gfa = ch_gfa.mix(miniasm.out.gfa)
    } else if (params.assembler == 'lr_unicycler' ) {
        unicycler(ch_fastq)
        ch_fasta = unicycler.out.fasta
    }

    gfastats(ch_gfa)
    circlator(gfastats.out.fasta)
    medaka(circlator.out.fasta.mix(ch_fasta).join(ch_fastq, by:0))

    emit:
    fasta = medaka.out.fasta
}