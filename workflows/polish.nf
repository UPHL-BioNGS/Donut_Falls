include { medaka } from '../modules/medaka'  addParams(params)
include { polca }  from '../modules/masurca' addParams(params)
//include { pilon }  from '../modules/pilon' addParams(params)

workflow polish {
    take:
    ch_fastq
    ch_fasta
    ch_illumina

    main:
    medaka(ch_fasta.join(ch_fastq, by:0))
    polca(medaka.out.fasta.join(ch_illumina, by:0))

    // mostly so that everything goes through busco
    ch_medaka_polished   = medaka.out.fasta.map{it -> tuple(it[0] + "_medaka"   , it[1])}
    ch_illumina_polished = polca.out.fasta.map{it  -> tuple(it[0] + "_polca"    , it[1])}
    ch_fasta.map{it -> tuple(it[0] + "_assembled", it[1])}
        .mix(ch_medaka_polished)
        .mix(ch_illumina_polished)
        .set{ch_consensus}

    emit:
    fasta = ch_consensus

}