include { bwa }    from '../modules/bwa'     addParams(params)
include { medaka } from '../modules/medaka'  addParams(params)
include { polca }  from '../modules/masurca' addParams(params)
include { polypolish } from '../modules/polypolish' addParams(params)
//include { pilon }  from '../modules/pilon' addParams(params)

workflow polish {
    take:
    ch_fastq
    ch_fasta
    ch_illumina

    main:
    medaka(ch_fasta.join(ch_fastq, by:0))
    bwa(medaka.out.fasta.join(ch_illumina, by:0))
    polypolish(bwa.out.sam)
    polca(polypolish.out.fasta.join(ch_illumina, by:0))

    // mostly so that everything goes through busco
    ch_medaka_polished = medaka.out.fasta.map{it     -> tuple(it[0] + "_medaka"    , it[1])}
    ch_poly_polished   = polypolish.out.fasta.map{it -> tuple(it[0] + "_polypolish", it[1])}
    ch_polca_polished  = polca.out.fasta.map{it      -> tuple(it[0] + "_polca"     , it[1])}
    
    ch_fasta
        .mix(ch_medaka_polished)
        .mix(ch_poly_polished)
        .mix(ch_polca_polished)
        .set{ch_consensus}

    // add a process to compress fasta files for download from nf-tower?

    emit:
    fasta = ch_consensus

}