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

    emit:
    fasta = polca.out.fasta
}