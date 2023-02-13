include { medaka } from '../modules/medaka'     addParams(params)
include { polca } from '../modules/masurca'    addParams(params)
include { pilon } from '../modules/pilon' addParams(params)

workflow polish {
    take:
    fastq
    illumina

    main:

    polca(medaka.out.fasta.join(fastp.out.reads, by:0))

    emit:
    fasta = fasta
}