include {raven} from '../modules/raven.nf' addParams(raven_options: params.raven_options)
include {flye} from '../modules/flye.nf' addParams(flye_options: params.flye_options)
include {miniasm} from '../modules/miniasm.nf' 
include {filtlong} from '../modules/filtlong.nf' addParams(filtlong_options: params.filtlong_options)
include {fastp} from '../modules/fastp.nf' addParams(fastp_options: params.fastp_options)
include {bgzip} from '../modules/bgzip.nf'

workflow donut_falls {
    take:
    fastq
    illumina

    main:
    fastp(illumina.join(fastq))
    filtlong(fastq.join(fastp.out.reads, by: 0, remainder: true))
    bgzip(filtlong.out.fastq)

    if (params.assembler == 'raven') {
        raven(bgzip.out.fastq)
        fastas = raven.out.fasta

    } else if (params.assembler == 'flye' ) {
        flye(bgzip.out.fastq)
        fastas = flye.out.fasta

    } else if (params.assembler == 'miniasm') {
        miniasm(bgzip.out.fastq)
        any2fasta(miniasm.out.gfa)
        fastas = any2fasta.out.fasta
    }

    emit:
    fastas = fastas
}