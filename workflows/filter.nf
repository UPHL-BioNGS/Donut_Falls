include { bgzip }    from '../modules/bgzip'    addParams(params)
include { fastp }    from '../modules/fastp'    addParams(params)
include { filtlong } from '../modules/filtlong' addParams(params)
include { porechop } from '../modules/porechop' addParams(params)

workflow filter {
    take:
    ch_input

    main:
    fastp(ch_input.map { it -> tuple (it[0], it[2])})
    porechop(ch_input.map {it -> tuple (it[0], it[1])})
    filtlong(porechop.out.fastq.join(fastp.out.reads, by: 0, remainder: true))
    bgzip(filtlong.out.fastq)

    emit:
    fastq = bgzip.out.fastq
    reads = fastp.out.reads
}