include { raven }       from '../modules/raven.nf'      addParams(raven_options: params.raven_options)
include { flye }        from '../modules/flye.nf'       addParams(flye_options: params.flye_options)
include { miniasm }     from '../modules/miniasm.nf' 
include { filtlong }    from '../modules/filtlong.nf'   addParams(filtlong_options: params.filtlong_options)
include { fastp }       from '../modules/fastp.nf'      addParams(fastp_options: params.fastp_options)
include { bgzip }       from '../modules/bgzip.nf'
include { any2fasta }   from '../modules/any2fasta.nf'
include { medaka }      from '../modules/medaka.nf'     addParams(medaka_options: params.medaka_options)
include { polca }       from '../modules/masurca.nf'    addParams(polca_options: params.polca_options)
include { circlator }   from '../modules/circlator.nf'  

//TODO : add circlator somewhere, but there needs to be a check for circular sequences, and it should only be for assemblers that don't already rotate their sequences

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
        fasta = raven.out.fasta

    } else if (params.assembler == 'flye' ) {
        flye(bgzip.out.fastq)
        fasta = flye.out.fasta

    } else if (params.assembler == 'miniasm') {
        miniasm(bgzip.out.fastq)
        any2fasta(miniasm.out.gfa)
        fasta = any2fasta.out.fasta
    }

    medaka(fasta.join(fastq, by:0))
    polca(medaka.out.fasta.join(fastp.out.reads, by:0))

    emit:
    fasta = fasta
}