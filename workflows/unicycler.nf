include { unicycler } from '../modules/unicycler.nf' addParams(unicycler_options: params.unicycler_options)

workflow unicycler {
    take:
    nanopore
    illumina

    main:
    unicycler(reads.join(illumina_fastqs, by: 0))
    
    emit:
    fasta = unicycler.out.fasta
}