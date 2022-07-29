include { unicycler } from '../modules/unicycler.nf' addParams(unicycler_options: params.unicycler_options)

workflow hybrid {
    take:
    reads

    main:
    reads
    unicycler(reads)
    
    emit:
    fasta = unicycler.out.fasta
}