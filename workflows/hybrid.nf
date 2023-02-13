include { unicycler } from '../modules/unicycler' addParams(params)
include { masurca }   from '../modules/masurca'   addParams(params)

workflow hybrid {
    take:
    reads

    main:
    if (params.assembler == "unicyler" ) {
        unicycler(reads)
        consensus = unicycler.out.fasta
    } else if (params.assembler == "masurca") {
        masurca(reads)
        consensus = masurca.out.fasta
    }    
 
    emit:
    fasta = consensus
}