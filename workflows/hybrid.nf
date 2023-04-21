include { unicycler } from '../modules/unicycler' addParams(params)
include { masurca }   from '../modules/masurca'   addParams(params)

workflow hybrid {
    take:
    ch_input

    main:
    ch_consensus = Channel.empty()
    ch_gfa       = Channel.empty()

    if (params.assembler == "unicycler" ) {
        unicycler(ch_input)
        ch_consensus = ch_consensus.mix(unicycler.out.fasta)
    } else if (params.assembler == "masurca") {
        masurca(ch_input)
        ch_consensus = ch_consensus.mix(masurca.out.fasta)
    }    
 
    emit:
    fasta = ch_consensus
}