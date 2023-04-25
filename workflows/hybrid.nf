include { bandage }   from '../modules/bandage'   addParams(params)
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
        ch_gfa = ch_gfa.mix(unicycler.out.gfa)
        ch_consensus = ch_consensus.mix(unicycler.out.fasta)
    } else if (params.assembler == "masurca") {
        masurca(ch_input)
        ch_gfa = ch_gfa.mix(masurca.out.gfa)
        ch_consensus = ch_consensus.mix(masurca.out.fasta)
    }    
 
    bandage(ch_gfa)

    emit:
    fasta   = ch_consensus
    summary = bandage.out.summary
}