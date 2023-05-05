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
        ch_gfa       = ch_gfa.mix(unicycler.out.gfa)
        ch_consensus = ch_consensus.mix(unicycler.out.fasta)
    } else if (params.assembler == "masurca") {
    //    masurca(ch_input)
    //    ch_consensus = ch_consensus.mix(masurca.out.fasta)
        println("We really wanted to support including of masurca, but it became too time consuming.")
        println("If this assembler is useful to you, please submit an issue at https://github.com/UPHL-BioNGS/Donut_Falls/issues")
        ch_consensus = Channel.empty()
    }    
 
    bandage(ch_gfa)

    emit:
    fasta   = ch_consensus
    summary = bandage.out.summary
}