include { assembly as flye_assembly }                                        from './assembly'           addParams(assembler:'flye' )
include { assembly as miniasm_assembly }                                     from './assembly'           addParams(assembler: 'minasm' )
include { assembly as raven_assembly }                                       from './assembly'           addParams(assembler: 'raven')
include { assembly as unicycler_assembly }                                   from './assembly'           addParams(assembler: 'unicycler')
include { bgzip }                                                            from '../modules/bgzip'     addParams(params)
include { cluster; consensus; dotplot; msa; partition; reconcile; subsample} from '../modules/trycycler' addParams(params)

workflow trycycler {
    take:
    ch_fastq

    main:
    subsample(ch_fastq)

    subsample.out.fastq
        .multiMap { it ->
           flye:      tuple (it[0], it[0] + '_flye',      [it[1][1], it[1][5], it[1][9]])
           miniasm:   tuple (it[0], it[0] + '_miniasm',   [it[1][2], it[1][6], it[1][10]])
           raven:     tuple (it[0], it[0] + '_raven',     [it[1][3], it[1][7], it[1][11]])
           unicycler: tuple (it[0], it[0] + '_unicycler', [it[1][4], it[1][8], it[1][0]])
        }
        .set { ch_subsampled }

    flye_assembly(ch_subsampled.flye.transpose().map           { it -> tuple( it[1] + it[2].toString().replaceAll(~/.+sample/,"").replaceAll(~/.fastq/,""), it[2] )})
    miniasm_assembly(ch_subsampled.miniasm.transpose().map     { it -> tuple( it[1] + it[2].toString().replaceAll(~/.+sample/,"").replaceAll(~/.fastq/,""), it[2] )})
    raven_assembly(ch_subsampled.raven.transpose().map         { it -> tuple( it[1] + it[2].toString().replaceAll(~/.+sample/,"").replaceAll(~/.fastq/,""), it[2] )})
    unicycler_assembly(ch_subsampled.unicycler.transpose().map { it -> tuple( it[1] + it[2].toString().replaceAll(~/.+sample/,"").replaceAll(~/.fastq/,""), it[2] )})

    // ch_fastq
    //     .join(flye_assembly.out.fasta, by:1)
    //     .join(miniasm_assembly.out.fasta, by:1)
    //     .join(raven_assembly.out.fasta, by: 1)
    //     .join(unicycler_assembly.out.fasta, by:1)
    //     .view()
//        .groupTuple()
//        .view()
//        .set { ch_assemblies }

//     cluster(ch_assemblies)
//     dotplot(cluster.out.cluster)
//     reconcile(cluster.out.cluster)
//     msa(reconcile.out.cluster)
//     partition(msa.out.cluster)
//     consensus(partition.out.cluster)

//     emit:
//     consensus = consensus.out.fasta
    
}
