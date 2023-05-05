include { assembly as flye_assembly }                                         from './assembly'           addParams(assembler: 'flye' )
include { assembly as miniasm_assembly }                                      from './assembly'           addParams(assembler: 'miniasm' )
include { assembly as raven_assembly }                                        from './assembly'           addParams(assembler: 'raven')
include { assembly as unicycler_assembly }                                    from './assembly'           addParams(assembler: 'lr_unicycler')
include { cluster; consensus; dotplot; msa; partition; reconcile; subsample } from '../modules/trycycler' addParams(params)
include { combine }                                                           from '../modules/trycycler' addParams(params)
include { rasusa }                                                            from '../modules/rasusa'    addParams(params)

workflow trycycler {
    take:
    ch_fastq
    ch_remove

    main:
    subsample(ch_fastq)
    rasusa(subsample.out.full)

    subsample.out.fastq
        .mix(rasusa.out.fastq)
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

    flye_assembly.out.assembly
        .mix(miniasm_assembly.out.assembly)
        .mix(raven_assembly.out.assembly)
        .mix(unicycler_assembly.out.assembly)
        .map { it -> tuple( it[0].replaceAll(~/_flye.+/,"").replaceAll(~/_miniasm.+/,"").replaceAll(~/_raven.+/,"").replaceAll(~/_unicycler.+/,""), it[1])}
        .groupTuple()
        .join(ch_fastq, by:0)
        .set { ch_assemblies }

    cluster(ch_assemblies)
    //dotplot(cluster.out.cluster.transpose())
    reconcile(cluster.out.cluster.join(ch_fastq, by: 0).transpose().combine(ch_remove))
    msa(reconcile.out.cluster)
    partition(msa.out.cluster.groupTuple().join(ch_fastq, by: 0).transpose())
    consensus(partition.out.cluster)
    combine(consensus.out.fasta.groupTuple())

    flye_assembly.out.summary
        .mix(raven_assembly.out.summary)
        .mix(unicycler_assembly.out.summary)
        .mix(miniasm_assembly.out.summary)
        .branch { it ->
            gfastats:  it =~ /gfastats/
            circlator: it =~ /circlator/
            other: true
        }
        .set { ch_for_summary }

    // ch_for_summary.gfastats
    //     .collectFile(
    //         storeDir: "${params.outdir}/gfastats/",
    //         keepHeader: true,
    //         name: "gfastats_summary.csv")
    //     .set { gfastats_summary }

    // ch_for_summary.circlator
    //     .collectFile(
    //         storeDir: "${params.outdir}/circlator/",
    //         keepHeader: true,
    //         name: "circlator_summary.csv")
    //     .set { circlator_summary }

    emit:
    fasta   = combine.out.fasta
    summary = ch_for_summary.other
}
