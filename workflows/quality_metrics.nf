include { nanoplot } from '../modules/nanoplot' addParams(params)
include { quast } from '../modules/quast' addParams(params)

workflow metrics {
    input:
    reads

    main:
    nanoplot(reads)
    quast(fasta)


    nanoplot.out.summary.collectFile(name: "NanoStats.csv",
    keepHeader: true,
    storeDir: "${params.outdir}/nanoplot")

    emit:
    summary = nanoplot.out.summary
}