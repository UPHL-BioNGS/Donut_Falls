include { nanoplot } from './modules/nanoplot.nf' addParams(nanoplot_options = params.nanoplot_options)

workflow metrics {
    input:
    reads

    main:
    nanoplot(reads)

    emit:
    summary = nanoplot.out.summary
}