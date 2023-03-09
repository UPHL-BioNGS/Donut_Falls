include { download } from './modules/download' addParams(params)

workflow download {
    input:
    ch_sra_accessions

    main:
    download(ch_sra_accessions)

    emit:
    fastq = download.out.fastq
}