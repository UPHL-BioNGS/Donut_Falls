include { download      } from '../modules/download' addParams(params)
include { great_dataset } from '../modules/download' addParams(params)

workflow test {
    input:

    main:
    if ( params.assembler == 'trycycler' ) { 
        great_dataset()
        fastq = great_dataset.out.fastq
    } else { 
        download()
        fastq = download.out.fastq
    }

    emit:
    fastq = fastq
}