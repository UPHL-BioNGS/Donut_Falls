#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
version = "0.0.20220810"

println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version:" + version)
println("")

params.sequencing_summary = workflow.launchDir + "/*sequencing_summary*txt"
params.reads = workflow.launchDir + '/reads'
params.illumina = workflow.launchDir + '/illumina'
params.assembler = 'flye'

params.nanoplot_summary_options = ''
params.nanoplot_options = ''
params.unicycler_options = ''
params.raven_options = '--polishing-rounds 2'
params.flye_options = ''
params.filtlong_options = '--min_length 1000 --keep_percent 95'
params.medaka_options = ''
params.polca_options = ''
params.fastp_options = ''

include { donut_falls } from './workflows/donut_falls.nf' addParams(assembler: params.assembler,
                                                                    raven_options: params.raven_options,
                                                                    flye_options: params.flye_options,
                                                                    polca_options: params.polca_options,
                                                                    filtlong_options: params.filtlong_options,
                                                                    fastp_options: params.fastp_options)
include { unicycler } from './workflows/unicycler.nf' addParams(unicycler_options: params.unicycler_options)

Channel
  .fromPath(params.sequencing_summary, type:'file')
  .view { "Summary File : $it" }
  .ifEmpty{
    println("Could not find sequencing summary file! Set with 'params.sequencing_summary'")
  }
  .set { sequencing_summary }

Channel
  .fromPath("${params.reads}/*.{fastq,fastq.gz,fq,fq.gz}", type:'file')
  .ifEmpty {
    println("Could not find fastq files for nanopore sequencing. Set with 'params.reads'")
    exit 1
  }
  .map { reads -> tuple(reads.simpleName, reads ) }
  .view { "Fastq file found : ${it[0]}" }
  .set { fastq }

Channel
  .fromFilePairs("${params.illumina}/*_R{1,2}*.{fastq,fastq.gz}", size: 2 )
  .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .view { "Illumina fastq files for for greater accuracy : ${it[0]}" }
  .set { illumina_fastq }

workflow {
  if ( params.assembler == 'flye' || params.assembler == 'raven' || params.assembler == 'miniasm' ) {
    donut_falls(fastq, illumina_fastq)
  } else if ( params.assembler == 'unicycler' ) {
    unicycler(fastq, illumina_fastq)
  }
}
