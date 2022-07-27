#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
version = "0.0.20220810"

println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version:" + version)
println("")

params.config_file                = false
if (params.config_file) {
  def src = new File("${workflow.projectDir}/configs/donut_falls_config_template.config")
  def dst = new File("${workflow.launchDir}/edit_me.config")
  dst << src.text
  println("A config file can be found at ${workflow.launchDir}/edit_me.config")
  exit 0
}


//# TODO : add circlator
//# TODO : add something to evaluate fastas
//# TODO : add summary file
//# TODO : add trycycler
//# TODO : add multiqc
//# TODO : metric subworkflow

params.sequencing_summary         = workflow.launchDir + "/*sequencing_summary*txt"
params.reads                      = workflow.launchDir + '/reads'
params.illumina                   = workflow.launchDir + '/illumina'
params.assembler                  = 'flye'
params.outdir                     = 'donut_falls'

params.nanoplot_summary_options   = '--barcoded'
params.nanoplot_options           = ''
params.nanoplot_illumina_options  = ''
params.unicycler_options          = ''
params.raven_options              = '--polishing-rounds 2'
params.flye_options               = ''
params.filtlong_options           = '--min_length 1000 --keep_percent 95'
params.medaka_options             = ''
params.polca_options              = ''
params.fastp_options              = ''

include { donut_falls }                                   from './workflows/donut_falls.nf' addParams(assembler: params.assembler,
                                                                                            raven_options: params.raven_options,
                                                                                            flye_options: params.flye_options,
                                                                                            polca_options: params.polca_options,
                                                                                            medaka_options: params.medaka_options,
                                                                                            filtlong_options: params.filtlong_options,
                                                                                            fastp_options: params.fastp_options)
include { unicycler }                                     from './workflows/unicycler.nf'   addParams(unicycler_options: params.unicycler_options)
include { nanoplot; nanoplot_summary; nanoplot_illumina } from './modules/nanoplot.nf'      addParams(nanoplot_options: params.nanoplot_options,
                                                                                            nanoplot_summary_options: params.nanoplot_summary_options,
                                                                                            nanoplot_illumina_options: params.nanoplot_illumina_options)

Channel
  .fromPath(params.sequencing_summary, type:'file')
  .view { "Summary File : $it" }
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
  nanoplot(fastq)
  nanoplot_illumina(illumina_fastq)
  nanoplot_summary(sequencing_summary)
  if ( params.assembler == 'flye' || params.assembler == 'raven' || params.assembler == 'miniasm' ) {
    donut_falls(fastq, illumina_fastq)
  } else if ( params.assembler == 'unicycler' ) {
    unicycler(fastq, illumina_fastq)
  }

  nanoplot.out.summary.collectFile(name: "NanoStats.csv",
    keepHeader: true,
    storeDir: "${params.outdir}/nanoplot")
}
