#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: ${workflow.manifest.version}")
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
params.sample_sheet               = 'sample_sheet.csv'
params.assembler                  = 'flye'
params.outdir                     = 'donut_falls'

params.circlator_options          = ''
params.fastp_options              = ''
params.filtlong_options           = '--min_length 1000 --keep_percent 95'
params.flye_options               = ''
params.gfastats_options           = ''
params.masurca_options            = ''
params.medaka_options             = ''
params.nanoplot_summary_options   = ''
params.nanoplot_options           = ''
params.polca_options              = ''
params.porechop_options           = ''
params.raven_options              = '--polishing-rounds 2'
params.trycycler_subsample_options = ''
params.unicycler_options          = ''

include { assembly }                     from './workflows/assembly'  addParams(params)
include { filter }                       from './workflows/filter'    addParams(params)
include { hybrid }                       from './workflows/hybrid'    addParams(params)
include { nanoplot_summary as nanoplot } from './modules/nanoplot'    addParams(params)
include { trycycler }                    from './workflows/trycycler' addParams(params)


Channel
  .fromPath(params.sequencing_summary, type:'file')
  .view { "Summary File : $it" }
  .set { ch_sequencing_summary }

Channel
  .fromPath("${params.sample_sheet}", type: "file")
  .splitCsv( header: true, sep: ',' )
  .map { row -> tuple( "${row.sample}", file("${row.fastq}"), [file("${row.fastq_1}"), file("${row.fastq_2}") ]) }
  .set { ch_input_files }

workflow {
  if ( params.assembler == 'flye' || params.assembler == 'raven' || params.assembler == 'miniasm' || params.assembler == 'lr_unicycler' ) {
    filter(ch_input_files)
    assembly(filter.out.fastq)
  } else if ( params.assembler == 'unicycler' || params.assembler == 'masurca' ) {
    hybrid(ch_input_files)
  } else if ( params.assembler == 'trycycler' ) {
    filter(ch_input_files)
    trycycler(filter.out.fastq)
  }

  nanoplot(ch_sequencing_summary)

}
