#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: ${workflow.manifest.version}")
println("")

params.config_file                 = false
if (params.config_file) {
  def src = new File("${workflow.projectDir}/configs/donut_falls_config_template.config")
  def dst = new File("${workflow.launchDir}/edit_me.config")
  dst << src.text
  println("A config file can be found at ${workflow.launchDir}/edit_me.config")
  exit 0
}

params.sequencing_summary          = workflow.launchDir + "/*sequencing_summary*txt"
params.sample_sheet                = ''
params.assembler                   = 'flye'
params.outdir                      = 'donut_falls'
params.remove                      = 'remove.txt'
params.reads                       = ''

params.busco_options               = ''
params.circlator_options           = ''
params.enable_porechop             = false
params.fastp_options               = ''
params.filtlong_options            = '--min_length 1000 --keep_percent 95'
params.flye_options                = ''
params.gfastats_options            = ''
params.masurca_options             = ''
params.medaka_options              = ''
params.multiqc_options             = ''
params.nanoplot_summary_options    = ''
params.nanoplot_options            = ''
params.polca_options               = ''
params.polypolish_options          = ''
params.porechop_options            = ''
params.quast_options               = ''
params.rasusa_options              = '--frac 80'
params.raven_options               = '--polishing-rounds 2'
params.trycycler_subsample_options = ''
params.trycycler_cluster_options   = ''
params.trycycler_consensus_options = ''
params.trycycler_dotplot_options   = ''
params.trycycler_msa_options       = ''
params.trycycler_partition_options = ''
params.trycycler_reconcile_options = ''
params.unicycler_options           = ''

include { assembly }                     from './workflows/assembly'  addParams(params)
include { filter }                       from './workflows/filter'    addParams(params)
include { hybrid }                       from './workflows/hybrid'    addParams(params)
include { nanoplot_summary as nanoplot } from './modules/nanoplot'    addParams(params)
include { polish }                       from './workflows/polish'    addParams(params)
include { trycycler }                    from './workflows/trycycler' addParams(params)
include { metrics }                      from './workflows/metrics'   addParams(params)

Channel
  .fromPath(params.sequencing_summary, type:'file')
  .view { "Summary File : $it" }
  .set { ch_sequencing_summary }

if ( params.sample_sheet) {
  Channel
    .fromPath("${params.sample_sheet}", type: "file")
    .splitCsv( header: true, sep: ',' )
    .map { it -> tuple( "${it.sample}", "${it.fastq}", "${it.fastq_1}", "${it.fastq_2}" ) }
    .branch { it ->
       sr:  it[2] != it[3]
       other: true
    }
    .set{ch_precheck}

    ch_precheck.sr.map { it -> tuple(it[0], file(it[1]), [file(it[2]), file(it[3])])}
      .mix(ch_precheck.other.map { it -> tuple(it[0], file(it[1]), null)})
      .set{ch_input_files}

} else if ( params.reads ) {
  Channel
    .fromPath("${params.reads}/*.{fastq,fastq.gz,fq,fq.gz}", type:'file')
    .ifEmpty {
      println("Could not find fastq files for nanopore sequencing. Set with 'params.reads'")
      exit 1
    }
    .map { reads -> tuple(reads.simpleName, reads, null ) }
    .view { "Fastq file found : ${it[0]}" }
    .set { ch_input_files }
} else {
  println("Thank you for using Donut Falls!\n")
  println("Please set a sample sheet with params.sample_sheet or a directory with params.reads!")
  exit 0
}

ch_remove = Channel.fromPath("${params.remove}", type: "file")

workflow {
  ch_illumina    = Channel.empty()
  ch_fastq       = Channel.empty()
  ch_fasta       = Channel.empty()
  ch_consensus   = Channel.empty()
  ch_summary     = Channel.empty()

  if ( params.assembler == 'flye' || params.assembler == 'raven' || params.assembler == 'miniasm' || params.assembler == 'lr_unicycler' ) {
    filter(ch_input_files)
    ch_fastq     = ch_fasta.mix(filter.out.fastq)
    ch_illumina  = ch_illumina.mix(filter.out.reads)

    assembly(filter.out.fastq)
    ch_fasta     = ch_fasta.mix(assembly.out.fasta)
    ch_summary   = ch_summary.mix(filter.out.summary)

  } else if ( params.assembler == 'unicycler' || params.assembler == 'masurca' ) {
    hybrid(ch_input_files.filter{it -> it[2]})
    ch_consensus = ch_consensus.mix(hybrid.out.fasta)

  } else if ( params.assembler == 'trycycler' ) {
    filter(ch_input_files)
    ch_fastq     = ch_fasta.mix(filter.out.fastq)
    ch_illumina  = ch_illumina.mix(filter.out.reads)
    ch_summary   = ch_summary.mix(filter.out.summary)

    trycycler(filter.out.fastq, ch_remove.ifEmpty([]))
    ch_fasta    = ch_fasta.mix(trycycler.out.fasta)
  }

  nanoplot(ch_sequencing_summary)
  polish(ch_fastq, ch_fasta, ch_illumina.ifEmpty([]))
  ch_consensus = ch_consensus.mix(polish.out.fasta)

  metrics(
    ch_input_files.map{it -> tuple (it[0], it[1])},
    ch_consensus,
    ch_summary.ifEmpty([]))
}
