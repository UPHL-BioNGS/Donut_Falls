#!/usr/bin/env nextflow

println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.0.20220316\n")

// TODO : Add longQC?
// TODO : Add summary file
// TODO : Add socru?


// Import processes
include {nanoplot} from './process/nanoplot.nf'
include {fastp} from './process/fastp.nf'
include {filtlong} from './process/filtlong.nf'
include {bgzip} from './process/bgzip.nf'
include {flye} from './process/flye.nf'
include {miniasm_and_minipolish} from './process/miniasm_and_minipolish'
include {any2fasta} from './process/any2fasta'
include {raven} from './process/raven'
include {medaka} from './process/medaka'
include {polca} from './process/polca'

params.outdir = workflow.launchDir + '/donut_falls'
println("The files and directory for results is " + params.outdir)

params.maxcpus = Runtime.runtime.availableProcessors()
println("Maximum number of CPUS used in this workflow : ${params.maxcpus}")
if ( params.maxcpus < 12 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 12
}

params.sequencing_summary = "${workflow.launchDir}/*sequencing_summary*txt"
Channel
  .fromPath(params.sequencing_summary, type:'file')
  .view { "Summary File : $it" }
  .ifEmpty{
    params.nanoplot = false
    println("Could not find sequencing summary file! Set with 'params.sequencing_summary'")
  }
  .set { sequencing_summary }

if ( params.nanoplot != false ) { params.nanoplot = true }

params.reads = workflow.launchDir + '/reads'
Channel
  .fromPath("${params.reads}/*.{fastq,fastq.gz,fq,fq.gz}", type:'file')
  .ifEmpty {
    println("Could not find fastq files for nanopore sequencing. Set with 'params.reads'")
    exit 1
  }
  .map { reads -> tuple(reads.simpleName, reads ) }
  .view { "Fastq file found : ${it[0]}" }
  .set { fastq }

params.illumina = workflow.launchDir + '/illumina'
Channel
  .fromFilePairs("${params.illumina}/*_R{1,2}*.{fastq,fastq.gz}", size: 2 )
  .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .view { "Illumina fastq files for for greater accuracy : ${it[0]}" }
  .into { illumina_fastqs ; illumina_fastqs_polishing }

params.assembler = 'flye'
//params.assembler = 'raven'
//params.assembler = 'miniasm'


// Workflows

workflow nanoplot {
  main:
  sequencing_summary
  nanoplot(sequencing_summary)
}
