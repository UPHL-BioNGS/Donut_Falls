#!/usr/bin/env nextflow

params.raven_options = '--polishing-rounds 2'
process raven {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/raven:latest'

  input:
  tuple val(sample), path(fastq) from filtered_fastq

  output:
  tuple val(sample), path("${task.process}/${sample}/${sample}.fasta") into raven_fastas
  path("${task.process}/${sample}/*")
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    raven --version 2>> $err_file >> $log_file

    raven !{params.raven_options} \
      --threads !{task.cpus} \
      --graphical-fragment-assembly !{task.process}/!{sample}/!{sample}.gfa \
      !{fastq} \
      2>> $err_file \
      > !{task.process}/!{sample}/!{sample}.fasta
  '''
}
