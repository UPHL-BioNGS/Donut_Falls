#!/usr/bin/env nextflow

params.any2fasta_options = ''

process any2fasta {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/any2fasta:latest'

  input:
  tuple val(sample), path(gfa) from miniasm_gfa

  output:
  tuple val(sample), path("${task.process}/${sample}.fasta") into assembled_fastas
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    any2fasta -v 2>> $err_file >> $log_file

    any2fasta !{params.any2fasta_options} \
      !{gfa} \
      2>> $err_file \
      > !{task.process}/!{sample}.fasta
  '''
}
