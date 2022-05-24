#!/usr/bin/env nextflow

params.flye_options = ''

process flye {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/flye:latest'
  //errorStrategy 'ignore'

  input:
  tuple val(sample), path(fastq) from filtered_fastq

  output:
  tuple val(sample), path("${task.process}/${sample}.fasta") into assembled_fastas
  path("${task.process}/${sample}/*")
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process} !{task.process}/!{sample}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    flye --version 2>> $err_file >> $log_file

    flye !{params.flye_options} \
      --nano-raw !{fastq} \
      --threads !{task.cpus} \
      --out-dir !{task.process}/!{sample} \
      2>> $err_file >> $log_file

    cp !{task.process}/!{sample}/assembly.fasta !{task.process}/!{sample}.fasta
  '''
}
