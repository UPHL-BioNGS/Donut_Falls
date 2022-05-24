#!/usr/bin/env nextflow

process bgzip {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/ivar:latest'

  input:
  tuple val(sample), path(fastq) from filtlong_fastq

  output:
  tuple val(sample), path("filtlong/${fastq}.gz") optional true into filtered_fastq, filtered_fastq_medaka
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p filtlong logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    bgzip --version 2>> $err_file >> $log_file

    bgzip -@ {task.cpus} !{fastq}
    mv !{fastq}.gz filtlong/.
  '''
}
