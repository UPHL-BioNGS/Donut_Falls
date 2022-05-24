#!/usr/bin/env nextflow

params.filtlong_options = "--min_length 1000 --keep_percent 95"

process filtlong {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.${workflow.sessionId}.{log,err}"
  tag "${sample}"
  cpus 1
  container 'staphb/filtlong:latest'

  input:
  tuple val(sample), file(fastq), file(short_reads) from fastq.join(clean_reads, by:0, remainder : true)

  output:
  tuple val(sample), path("${task.process}/${sample}_filtered.fastq") optional true into filtlong_fastq
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  if (short_reads[1] == null) {
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    filtlong --version >> $log_file

    filtlong !{params.filtlong_options} \
      !{fastq} 2>> $err_file > !{task.process}/!{sample}_filtered.fastq
  '''
  } else {
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    filtlong --version >> $log_file

    filtlong !{params.filtlong_options} \
      -1 !{short_reads[0]} \
      -2 !{short_reads[1]} \
      !{fastq} \
      2>> $err_file > !{task.process}/!{sample}_filtered.fastq
    '''
  }
}
