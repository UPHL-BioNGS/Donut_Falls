#!/usr/bin/env nextflow

params.fastp_options = ''

process fastp {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'bromberglab/fastp:latest'

  when:
  sample != null

  input:
  tuple val(sample), path(fastq) from illumina_fastqs

  output:
  tuple val(sample), path("${task.process}/${sample}_{R1,R2}.fastq.gz") into clean_reads
  path("${task.process}")
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastp --version >> $log_file 2>> $err_file

    fastp !{params.fastp_options} \
      --in1 !{fastq[0]} \
      --in2 !{fastq[1]} \
      --out1 !{task.process}/!{sample}_R1.fastq.gz \
      --out2 !{task.process}/!{sample}_R2.fastq.gz \
      --unpaired1 !{task.process}/!{sample}_u.fastq.gz \
      --unpaired2 !{task.process}/!{sample}_u.fastq.gz \
      2>> $err_file >> $log_file

    cp *.html !{task.process}/.
    cp *.json !{task.process}/.
  '''
}
