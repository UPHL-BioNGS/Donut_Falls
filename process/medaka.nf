#!/usr/bin/env nextflow

params.medaka_options = ''

process medaka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 6
  container 'ontresearch/medaka:latest'

  input:
  tuple val(sample), path(fasta), path(fastq) from assembled_fastas.join(filtered_fastq_medaka, by:0 )

  output:
  path("${task.process}/${sample}/")
  tuple val(sample), path("${task.process}/${sample}/${sample}_medaka_consensus.fasta") into medaka_fastas
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process} !{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    medaka --version >> $log_file

    medaka_consensus !{params.medaka_options} \
      -i !{fastq} \
      -d !{fasta} \
      -o !{task.process}/!{sample} \
      -t 2

    cp !{task.process}/!{sample}/consensus.fasta !{task.process}/!{sample}/!{sample}_medaka_consensus.fasta
  '''
}
