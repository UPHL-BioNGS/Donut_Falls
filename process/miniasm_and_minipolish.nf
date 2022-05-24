#!/usr/bin/env nextflow

params.miniasm_and_minipolish_options = ''
process miniasm_and_minipolish {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/minipolish:latest'

  input:
  tuple val(sample), path(fastq) from filtered_fastq

  output:
  tuple val(sample), path("${task.process}/${sample}/*gfa") into miniasm_gfa
  path("${task.process}/${sample}/*")
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "miniasm version : $(miniasm -V)" 2>> $err_file >> $log_file
    minimap2 --version 2>> $err_file >> $log_file
    minipolish --version 2>> $err_file >> $log_file

    miniasm_and_minipolish.sh \
      !{fastq} \
      !{task.cpus} \
      2>> $err_file > !{task.process}/!{sample}/!{sample}.gfa
  '''
}
