#!/usr/bin/env nextflow

params.nanoplot_options = '--barcoded'

process nanoplot {
  publishDir "${params.outdir}", mode: 'copy'
  tag "nanoplot"
  cpus params.medcpus
  container 'staphb/nanoplot:latest'

  when:
  params.nanoplot

  input:
  file(sequencing_summary) from sequencing_summary.view()

  output:
  path("${task.process}")
  path("logs/${task.process}/${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    NanoPlot --version | tee -a $log_file $err_file > /dev/null

    NanoPlot !{params.nanoplot_options} \
      --summary !{sequencing_summary} \
      --threads !{task.cpus} \
      --outdir !{task.process} \
      --raw \
      2>> $err_file >> $log_file
  '''
}
