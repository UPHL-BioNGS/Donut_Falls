process nanoplot_summary {
  publishDir "donut_falls", mode: 'copy'
  tag "${sequencing_summary}"
  cpus 6
  container 'staphb/nanoplot:latest'

  input:
  file(sequencing_summary)

  output:
  path "nanoplot_summary",                                      emit: final_directory                                  
  path "logs/nanoplot_summary/${workflow.sessionId}.{log,err}", emit: logs

  shell:
  '''
    mkdir -p nanoplot_summary/basic logs/nanoplot_summary
    log_file=logs/nanoplot_summary/!{workflow.sessionId}.log
    err_file=logs/nanoplot_summary/!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    NanoPlot --version | tee -a $log_file $err_file > /dev/null

    NanoPlot \
      --summary !{sequencing_summary} \
      --threads !{task.cpus} \
      --outdir nanoplot_summary \
      --tsv_stats \
      2>> $err_file >> $log_file
  '''
}

process nanoplot {
  publishDir "donut_falls", mode: 'copy'
  tag "${sample}"
  cpus 6
  container 'staphb/nanoplot:latest'

  input:
  tuple val(sample), file(fastq)

  output:
  path "nanoplot/${sample}",                                      emit: final_directory
  tuple val(sample), file("nanoplot/summary.tsv"),                emit: summary
  path "logs/nanoplot/${sample}.${workflow.sessionId}.{log,err}", emit: logs

  shell:
  '''
    mkdir -p nanoplot/!{sample} logs/nanoplot
    log_file=logs/nanoplot/!{sample}.!{workflow.sessionId}.log
    err_file=logs/nanoplot/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    NanoPlot --version | tee -a $log_file $err_file > /dev/null

    NanoPlot \
      --fastq !{fastq} \
      --threads !{task.cpus} \
      --tsv_stats \
      --outdir nanoplot/!{sample} \
      2>> $err_file >> $log_file
  '''
}