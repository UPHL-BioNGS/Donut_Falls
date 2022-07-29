process nanoplot_summary {
  tag "${sequencing_summary}"
  cpus 6

  input:
  file(sequencing_summary)

  output:
  path "nanoplot/summary",                                      emit: final_directory                                  
  path "logs/nanoplot_summary/${workflow.sessionId}.{log,err}", emit: logs

  shell:
  '''
    mkdir -p nanoplot/summary logs/nanoplot_summary
    log_file=logs/nanoplot_summary/!{workflow.sessionId}.log
    err_file=logs/nanoplot_summary/!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    NanoPlot --version | tee -a $log_file $err_file > /dev/null

    NanoPlot !{params.nanoplot_summary_options} \
      --summary !{sequencing_summary} \
      --threads !{task.cpus} \
      --outdir nanoplot/summary \
      --tsv_stats \
      2>> $err_file >> $log_file
  '''
}

process nanoplot {
  tag "${sample}"
  cpus 6

  input:
  tuple val(sample), file(fastq)

  output:
  path "nanoplot/${sample}",                                      emit: final_directory
  path "nanoplot/${sample}/${sample}_NanoStats.csv",              emit: summary
  path "logs/nanoplot/${sample}.${workflow.sessionId}.{log,err}", emit: logs

  shell:
  '''
    mkdir -p nanoplot/!{sample} logs/nanoplot
    log_file=logs/nanoplot/!{sample}.!{workflow.sessionId}.log
    err_file=logs/nanoplot/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    NanoPlot --version | tee -a $log_file $err_file > /dev/null

    NanoPlot !{params.nanoplot_options} \
      --fastq !{fastq} \
      --threads !{task.cpus} \
      --tsv_stats \
      --outdir nanoplot/!{sample} \
      2>> $err_file >> $log_file

    cp nanoplot/!{sample}/NanoStats.txt nanoplot/!{sample}/!{sample}_NanoStats.txt

    echo "sample,$(cut -f 1 nanoplot/!{sample}/!{sample}_NanoStats.txt | tr '\\n' ',' )" >  nanoplot/!{sample}/!{sample}_NanoStats.csv
    echo "!{sample},$(cut -f 2 nanoplot/!{sample}/!{sample}_NanoStats.txt | tr '\\n' ',' )" >> nanoplot/!{sample}/!{sample}_NanoStats.csv
  '''
}

process nanoplot_illumina {
  tag "${sample}"
  cpus 6

  input:
  tuple val(sample), file(fastq)

  output:
  path "nanoplot/${sample}_il",                                            emit: final_directory
  tuple val(sample), file("nanoplot/${sample}_il/NanoStats.txt"),          emit: summary
  path "logs/nanoplot_illumina/${sample}.${workflow.sessionId}.{log,err}", emit: logs

  shell:
  '''
    mkdir -p nanoplot/!{sample}_il logs/nanoplot_illumina
    log_file=logs/nanoplot_illumina/!{sample}.!{workflow.sessionId}.log
    err_file=logs/nanoplot_illumina/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    NanoPlot --version | tee -a $log_file $err_file > /dev/null

    NanoPlot !{params.nanoplot_options} \
      --fastq !{fastq} \
      --threads !{task.cpus} \
      --tsv_stats \
      --outdir nanoplot/!{sample}_il \
      2>> $err_file >> $log_file

  '''
}