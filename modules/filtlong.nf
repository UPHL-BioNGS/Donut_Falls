process filtlong {
  publishDir "donut_falls", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/filtlong:latest'

  input:
  tuple val(sample), file(fastq), file(short_reads)

  output:
  tuple val(sample), file("filtlong/${sample}_filtered.fastq"), optional: true, emit: fastq
  path "logs/filtlong/${sample}.${workflow.sessionId}.{log,err}",               emit: logs

  shell:
  if (short_reads[1] == null) {
  '''
    mkdir -p filtlong logs/filtlong
    log_file=logs/filtlong/!{sample}.!{workflow.sessionId}.log
    err_file=logs/filtlong/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    filtlong --version >> $log_file

    filtlong !{params.filtlong_options} \
      !{fastq} 2>> $err_file > filtlong/!{sample}_filtered.fastq
  '''
  } else {
  '''
    mkdir -p filtlong logs/filtlong
    log_file=logs/filtlong/!{sample}.!{workflow.sessionId}.log
    err_file=logs/filtlong/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    filtlong --version >> $log_file

    filtlong !{params.filtlong_options} \
      -1 !{short_reads[0]} \
      -2 !{short_reads[1]} \
      !{fastq} \
      2>> $err_file > filtlong/!{sample}_filtered.fastq
    '''
  }
}