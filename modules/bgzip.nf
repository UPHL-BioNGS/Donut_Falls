process bgzip {
  tag "${sample}"
  cpus 1

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), path("filtlong/${fastq}.gz"),                       emit: fastq
  path "logs/bgzip/${sample}.${workflow.sessionId}.{log,err}", emit: logs

  shell:
  '''
    mkdir -p filtlong logs/bgzip
    log_file=logs/bgzip/!{sample}.!{workflow.sessionId}.log
    err_file=logs/bgzip/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    bgzip --version 2>> $err_file >> $log_file

    bgzip -@ !{task.cpus} !{fastq}
    mv !{fastq}.gz filtlong/.
  '''
}