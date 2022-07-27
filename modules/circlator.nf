process circlator {
  tag "${sample}"
  cpus 1

  input:
  tuple val(sample), file(fasta)

  output:
  tuple val(sample), file("circlator/${sample}/${sample}.fasta"),   emit: fasta
  path("circlator/${sample}/*"),                                    emit: directory
  path("logs/circlator/${sample}.${workflow.sessionId}.{log,err}"), emit: logs

  shell:
  '''
    mkdir -p circlator/!{sample} logs/circlator
    log_file=logs/circlator/!{sample}.!{workflow.sessionId}.log
    err_file=logs/circlator/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    circlator --version 2>> $err_file >> $log_file

    circlator fixstart  
        !{fasta} \
        circlator/!{sample} \
        2>> $err_file >> $log_file

    exit 1
  '''
}