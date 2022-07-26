  process miniasm {
  publishDir "donut_falls", mode: 'copy'
  tag "${sample}"
  cpus 12
  container 'staphb/minipolish:latest'

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), path("miniasm/${sample}/*gfa"),             emit: gfa
  path "miniasm/${sample}/*",                                    emit: directory
  path "logs/miniasm/${sample}.${workflow.sessionId}.{log,err}", emit: logs

  shell:
  '''
    mkdir -p miniasm/!{sample} logs/miniasm
    log_file=logs/miniasm/!{sample}.!{workflow.sessionId}.log
    err_file=logs/miniasm/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "miniasm version : $(miniasm -V)" 2>> $err_file >> $log_file
    minimap2 --version 2>> $err_file >> $log_file
    minipolish --version 2>> $err_file >> $log_file

    miniasm_and_minipolish.sh !{params.miniasm_options} \
      !{fastq} \
      !{task.cpus} \
      2>> $err_file > miniasm/!{sample}/!{sample}.gfa
  '''
}

