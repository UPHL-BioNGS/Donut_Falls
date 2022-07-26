process any2fasta {
  publishDir "donut_falls", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/any2fasta:latest'

  input:
  tuple val(sample), file(gfa)

  output:
  tuple val(sample), file("miniasm/${sample}.fasta"),               emit: fasta
  path("logs/any2fasta/${sample}.${workflow.sessionId}.{log,err}"), emit: logs

  shell:
  '''
    mkdir -p miniasm/!{sample} logs/any2fasta
    log_file=logs/any2fasta/!{sample}.!{workflow.sessionId}.log
    err_file=logs/any2fasta/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    any2fasta -v 2>> $err_file >> $log_file

    any2fasta \
      !{gfa} \
      2>> $err_file \
      > miniasm/!{sample}.fasta
  '''
}