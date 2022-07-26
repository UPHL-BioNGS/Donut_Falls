process medaka {
  publishDir "donut_falls", mode: 'copy'
  tag "${sample}"
  cpus 6
  container 'ontresearch/medaka:latest'

  input:
  tuple val(sample), path(fasta), path(fastq)

  output:
  path "medaka/${sample}/",                                                     emit: directory
  tuple val(sample), path("medaka/${sample}/${sample}_medaka_consensus.fasta"), emit: fasta
  path "logs/medaka/${sample}.${workflow.sessionId}.{log,err}",                 emit: logs

  shell:
  '''
    mkdir -p logs/medaka medaka
    log_file=logs/medaka/!{sample}.!{workflow.sessionId}.log
    err_file=logs/medaka/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    medaka --version >> $log_file

    medaka_consensus !{params.medaka_options} \
      -i !{fastq} \
      -d !{fasta} \
      -o medaka/!{sample} \
      -t 2

    cp medaka/!{sample}/consensus.fasta medaka/!{sample}/!{sample}_medaka_consensus.fasta
  '''
}