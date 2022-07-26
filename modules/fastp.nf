process fastp {
  publishDir "donut_falls", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/fastp:latest'

  input:
  tuple val(sample), file(fastq), file(nanopore)

  output:
  tuple val(sample), file("fastp/${sample}_{R1,R2}.fastq.gz"), emit: reads
  path "fastp/*"                                             , emit: directory
  path "logs/fastp/${sample}.${workflow.sessionId}.{log,err}", emit: logs

  shell:
  '''
    mkdir -p fastp logs/fastp
    log_file=logs/fastp/!{sample}.!{workflow.sessionId}.log
    err_file=logs/fastp/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastp --version >> $log_file 2>> $err_file

    fastp !{params.fastp_options} \
      --in1 !{fastq[0]} \
      --in2 !{fastq[1]} \
      --out1 fastp/!{sample}_R1.fastq.gz \
      --out2 fastp/!{sample}_R2.fastq.gz \
      --unpaired1 fastp/!{sample}_u.fastq.gz \
      --unpaired2 fastp/!{sample}_u.fastq.gz \
      2>> $err_file >> $log_file

    cp *.html fastp/.
    cp *.json fastp/.
  '''
}