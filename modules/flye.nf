process flye {
  publishDir "donut_falls", mode: 'copy'
  tag "${sample}"
  cpus 12
  container 'staphb/flye:latest'
//  errorStrategy 'ignore'

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), file("flye/${sample}.fasta"), optional: true, emit: fasta
  tuple val(sample), file("flye/${sample}/assembly_info.txt"),     emit: into
  path "flye/${sample}",                                           emit: directory
  path "logs/flye/${sample}.${workflow.sessionId}.{log,err}",      emit: logs

  shell:
  '''
    mkdir -p logs/flye flye/!{sample}
    log_file=logs/flye/!{sample}.!{workflow.sessionId}.log
    err_file=logs/flye/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    flye --version 2>> $err_file >> $log_file

    flye !{params.flye_options} \
      --nano-raw !{fastq} \
      --threads !{task.cpus} \
      --out-dir flye/!{sample} \
      2>> $err_file >> $log_file

    if [ -f "flye/!{sample}/assembly.fasta" ]
    then
      cp flye/!{sample}/assembly.fasta flye/!{sample}.fasta
    fi
  '''
}