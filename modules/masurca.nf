process polca {
  tag "${sample}"
  cpus 6

  input:
  tuple val(sample), file(fasta), file(fastq)

  output:
  tuple val(sample), file("polca/${sample}_final.fa"),      emit: fasta
  path "polca/${sample}",                                      emit: directory
  path "logs/polca/${sample}.${workflow.sessionId}.{log,err}", emit: logs

  shell:
  '''
    mkdir -p logs/polca polca/!{sample}
    log_file=logs/polca/!{sample}.!{workflow.sessionId}.log
    err_file=logs/polca/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    masurca --version >> $log_file

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{fasta} \
      -t !{task.cpus} \
      2>> $err_file >> $log_file

    # there was going to be something fancy, but this is here instead
    mkdir round_1
    mv !{fasta}.* round_1/.
    cp round_1/!{fasta}.PolcaCorrected.fa !{sample}_round1.fa

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{sample}_round1.fa \
      -t !{task.cpus} \
      2>> $err_file >> $log_file

    mkdir round_2
    mv !{sample}_round1.fa.* round_2/.
    cp round_2/!{sample}_round1.fa.PolcaCorrected.fa !{sample}_round2.fa

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{sample}_round2.fa \
      -t !{task.cpus} \
      2>> $err_file >> $log_file

    mkdir round_3
    mv !{sample}_round2.fa.* round_3/.
    cp round_3/!{sample}_round2.fa.PolcaCorrected.fa !{sample}_round3.fa

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{sample}_round3.fa \
      -t !{task.cpus} \
      2>> $err_file >> $log_file

    mkdir round_4
    mv !{sample}_round3.fa.* round_4/.
    cp round_4/!{sample}_round3.fa.PolcaCorrected.fa !{sample}_round4.fa

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{sample}_round4.fa \
      -t !{task.cpus} \
      2>> $err_file >> $log_file

    mkdir round_5
    mv !{sample}_round4.fa.* round_5/.
    cp round_5/!{sample}_round4.fa.PolcaCorrected.fa !{sample}_round5.fa
    cp round_5/!{sample}_round4.fa.PolcaCorrected.fa polca/!{sample}_final.fa
    
    mv round_* polca/!{sample}/.
  '''
}