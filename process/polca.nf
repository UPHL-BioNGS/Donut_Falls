#!/usr/bin/env nextflow

params.max_polish_rounds = 10
params.polca_options = ''
process polca {
  publishDir "${params.outdir}/${task.process}/${sample}/round_${round}", mode: 'copy', pattern: "*{err,batches,success,vcf,names,report,PolcaCorrected.fa}"
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/${sample}_${round}.${workflow.sessionId}.{log,err}"
  publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/${sample}/${sample}_*.fasta"
  tag "${sample} : round ${round} : changes ${changes}"
  cpus params.medcpus
//  container 'quay.io/uphl/masurca:latest'

  input:
  tuple val(sample), path(fasta), path(fastq), val(round), val(changes) from for_polca

  output:
  tuple val(sample), path("next/${sample}_${round}.fasta"), path(fastq), env(next_round), env(change_test) optional true into new_rounds
  path("*")
  path("${task.process}/${sample}/${sample}_{${round},final}.fasta") optional true
  path("logs/${task.process}/${sample}_${round}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process} !{task.process}/!{sample} next
    log_file=logs/!{task.process}/!{sample}_!{round}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}_!{round}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    masurca --version >> $log_file

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{fasta} \
      -t !{task.cpus}

    if [ -f "!{fasta}.report" ]
    then
      sub_err=$(grep "Substitution Errors" !{fasta}.report | awk '{print $3}')
      del_err=$(grep "Deletion Errors" !{fasta}.report | awk '{print $3}')
      if [ -z "$sub_err" ] ; then sub_err=0 ; fi
      if [ -z "$del_err" ] ; then del_err=0 ; fi
      change_test=$(( sub_err + del_err ))
    else
      echo "WARNING : report not found" >> $err_file
      sub_err=0
      del_err=0
      change_test=0
    fi

    echo "The number of changes from this round was $change_test" >> $log_file
    echo "$sub_err changes were from substitution errors" >> $log_file
    echo "$del_err changes were from insertion/deletion errors" >> $log_file

    if [ "$change_test" -lt "1" ]
    then
      next_round=1000
      cp !{fasta}.PolcaCorrected.fa !{task.process}/!{sample}/!{sample}_!{round}.fasta
      cp !{fasta}.PolcaCorrected.fa !{task.process}/!{sample}/!{sample}_final.fasta
    elif [ "$change_test" -eq "!{changes}" ] && [ "$change_test" -lt "10000" ]
    then
      next_round=1000
      cp !{fasta}.PolcaCorrected.fa !{task.process}/!{sample}/!{sample}_!{round}.fasta
      cp !{fasta}.PolcaCorrected.fa !{task.process}/!{sample}/!{sample}_final.fasta
    else
      cp !{fasta}.PolcaCorrected.fa !{task.process}/!{sample}/!{sample}_!{round}.fasta

      next_round=$(( !{round} + 1 ))
      if [ "$next_round" -le "!{params.max_polish_rounds}" ] ; then cp !{fasta}.PolcaCorrected.fa next/!{sample}_!{round}.fasta ; fi
    fi
  '''
}
