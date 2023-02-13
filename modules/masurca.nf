process masurca {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 12
  container 'staphb/masurca:latest'

  input:
  tuple val(sample), file(fastq), file(nanopore)

  output:
  path "masura/${sample}"

  shell:
  '''
  masurca --version

  masurca !{params.masurca_options} \
    -t !{task.cpus} \
    -i !{fastq[0]},!{fastq[1]} \
    -r !{nanopore} 
  
  exit 1
  '''
}

process polca {
  tag "${sample}"
  cpus 6
  container 'staphb/masurca:latest'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(sample), file(fasta), file(fastq)

  output:
  tuple val(sample), file("polca/${sample}_final.fa"),      emit: fasta
  path "polca/${sample}",                                      emit: directory

  shell:
  '''
    masurca --version

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{fasta} \
      -t !{task.cpus} 

    # there was going to be something fancy, but this is here instead
    mkdir round_1
    mv !{fasta}.* round_1/.
    cp round_1/!{fasta}.PolcaCorrected.fa !{sample}_round1.fa

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{sample}_round1.fa \
      -t !{task.cpus}

    mkdir round_2
    mv !{sample}_round1.fa.* round_2/.
    cp round_2/!{sample}_round1.fa.PolcaCorrected.fa !{sample}_round2.fa

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{sample}_round2.fa \
      -t !{task.cpus} 

    mkdir round_3
    mv !{sample}_round2.fa.* round_3/.
    cp round_3/!{sample}_round2.fa.PolcaCorrected.fa !{sample}_round3.fa

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{sample}_round3.fa \
      -t !{task.cpus} 

    mkdir round_4
    mv !{sample}_round3.fa.* round_4/.
    cp round_4/!{sample}_round3.fa.PolcaCorrected.fa !{sample}_round4.fa

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{sample}_round4.fa \
      -t !{task.cpus} 

    mkdir round_5
    mv !{sample}_round4.fa.* round_5/.
    cp round_5/!{sample}_round4.fa.PolcaCorrected.fa !{sample}_round5.fa
    cp round_5/!{sample}_round4.fa.PolcaCorrected.fa polca/!{sample}_final.fa
    
    mv round_* polca/!{sample}/.
  '''
}