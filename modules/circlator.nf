process circlator {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  //container 'quay.io/biocontainers/circlator:1.5.5--py_3'
  //container 'staphb/circlator:latest'
  //container 'sangerpathogens/circlator:latest'

  input:
  tuple val(sample), file(fasta)

  output:
  tuple val(sample), file("circlator/${sample}_unpolished.fasta"),   emit: fasta
  path("circlator/${sample}*"),                                    emit: directory

  shell:
  '''
    mkdir -p circlator

    circlator version

    touch test_circular.fasta
    cat *circular.fasta > circular.fasta

    circlator fixstart !{params.circlator_options} \
        circular.fasta \
        circlator/!{sample}_fixstart

    cp circlator/!{sample}_fixstart.fasta circlator/!{sample}_unpolished.fasta

    touch test_open.fasta
    cat *open.fasta >> circlator/!{sample}_unpolished.fasta
  '''
}