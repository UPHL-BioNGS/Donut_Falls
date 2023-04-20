process circlator {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  //stageInMode 'copy'
  container 'staphb/circlator:1.5.5'
  //errorStrategy 'ignore'

  input:
  tuple val(sample), file(fasta)

  output:
  tuple val(sample), file("circlator/${sample}_unpolished.fasta"), emit: fasta
  path "circlator/${sample}*",                                     emit: directory
  path "circlator/${sample}_fixstart_summary.csv",                 emit: summary

  shell:
  '''
    mkdir -p circlator

    circlator version

    nucmer --version

    touch test_circular.fasta
    cat *circular.fasta > circular.fa

    circlator fixstart !{params.circlator_options} \
        circular.fa \
        circlator/!{sample}_fixstart

    cp circlator/!{sample}_fixstart.fasta circlator/!{sample}_unpolished.fasta

    touch test_open.fasta
    cat *open.fasta >> circlator/!{sample}_unpolished.fasta

    head -n 1 circlator/!{sample}_fixstart.log | tr "\\t" "," | awk '{print "sample," $0 }' > circlator/!{sample}_fixstart_summary.csv
    tail -n+2 circlator/!{sample}_fixstart.log | tr "\\t" "," | awk -v sample=!{sample} '{print sample "," $0 }' >> circlator/!{sample}_fixstart_summary.csv
  '''
}