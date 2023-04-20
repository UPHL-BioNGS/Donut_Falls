process bgzip {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/htslib:1.17'

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), path("filtlong/${fastq}.gz"), emit: fastq

  shell:
  '''
    mkdir -p filtlong 
    bgzip --version 

    bgzip -@ !{task.cpus} !{fastq}
    mv !{fastq}.gz filtlong/.
  '''
}