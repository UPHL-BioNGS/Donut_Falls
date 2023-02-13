process porechop {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 6
  container 'quay.io/biocontainers/porechop:0.2.4--py310h30d9df9_3'

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), file("porechop/${sample}_chopped.fastq.gz"), emit: fastq

  shell:
  '''
    mkdir -p porechop 

    porechop --version

    porechop !{params.porechop_options} \
        --threads !{task.cpus} \
        -i !{fastq} \
        -o porechop/!{sample}_chopped.fastq.gz  
  '''
}