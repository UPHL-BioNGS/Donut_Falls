process busco {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'ezlabgva/busco:v5.4.5_cv1'

  input:
  tuple val(sample), file(fasta)

  output:
  path("busco/${sample}/*")
  path("busco/${sample}/short_summary_*.txt"), emit: summary

  shell:
  '''
    mkdir -p busco/!{sample}
    busco --version 

    busco !{params.busco_options} \
        -m genome \
        -i !{fasta} \
        -o busco/!{sample} \
        --cpu !{task.cpus} \
        --auto-lineage-prok
  '''
}