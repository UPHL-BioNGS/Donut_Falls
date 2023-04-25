process bandage {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'quay.io/biocontainers/bandage:0.8.1--hc9558a2_2'

  input:
  tuple val(sample), file(gfa)

  output:
  tuple val(sample), path("bandage/${sample}.{png,svg}"),   emit: fastq
  path "bandage/${sample}_mqc.png",                         emit: summary

  shell:
  '''
    mkdir -p bandage

    Bandage image !{gfa} bandage/!{sample}.png !{params.bandage_options}
    Bandage image !{gfa} bandage/!{sample}.svg !{params.bandage_options}

    cp bandage/!{sample}.png bandage/!{sample}_mqc.png
  '''
}