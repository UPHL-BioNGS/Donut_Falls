process polypolish {
  tag "${sample}"
  cpus 6
  container 'quay.io/biocontainers/polypolish:0.5.0'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(sample), file(fasta), file(sam)

  output:
  tuple val(sample), file("polypolish/${sample}_polypolish.fasta"), emit: fasta

  shell:
  '''
    mkdir -p polypolish

    polypolish_insert_filter.py --in1 !{sam[0]} --in2 !{sam[1]} --out1 !{sample}_filtered_1.sam --out2 !{sample}_filtered_2.sam 
    
    polypolish !{params.polypolish_options} \
      !{fasta} \
      !{sample}_filtered_1.sam \
      !{sample}_filtered_2.sam > polypolish/!{sample}_polypolish.fasta
  '''
}

// in polypolish container python3 -m pip install edlib mappy