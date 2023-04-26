process medaka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 6
  container 'ontresearch/medaka:v1.7.3'

  input:
  tuple val(sample), path(fasta), path(fastq)

  output:
  path "medaka/${sample}/",                                                     emit: directory
  tuple val(sample), path("medaka/${sample}/${sample}_medaka_consensus.fasta"), emit: fasta

  shell:
  '''
    mkdir -p medaka
    
    medaka --version

    cat !{fasta} | sed 's/ /_/g' > !{fasta}.fasta

    medaka_consensus !{params.medaka_options} \
      -i !{fastq} \
      -d !{fasta}.fasta \
      -o medaka/!{sample} \
      -t !{task.cpus}

    cp medaka/!{sample}/consensus.fasta medaka/!{sample}/!{sample}_medaka_consensus.fasta
  '''
}