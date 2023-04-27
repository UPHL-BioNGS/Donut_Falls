process copy {
  publishDir "${params.outdir}", mode: 'copy'
  tag       "putting all fasta files in ${params.outdir}/consensus"
  container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'

  input:
  file(fasta)

  output:
  path "consensus/"

  shell:
  '''
    mkdir consensus
    
    for fasta in !{fasta}
    do
        cat $fasta | sed 's/_length/ /g' | sed 's/_circular/ /g' > consensus/$fasta
    done
  '''
}
