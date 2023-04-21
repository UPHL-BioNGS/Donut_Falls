process rasusa {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 6
  container 'staphb/rasusa:0.7.0'

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), file("rasusa/${sample}/*.fastq"), emit: fastq

  shell:
  '''
    mkdir -p rasusa/!{sample}

    rasusa --version

    for i in 01 02 03 04 05 06 07 08 09 10 11 12
    do
      rasusa !{params.rasusa_options} \
        -i !{fastq} \
        --output rasusa/!{sample}/sample_$i.fastq
    done


  '''
}