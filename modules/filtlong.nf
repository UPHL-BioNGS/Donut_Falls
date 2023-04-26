process filtlong {
  tag "${sample}"
  cpus 1
  container 'staphb/filtlong:0.2.1'

  input:
  tuple val(sample), file(fastq), file(short_reads)

  output:
  tuple val(sample), file("filtlong/${sample}_filtered.fastq"), optional: true, emit: fastq

  shell:
  if (short_reads[1] == null) {
  '''
    mkdir -p filtlong 

    filtlong --version

    filtlong !{params.filtlong_options} \
      !{fastq} \
      > filtlong/!{sample}_filtered.fastq
  '''
  } else {
  '''
    mkdir -p filtlong 
    
    filtlong --version

    filtlong !{params.filtlong_options} \
      -1 !{short_reads[0]} \
      -2 !{short_reads[1]} \
      !{fastq} \
      > filtlong/!{sample}_filtered.fastq
    '''
  }
}