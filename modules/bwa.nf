process bwa {
  tag "${sample}"
  label "maxcpus"
  container 'staphb/bwa:latest'
  cpus 6

  input:
  tuple val(sample), file(fasta), file(fastq)

  output:
  tuple val(sample), file(fasta), file("bwa/${sample}_{1,2}.sam"), emit: sam

  shell:
  '''
    mkdir -p bwa

    bwa index !{fasta}
    bwa mem -t !{task.cpus} -a !{fasta} !{fastq[0]} > bwa/!{sample}_1.sam
    bwa mem -t !{task.cpus} -a !{fasta} !{fastq[1]} > bwa/!{sample}_2.sam
  '''
}
