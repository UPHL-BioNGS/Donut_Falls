process fastp {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/fastp:0.23.2'

  when:
  sample != null

  input:
  tuple val(sample), file(reads)

  output:
  tuple val(sample), file("fastp/${sample}_fastp_{R1,R2}.fastq.gz"), emit: reads
  path "fastp/*"
  path "fastp/${sample}_fastp.json",                                 emit: summary

  shell:
  '''
  mkdir -p fastp 
  fastp --version 

  fastp !{params.fastp_options} \
    --in1 !{reads[0]} \
    --in2 !{reads[1]} \
    --out1 fastp/!{sample}_fastp_R1.fastq.gz \
    --out2 fastp/!{sample}_fastp_R2.fastq.gz \
    --unpaired1 fastp/!{sample}_u.fastq.gz \
    --unpaired2 fastp/!{sample}_u.fastq.gz \
    -h fastp/!{sample}_fastp.html \
    -j fastp/!{sample}_fastp.json
  '''
}