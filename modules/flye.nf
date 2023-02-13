process flye {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 12
  container 'staphb/flye:latest'
  //errorStrategy 'ignore'

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), file("flye/${sample}/${sample}_flye.fasta"), optional: true,  emit: fasta
  tuple val(sample), file("flye/${sample}/${sample}_flye.gfa"),   optional: true,  emit: gfa
  tuple val(sample), file("flye/${sample}/assembly_info.txt"),                emit: info
  path "flye/${sample}",                                                      emit: directory

  shell:
  '''
    mkdir -p flye/!{sample}

    flye --version 

    flye !{params.flye_options} \
      --nano-raw !{fastq} \
      --threads !{task.cpus} \
      --out-dir flye/!{sample}

    if [ -f "flye/!{sample}/assembly.fasta" ] ; then cp flye/!{sample}/assembly.fasta flye/!{sample}/!{sample}_flye.fasta ; fi
    if [ -f "flye/!{sample}/assembly.gfa" ]   ; then cp flye/!{sample}/assembly.fasta flye/!{sample}/!{sample}_flye.gfa ; fi
  '''
}