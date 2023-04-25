process flye {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 12
  container 'staphb/flye:2.9.2'
  errorStrategy 'ignore'

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), file("flye/${sample}/${sample}_flye.fasta"), optional: true,  emit: fasta
  tuple val(sample), file("flye/${sample}/${sample}_flye.gfa"),   optional: true,  emit: gfa
  path "flye/${sample}/${sample}_assembly_info.tsv",                               emit: summary
  path "flye/${sample}/*"                                                     

  shell:
  '''
    mkdir -p flye/!{sample}

    flye --version 

    flye !{params.flye_options} \
      --nano-raw !{fastq} \
      --threads !{task.cpus} \
      --out-dir flye/!{sample}

    # renaming final files
    if [ -f "flye/!{sample}/assembly.fasta" ]     ; then cp flye/!{sample}/assembly.fasta     flye/!{sample}/!{sample}_flye.fasta ; fi
    if [ -f "flye/!{sample}/assembly_graph.gfa" ] ; then cp flye/!{sample}/assembly_graph.gfa flye/!{sample}/!{sample}_flye.gfa   ; fi

    # getting a summary file
    head -n 1 flye/!{sample}/assembly_info.txt | awk '{print "sample\\t" $0}' > flye/!{sample}/!{sample}_assembly_info.tsv
    tail -n+2 flye/!{sample}/assembly_info.txt | awk -v sample=!{sample} '{print sample "\\t" $0}' >> flye/!{sample}/!{sample}_assembly_info.tsv
  '''
}