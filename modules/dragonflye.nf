process dragonflye {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 12
  container 'staphb/dragonflye:1.0.14'

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), file("dragonflye/${sample}/${sample}_dragonflye.fasta"), optional: true,  emit: fasta
  tuple val(sample), file("dragonflye/${sample}/${sample}_dragonflye.gfa"),   optional: true,  emit: gfa
  path "dragonflye/${sample}/${sample}_assembly_info.csv",                                     emit: summary
  path "dragonflye/${sample}/*"                                                     

  shell:
  '''
    mkdir -p dragonflye

    dragonflye --version 

    dragonflye !{params.dragonflye_options} \
      --reads !{fastq} \
      --cpus !{task.cpus} \
      --outdir dragonflye/!{sample} \
      --prefix !{sample}

    # renaming final files
    if [ -f "dragonflye/!{sample}/flye-unpolished.gfa" ] ; then cp dragonflye/!{sample}/flye-unpolished.gfa dragonflye/!{sample}/!{sample}_dragonflye.gfa ; fi
    if [ -f "dragonflye/!{sample}/flye.fasta" ] ; then cp dragonflye/!{sample}/flye.fasta dragonflye/!{sample}/!{sample}_dragonflye.fasta ; fi
 
    # getting a summary file
    head -n 1 dragonflye/!{sample}/flye-info.txt | tr "\\t" "," | awk '{print "sample," $0}' > dragonflye/!{sample}/!{sample}_assembly_info.csv
    tail -n+2 dragonflye/!{sample}/flye-info.txt | tr "\\t" "," | awk -v sample=!{sample} '{print sample "," $0}' >> dragonflye/!{sample}/!{sample}_assembly_info.csv
  '''
}