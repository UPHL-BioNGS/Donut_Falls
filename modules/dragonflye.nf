process dragonflye {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 12
  container 'staphb/dragonflye:2.9.2'

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), file("flye/${sample}/${sample}_flye.fasta"), optional: true,  emit: fasta
  tuple val(sample), file("flye/${sample}/${sample}_flye.gfa"),   optional: true,  emit: gfa
  path "flye/${sample}/${sample}_assembly_info.csv",                               emit: summary
  path "flye/${sample}/*"                                                     

  shell:
  '''
    mkdir -p dragonflye/!{sample}

    dragonflye --version 

    dragonflye !{params.dragonflye_options} \
      --reads !{fastq} \
      --cpus !{task.cpus} \
      --outdir flye/!{sample} \
      --prefix !{sample}

    # renaming final files
    if [ -f "dragonflye/!{sample}/assembly.fasta" ]     ; then cp dragonflye/!{sample}/assembly.fasta     dragonflye/!{sample}/!{sample}_flye.fasta ; fi
    if [ -f "dragonflye/!{sample}/assembly_graph.gfa" ] ; then cp dragonflye/!{sample}/assembly_graph.gfa dragonflye/!{sample}/!{sample}_flye.gfa   ; fi

    # getting a summary file
    head -n 1 dragonflye/!{sample}/assembly_info.txt | tr "\\t" "," | awk '{print "sample," $0}' > flye/!{sample}/!{sample}_assembly_info.csv
    tail -n+2 dragonflye/!{sample}/assembly_info.txt | tr "\\t" "," | awk -v sample=!{sample} '{print sample "," $0}' >> flye/!{sample}/!{sample}_assembly_info.csv
  '''
}