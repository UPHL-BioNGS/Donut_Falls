process unicycler {
  tag "${sample}"
  cpus 12
  publishDir "${params.outdir}", mode: 'copy'
  container 'staphb/unicycler:latest'

  input:
  tuple val(sample), file(nanopore), file(illumina)

  output:
  path "unicycler/${sample}",                                    emit: directory
  tuple val(sample), file("unicycler/${sample}/assembly.fasta"), emit: fasta

  shell:
  '''
    mkdir -p unicycler 

    unicycler --version

    unicycler !{params.unicycler_options} \
      -1 !{illumina[0]} \
      -2 !{illumina[1]} \
      -l !{nanopore} \
      -o unicycler/!{sample} \
      -t 20 

    if [ -f "unicycler/!{sample}/assembly.fasta" ] ; then cp unicycler/!{sample}/assembly.fasta unicycler/!{sample}/!{sample}.fasta ; fi 
    '''
}

process unicycler_long {
  tag "${sample}"
  cpus 12
  publishDir "${params.outdir}", mode: 'copy'
  container 'staphb/unicycler:latest'

  input:
  tuple val(sample), file(nanopore)

  output:
  path "unicycler/${sample}",                                    emit: directory
  tuple val(sample), file("unicycler/${sample}/assembly.fasta"), emit: fasta
  tuple val(sample), file("unicycler/${sample}/assembly.gfa"),   emit: gfa

  shell:
  '''
    mkdir -p unicycler 

    unicycler --version 

    unicycler !{params.unicycler_options} \
      -l !{nanopore} \
      -o unicycler/!{sample} \
      -t !{task.cpus}

    if [ -f "unicycler/!{sample}/assembly.fasta" ] ; then cp unicycler/!{sample}/assembly.fasta unicycler/!{sample}/!{sample}.fasta ; fi
    if [ -f "unicycler/!{sample}/assembly.gfa" ] ; then cp unicycler/!{sample}/assembly.gfa unicycler/!{sample}/!{sample}.gfa ; fi
    '''
}