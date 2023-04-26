process masurca {
  publishDir "${params.outdir}", mode: 'copy'
  stageInMode 'copy'

  tag "${sample}"
  cpus 12
  container 'staphb/masurca:4.1.0'

  input:
  tuple val(sample), file(nanopore), file(fastq)

  output:
  path "masurca/${sample}"
  tuple val(sample), file("masurca/${sample}/${sample}_primary.genome.scf.fasta"), emit: fasta

  shell:
  '''
    mkdir masurca/!{sample}

    masurca --version

    masurca !{params.masurca_options} \
      -t !{task.cpus} \
      -i !{fastq[0]},!{fastq[1]} \
      -r !{nanopore}

    for dir in ls -d CA*
    do
      mv $dir masurca/!{sample}/.
    done

    cp $(ls masurca/!{sample}/CA*/primary.genome.scf.fasta) masurca/!{sample}/!{sample}_primary.genome.scf.fasta 
  '''
}

process polca {
  tag "${sample}"
  cpus 6
  container 'staphb/masurca:latest'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(sample), file(fasta), file(fastq)

  output:
  tuple val(sample), file("polca/${sample}/${sample}_polca_polished.fa"), optional: true, emit: fasta
  path "polca/${sample}/*",                                                               emit: directory

  shell:
  '''
    mkdir -p polca/!{sample}

    masurca --version

    cp !{fasta} !{sample}.fasta

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{sample}.fasta \
      -t !{task.cpus} 

    mv !{sample}.fasta* polca/!{sample}/.

    if [ -f "polca/!{sample}/!{sample}.fasta.PolcaCorrected.fa" ] ; then cp polca/!{sample}/!{sample}.fasta.PolcaCorrected.fa polca/!{sample}/!{sample}_polca_polished.fa ; fi
  '''
}