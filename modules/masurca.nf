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
  tuple val(sample), file("masruca/${sample}.fasta"), emit: fasta
  tuple val(sample), file("masruca/${sample}.gfa"), emit: gfa

  shell:
  '''
  masurca --version

# The same error with 4.0.4 with Illumina and nanopore data
# Upd package numactl must be installed.
# see https://github.com/alekseyzimin/masurca/issues/239

  gunzip !{fastq[0]}
  gunzip !{fastq[1]}

  masurca !{params.masurca_options} \
    -t !{task.cpus} \
    -i $(echo !{fastq[0]} | sed 's/.gz$//g'),$(echo !{fastq[1]}  | sed 's/.gz$//g') \
    -r !{nanopore}

  exit 1
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