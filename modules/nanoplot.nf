process nanoplot_summary {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sequencing_summary}"
  cpus 6
  container 'staphb/nanoplot:latest'

  input:
  file(sequencing_summary)

  output:
  path "nanoplot/summary",                                      emit: final_directory                                  

  shell:
  '''
    mkdir -p nanoplot/summary 
    NanoPlot --version

    NanoPlot !{params.nanoplot_summary_options} \
      --summary !{sequencing_summary} \
      --threads !{task.cpus} \
      --outdir nanoplot/summary \
      --tsv_stats
  '''
}

process nanoplot {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 6
  container 'staphb/nanoplot:latest'

  input:
  tuple val(sample), file(fastq)

  output:
  path "nanoplot/${sample}",                                      emit: final_directory
  path "nanoplot/${sample}/${sample}_NanoStats.csv",              emit: summary

  shell:
  '''
    mkdir -p nanoplot/!{sample}

    NanoPlot --version

    NanoPlot !{params.nanoplot_options} \
      --fastq !{fastq} \
      --threads !{task.cpus} \
      --tsv_stats \
      --outdir nanoplot/!{sample}

    cp nanoplot/!{sample}/NanoStats.txt nanoplot/!{sample}/!{sample}_NanoStats.txt

    echo "sample,$(cut -f 1 nanoplot/!{sample}/!{sample}_NanoStats.txt | tr '\\n' ',' )" >  nanoplot/!{sample}/!{sample}_NanoStats.csv
    echo "!{sample},$(cut -f 2 nanoplot/!{sample}/!{sample}_NanoStats.txt | tr '\\n' ',' )" >> nanoplot/!{sample}/!{sample}_NanoStats.csv
  '''
}
