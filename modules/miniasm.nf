  process miniasm {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 12
  container 'staphb/minipolish:0.1.3'

  input:
  tuple val(sample), file(fastq)

  output:
  tuple val(sample), path("miniasm/${sample}/${sample}_miniasm.gfa"), emit: gfa
  path "miniasm/${sample}/*",                                    emit: directory

  shell:
  '''
    mkdir -p miniasm/!{sample} 
  
    echo "miniasm version : $(miniasm -V)" 
    minimap2 --version 
    minipolish --version 

    miniasm_and_minipolish.sh \
      !{fastq} \
      !{task.cpus} \
      > miniasm/!{sample}/!{sample}_miniasm.gfa
  '''
}

