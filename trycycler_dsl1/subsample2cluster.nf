#!/usr/bin/env nextflow

println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.0.20220430")
println("")

params.outdir = workflow.launchDir + '/trycycler_1'

if ( Runtime.runtime.availableProcessors() < 12 ) {
  params.maxcpus = Runtime.runtime.availableProcessors()
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 6
  params.maxcpus = 12
}

params.reads = workflow.launchDir + '/reads'
Channel
  .fromPath("${params.reads}/*.{fastq,fastq.gz,fq,fq.gz}", type:'file')
  .ifEmpty {
    println("Could not find fastq files for nanopore sequencing. Set with 'params.reads'")
    exit 1
  }
  .map { reads -> tuple(reads.simpleName, reads ) }
  .view { "Fastq file found : ${it[0]}" }
  .set { fastq }

params.filtlong_options = "--min_length 1000 --keep_percent 95"
process filtlong {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/filtlong:latest'

  input:
  tuple val(sample), file(fastq) from fastq

  output:
  tuple val(sample), path("${task.process}/${sample}_filtered.fastq") into filtered_fastq, filtered_fastq_cluster
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
    """
      mkdir -p ${task.process} logs/${task.process}
      log_file=logs/${task.process}/${sample}.${workflow.sessionId}.log
      err_file=logs/${task.process}/${sample}.${workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a \$log_file \$err_file > /dev/null
      filtlong --version >> \$log_file

      filtlong ${params.filtlong_options} \\
        ${fastq} \\
        2>> \$err_file > ${task.process}/${sample}_filtered.fastq
      """
    }

params.trycycler_subsample_options = ""
params.trycycler_subsample_count = "12"
if ( params.trycycler_subsample_count.toInteger() % 3 > 0 ) {
  println("For this workflow, 'params.trycycler_subsample_count' (which you set as ${params.trycycler_subsample_count}) needs to be divisible by three. Sorry!")
  exit 1
}
subsample_count = params.trycycler_subsample_count.toInteger() / 3
twice_count = subsample_count * 2
process trycycler_subsample {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.maxcpus
  container 'staphb/trycycler:latest'
  //errorStrategy 'ignore' // because this will attempt to subsample even when there aren't enough reads

  input:
  set val(sample), file(fastq) from filtered_fastq

  output:
  tuple sample, file("${task.process}/${sample}/miniasm/sample*.fastq") into subsampled_miniasm_fastq
  tuple sample, file("${task.process}/${sample}/flye/sample*fastq") into subsampled_flye_fastq
  tuple sample, file("${task.process}/${sample}/raven/sample*fastq") into subsampled_raven_fastq
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample}/{miniasm,flye,raven} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version 2>> $err_file >> $log_file
    trycycler subsample !{params.trycycler_subsample_options} \
      --reads !{fastq} \
      --count !{params.trycycler_subsample_count} \
      --threads !{task.cpus} \
      --out_dir !{task.process}/!{sample} \
      2>> $err_file >> $log_file

    ls !{task.process}/!{sample}/*fastq  | head -n !{subsample_count} | xargs -I {} cp {} !{task.process}/!{sample}/miniasm/.
    ls !{task.process}/!{sample}/*fastq  | tail -n !{subsample_count} | xargs -I {} cp {} !{task.process}/!{sample}/flye/.
    ls !{task.process}/!{sample}/*fastq  | tail -n !{twice_count} | head -n !{subsample_count} | xargs -I {} cp {} !{task.process}/!{sample}/raven/.
  '''
}

params.flye_options = '--plasmids'
process flye {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/flye:latest'

  input:
  set val(sample), file(fastq) from subsampled_flye_fastq.transpose()

  output:
  tuple sample, file("${task.process}/${sample}/${sample}_*_${task.process}.fasta") into flye_fastas
  file("${task.process}/${sample}/*/{assembly_graph.gfa,assembly_graph.gv,assembly_info.txt,flye.log,params.json}")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    flye --version 2>> $err_file >> $log_file

    name=$(echo !{fastq} | cut -f 2 -d "_" | cut -f 1 -d ".")
    mkdir -p !{task.process}/!{sample}/$name
    flye !{params.flye_options} \
      --nano-raw !{fastq} \
      --threads !{task.cpus} \
      --out-dir !{task.process}/!{sample}/$name \
      2>> $err_file >> $log_file
      cp !{task.process}/!{sample}/$name/assembly.fasta !{task.process}/!{sample}/!{sample}_${name}_!{task.process}.fasta
  '''
}

params.miniasm_and_minipolish_options = ''
process miniasm_and_minipolish {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/minipolish:latest'

  input:
  set val(sample), file(fastq) from subsampled_miniasm_fastq.transpose()

  output:
  tuple sample, file("${task.process}/${sample}/*gfa") into miniasm_gfa
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "miniasm version : $(miniasm -V)" 2>> $err_file >> $log_file

    miniasm_and_minipolish.sh --version 2>> $err_file >> $log_file

    name=$(echo !{fastq} | cut -f 2 -d "_" | cut -f 1 -d ".")
    miniasm_and_minipolish.sh \
      !{fastq} \
      !{task.cpus} \
      2>> $err_file \
      > !{task.process}/!{sample}/!{fastq}_${name}.gfa
  '''
}

params.any2fasta_options = ''
process any2fasta {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  container 'staphb/any2fasta:latest'

  input:
  tuple val(sample), file(gfa) from miniasm_gfa

  output:
  tuple sample, file("${task.process}/${sample}/*_miniasm.fasta") into miniasm_fastas
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    any2fasta -v 2>> $err_file >> $log_file

    name=$(echo !{gfa} | cut -f 2 -d "_" | cut -f 1 -d ".")
    any2fasta !{params.any2fasta_options} \
      !{gfa} \
      2>> $err_file \
      > !{task.process}/!{sample}/!{sample}_${name}_miniasm.fasta
  '''
}

params.raven_options = '--polishing-rounds 2'
process raven {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/raven:latest'

  input:
  set val(sample), file(fastq) from subsampled_raven_fastq.transpose()

  output:
  tuple sample, file("${task.process}/${sample}/*.fasta") into raven_fastas
  file("${task.process}/${sample}/*.gfa")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    raven --version 2>> $err_file >> $log_file

    name=$(echo !{fastq} | cut -f 2 -d "_" | cut -f 1 -d ".")
    raven !{params.raven_options} \
      --threads !{task.cpus} \
      --graphical-fragment-assembly !{task.process}/!{sample}/!{sample}_${name}.gfa \
      !{fastq} \
      2>> $err_file \
      > !{task.process}/!{sample}/!{sample}_${name}_!{task.process}.fasta
  '''
}

miniasm_fastas
  .concat(flye_fastas)
  .concat(raven_fastas)
  .groupTuple(by : 0)
  .join(filtered_fastq_cluster, by:0)
  .set {assembled_fastas}

params.trycycler_cluster_options = ""
process trycycler_cluster {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/trycycler:latest'

  input:
  set val(sample), file(fastas), file(fastq) from assembled_fastas

  output:
  file("${task.process}/${sample}/contigs.{newick,phylip}")
  file("${task.process}/${sample}/cluster*/*contigs/*.fasta")
  tuple sample, path("${task.process}/${sample}") into inital_clusters
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version 2>> $err_file >> $log_file

    trycycler cluster !{params.trycycler_cluster_options} \
      --threads !{task.cpus} \
      --assemblies *.fasta \
      --reads !{fastq} \
      --out_dir !{task.process}/!{sample} \
      2>> $err_file >> $log_file
  '''
}

params.trycycler_dotplot_options = ''
process trycycler_dotplot {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${cluster}"
  echo false
  cpus 1
  container 'staphb/trycycler:latest'

  input:
  set val(sample), path(cluster) from inital_clusters.transpose()

  output:
  path("${cluster}/*/dotplots.png")
  file("logs/${task.process}/${cluster}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{cluster}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{cluster}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file


    trycycler dotplot !{params.trycycler_dotplot_options} \
      --cluster_dir !{cluster} \
      2>> $err_file >> $log_file
  '''
}
