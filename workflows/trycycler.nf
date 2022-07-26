#!/usr/bin/env nextflow

println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.20220430")
println("")

params.clustered_assemblies = workflow.launchDir + '/donut_falls/trycycler_cluster'
println("Directory for clustered assemblies:" + params.clustered_assemblies)

params.outdir = workflow.launchDir + '/three_sisters'

Channel
  .fromPath("${params.clustered_assemblies}/*", type:'dir')
  .ifEmpty {
    println("Could not find fastq pass directories. Set with 'params.fastq_pass_directory'")
    exit 1
  }
//  .map()
  .view()
  .into {clustered_assemblies, clustered_assemblies_dotplot}

params.filtered_fastq = workflow.launchDir + '/donut_falls/rename'
Channel
  .fromPath("params.filtered_fastq\*fastq", type:'file')
  .ifEmpty {
    println("Could not find fastq files. Set with 'params.fastq_pass_directory'")
    exit 1
  }
//  .map()
  .view()
  .set{filtered_fastq}

params.trycycler_dotplot_options = ''
process trycycler_dotplot {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo true
  cpus 1
  //  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), path(cluster), from clustered_assemblies_dotplot

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into dotplots
  file("logs/${task.process}/*.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    trycycler dotplot !{params.trycycler_dotplot_options} \
      --cluster_dir !{cluster} \
      2>> $err_file >> $log_file

    exit 1
  '''
}

clustered_assemblies
  .join(filtered_fastq, by: 0)
  .join(dotplots, by:0)
  .view()
  .set{assemblies_and_fastq}

process trycycler_reconcile {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo true
  cpus 1
//  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), path(cluster), file(fastq), file(dotplots) from fastq_pass_directory

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into reconciled_fastas
  file("logs/${task.process}/*.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    trycycler reconcile \
      --reads !{fastq} \
      --cluster_dir !{cluster} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file

    exit 1
  '''
}

params.trycycler_msa_options = ''
process trycycler_msa {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo true
  cpus 1
//  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), file(fasta) from reconciled_fastas

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into msa
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    trycycler msa !{params.trycycler_msa_options} \
      --cluster_dir !{cluster} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file

      exit 1
  '''
}

params.trycycler_partition_options = ''
process trycycler_partition {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo true
  cpus 1
//  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), file(fasta), file(fastq) from reconciled_fastas

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into msa
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    trycycler partition !{trycycler_partition_options} \
      --reads !{fastq} \
      --cluster_dirs !{cluster} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file

    exit 1
  '''
}


params.trycycler_consensus_options = ''
process trycycler_consensus {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo true
  cpus 1
//  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), file(fasta), file(fastq) from reconciled_fastas

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into msa
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    trycycler consensus !{params.trycycler_consensus_options} \
      --cluster_dir !{cluster} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file

    exit 1
  '''
}

params.medaka_options = ''
process medaka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo true
  cpus 1
//  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), file(fasta), file(fastq) from reconciled_fastas

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into medaka_fastas
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    medaka --version >> $log_file

    exit 1
  '''
}


workflow.onComplete {
    println("Pipeline completed at: $workflow.complete")
    println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
