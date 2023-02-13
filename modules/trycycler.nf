process subsample {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus 12
    container 'quay.io/biocontainers/trycycler:0.5.3--pyhdfd78af_0'

    input:
    tuple val(sample), file(fastq)

    output:
    tuple val(sample), file("trycycler/subsample/${sample}/*fastq"), emit: fastq

    shell:
    '''
    mkdir -p trycycler/subsample/!{sample}

    trycycler --version

    trycycler subsample !{params.trycycler_subsample_options} \
      --reads !{fastq} \
      --threads !{task.cpus} \
      --out_dir trycycler/subsample/!{sample}
    '''
}

    // num_lines=$(zcat !{fastq} | wc -l )
    // if [ "$num_lines" -gt 40000 ]
    // then
    //   trycycler subsample !{params.trycycler_subsample_options} \
    //     --reads !{fastq} \
    //     --threads !{task.cpus} \
    //     --out_dir trycycler/subsample/!{sample}
    // else
    //   echo "Not enough reads for Trycycler subsample"
    //   cp !{fastq} trycycler/subsample/!{sample}_1.fastq.gz
    //   cp !{fastq} trycycler/subsample/!{sample}_2.fastq.gz
    //   cp !{fastq} trycycler/subsample/!{sample}_3.fastq.gz
    //   cp !{fastq} trycycler/subsample/!{sample}_4.fastq.gz
    //   cp !{fastq} trycycler/subsample/!{sample}_5.fastq.gz
    //   cp !{fastq} trycycler/subsample/!{sample}_6.fastq.gz
    //   cp !{fastq} trycycler/subsample/!{sample}_7.fastq.gz
    //   cp !{fastq} trycycler/subsample/!{sample}_8.fastq.gz
    //   fi

process cluster {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 12
  container 'quay.io/biocontainers/trycycler:0.5.3--pyhdfd78af_0'

  input:
  tuple val(sample), file(fastas), file(fastq)

  output:
  file("${task.process}/${sample}/contigs.{newick,phylip}")
  file("${task.process}/${sample}/cluster*/*contigs/*.fasta")
  tuple sample, path("${task.process}/${sample}") into inital_clusters

  shell:
  '''
    mkdir -p !{task.process}/!{sample}
    trycycler --version

    trycycler cluster !{params.trycycler_cluster_options} \
      --threads !{task.cpus} \
      --assemblies *.fasta \
      --reads !{fastq} \
      --out_dir !{task.process}/!{sample}
  '''
}

process dotplot {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${cluster}"
  cpus 1
  container 'quay.io/biocontainers/trycycler:0.5.3--pyhdfd78af_0'

  input:
  set val(sample), path(cluster)

  output:
  path("${cluster}/*/dotplots.png")

  shell:
  '''
    mkdir -p dotplot

    trycycler --version

    trycycler dotplot !{params.trycycler_dotplot_options} \
      --cluster_dir !{cluster}
  '''
}

process reconcile {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 12
  //errorStrategy 'finish'
  container 'quay.io/biocontainers/trycycler:0.5.3--pyhdfd78af_0'

  input:
  set val(sample), path(cluster), file(fastq) from clusters.join(fastq, by: 0).view()

  output:
  path("${cluster}")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    ls -d !{cluster}/cluster*/
    for cluster in ls -d !{cluster}/cluster*/
    do
      num_fasta=$(ls $cluster/1_contigs/*.fasta | wc -l)
      if [ "$num_fasta" -ge "!{params.trycycler_reconcile_minimum}" ]
      then
        echo "Working on cluster $cluster" | tee -a $err_file $log_file
        trycycler reconcile !{params.trycycler_reconcile_options} \
          --reads !{fastq} \
          --cluster_dir $cluster \
          --threads !{task.cpus} \
          2>> $err_file >> $log_file
      else
        echo "cluster $cluster only had $num_fasta fastas" >> $log_file
      fi
    done
  '''
}

process msa {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}_${trycycler}"
  cpus 12
  container 'quay.io/biocontainers/trycycler:0.5.3--pyhdfd78af_0'

  input:
  tuple val(sample), path(trycycler)

  output:
  tuple sample, env(cluster), path(trycycler) into msa

  shell:
  '''
    trycycler --version

    cluster=!{trycycler}

    trycycler msa !{params.trycycler_msa_options} \
      --cluster_dir !{trycycler} \
      --threads !{task.cpus}
  '''
}

process partition {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}_${cluster}"
  cpus 12
  container 'quay.io/biocontainers/trycycler:0.5.3--pyhdfd78af_0'

  input:
  set val(sample), val(cluster), path(trycycler), file(fastq) from msa.combine(fastq, by: 0)

  output:
  tuple sample, cluster, path(trycycler) into partition
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file
    trycycler partition !{params.trycycler_partition_options} \
      --reads !{fastq} \
      --cluster_dirs !{trycycler} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file
  '''
}

process consensus {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}_${cluster}"
  cpus 12
  container 'quay.io/biocontainers/trycycler:0.5.3--pyhdfd78af_0'

  input:
  set val(sample), val(cluster), path(trycycler) from partition

  output:
  tuple sample, cluster, path(trycycler) into consensus
  file("${trycycler}/${sample}_${cluster}_trycycler_consensus.fasta") into consensus_all
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    trycycler --version

    trycycler consensus !{params.trycycler_consensus_options} \
      --cluster_dir !{trycycler} \
      --threads !{task.cpus} 
      
    cp !{trycycler}/7_final_consensus.fasta !{trycycler}/!{sample}_!{cluster}_trycycler_consensus.fasta
  '''
}