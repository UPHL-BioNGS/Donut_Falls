process subsample {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus 12
    container 'staphb/trycycler:0.5.4'

    input:
    tuple val(sample), file(fastq)

    output:
    tuple val(sample), file("trycycler/subsample/${sample}/*fastq"),       optional: true, emit: fastq
    tuple val(sample), file("trycycler/notsubsample/${sample}/*fastq.gz"), optional: true, emit: full

    shell:
    '''
    mkdir -p trycycler/subsample/!{sample} trycycler/notsubsample/!{sample}

    trycycler --version

    trycycler subsample !{params.trycycler_subsample_options} \
      --reads !{fastq} \
      --threads !{task.cpus} \
      --out_dir trycycler/subsample/!{sample} || \
      cp !{fastq} trycycler/notsubsample/!{sample}/.
    '''
}

process cluster {
  publishDir "${params.outdir}/trycycler", mode: 'copy'
  tag "${sample}"
  cpus 12
  container 'staphb/trycycler:0.5.4'

  input:
  tuple val(sample), file(fasta), file(fastq)

  output:
  path "trycycler/${sample}/contigs.{newick,phylip}"
  path "trycycler/${sample}/cluster*/*contigs/*.fasta"
  tuple val(sample), path("trycycler/${sample}/cluster*"), emit: cluster

  shell:
  '''
    mkdir -p trycycler/!{sample}
    trycycler --version

    trycycler cluster !{params.trycycler_cluster_options} \
      --threads !{task.cpus} \
      --assemblies !{fasta} \
      --reads !{fastq} \
      --out_dir trycycler/!{sample}
  '''
}

process dotplot {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}_${cluster}"
  cpus 1
  container 'staphb/trycycler:0.5.4'

  input:
  tuple val(sample), path(cluster)

  output:
  path("trycycler/dotplot/${sample}_${cluster}_dotplots.png")

  shell:
  '''
    mkdir -p trycycler/dotplot

    trycycler --version

    trycycler dotplot !{params.trycycler_dotplot_options} \
      --cluster_dir !{cluster}

    cp !{cluster}/dotplots.png trycycler/dotplot/!{sample}_!{cluster}_dotplots.png
  '''
}

process reconcile {
  tag "${sample}_${cluster}"
  cpus 12
  //errorStrategy 'finish'
  container 'staphb/trycycler:0.5.4'

  input:
  tuple val(sample), path(cluster), file(fastq), file(remove)

  output:
  tuple val(sample), path(cluster), optional:true, emit: cluster

  shell:
  '''
    trycycler --version 

    if [ -f "!{remove}" ]
    then
      while read line
      do
        cluster=$(echo $line | cut -f 2 -d ,)
        file=$(echo $line | cut -f 3 -d ,)
        if [ -f "$cluster/1_contigs/$file.fasta" ] ; then mv $cluster/1_contigs/$file.fasta $cluster/1_contigs/$file.fasta_remove ; fi
      done < <(grep ^!{sample}, !{remove})
    fi

    num_fasta=$(ls !{cluster}/1_contigs/*.fasta | wc -l)
    echo "There are $num_fasta in !{cluster} for !{sample}"
    if [ "$num_fasta" -ge "!{params.trycycler_min_fasta}" ]
    then
      trycycler reconcile !{params.trycycler_reconcile_options} \
        --reads !{fastq} \
        --cluster_dir !{cluster} \
        --threads !{task.cpus}

        ls

        echo "There are $num_fasta in !{cluster} for !{sample}"

        ls !{cluster}/2_all_seqs.fasta
    else
      echo "!{sample} cluster !{cluster} only had $num_fasta fastas"
      mv !{cluster} !{cluster}_cluster_too_small
    fi
  '''
}

process msa {
  tag "${sample}_${cluster}"
  cpus 12
  container 'staphb/trycycler:0.5.4'

  input:
  tuple val(sample), path(cluster)

  output:
  tuple val(sample), path(cluster), emit: cluster

  shell:
  '''
    trycycler --version

    trycycler msa !{params.trycycler_msa_options} \
      --cluster_dir !{cluster} \
      --threads !{task.cpus}
  '''
}

process partition {
  tag "${sample}_${cluster}"
  cpus 12
  container 'staphb/trycycler:0.5.4'

  input:
  tuple val(sample), path(cluster), file(fastq)

  output:
  tuple val(sample), path(cluster), emit: cluster

  shell:
  '''
    trycycler --version

    trycycler partition !{params.trycycler_partition_options} \
      --reads !{fastq} \
      --cluster_dirs !{cluster} \
      --threads !{task.cpus} 
  '''
}

process consensus {
  publishDir "${params.outdir}/trycycler/${sample}", mode: 'copy'
  tag "${sample}_${cluster}"
  cpus 12
  container 'staphb/trycycler:0.5.4'

  input:
  tuple val(sample), path(cluster)

  output:
  tuple val(sample), path("${cluster}/${sample}_*_trycycler_consensus.fasta"), emit: fasta
  path "${cluster}"

  shell:
  '''
    trycycler --version

    trycycler consensus !{params.trycycler_consensus_options} \
      --cluster_dir !{cluster} \
      --threads !{task.cpus} 
      
    cp !{cluster}/7_final_consensus.fasta !{cluster}/!{sample}_!{cluster}_trycycler_consensus.fasta
  '''
}

process combine {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/trycycler:0.5.4'

  input:
  tuple val(sample), file(fasta)

  output:
  tuple val(sample), path("trycycler/${sample}_trycycler_consensus.fasta"), emit: fasta

  shell:
  '''
    mkdir -p trycycler
    cat !{fasta} > trycycler/!{sample}_trycycler_consensus.fasta
  '''
}