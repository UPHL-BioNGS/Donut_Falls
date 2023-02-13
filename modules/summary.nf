process summary {
  tag "${sample}"
  cpus 1
  publishDir "${params.outdir}", mode: 'copy'

  input:
  tuple val(sample),
    val(num_contigs),
    val(num_closed_contigs),
    val(chr_closed),
    val(chr_cov),
    val(genome_size),
    val(nanoplot_mean_read_length),
    val(nanoplot_mean_read_quality),
    val(nanoplot_median_read_length),
    val(nanoplot_median_read_quality),
    val(nanoplot_number_of_reads),
    val(nanoplot_read_length_N),
    val(nanoplot_total_bases) from results

  output:
  path "summary/${sample}.summary.tsv",                          emit: summary
  path "logs/summary/${sample}.${workflow.sessionId}.{log,err}", emit: logs

  shell:
  '''
    mkdir -p logs/summary summary
    log_file=logs/summary/!{sample}.!{workflow.sessionId}.log
    err_file=logs/summary/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null

    header="sample\tnumber_of_contigs\tclosed_contigs\tchr_cov\tchr_closed\tgenome_size"
    result="!{sample}\t!{num_contigs}\t!{num_closed_contigs}\t!{chr_cov}\t!{chr_closed}\t!{genome_size}"

    if [ "!{params.nanoplot}" != "false" ]
    then
      header="$header\tmean_read_length\tmean_read_quality\tmedian_read_length\tmediean_read_quality\tnumber_of_reads\tread_length_N50\ttotal_bases"
      result="$result\t!{nanoplot_mean_read_length}\t!{nanoplot_mean_read_quality}\t!{nanoplot_median_read_length}\t!{nanoplot_median_read_quality}\t!{nanoplot_number_of_reads}\t!{nanoplot_read_length_N}\t!{nanoplot_total_bases}"
    fi

    echo -e "$header" > summary/!{sample}.summary.tsv
    echo -e "$result" >> summary/!{sample}.summary.tsv
  '''
}