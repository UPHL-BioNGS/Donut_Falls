process multiqc {
  publishDir "${params.outdir}", mode: 'copy'
  tag       "multiqc"
  container 'quay.io/biocontainers/multiqc:1.14--pyhdfd78af_0'
  

  //fastp
  //filtlong
  //porechop
  //quast
  //busco

  input:
  file(input)

  output:
  path "multiqc/multiqc_report.html"
  path "multiqc/multiqc_data/*"

  shell:
  '''
    mkdir -p multiqc 
    multiqc --version

    if [ -f "circlator_summary.csv" ]   ; then mv circlator_summary.csv circlator_summary_mqc.csv    ; fi
    if [ -f "flye_summary.tsv" ]        ; then mv flye_summary.tsv flye_summary_mqc.tsv              ; fi
    if [ -f "dragonflye_summary.tsv" ]  ; then mv dragonflye_summary.tsv dragonflye_summary_mqc.tsv  ; fi
    if [ -f "gfastats_summary.csv" ]    ; then mv gfastats_summary.csv gfastats_summary_mqc.csv      ; fi
    if [ -f "NanoStats.csv" ]           ; then mv NanoStats.csv NanoStats_mqc.csv                    ; fi

    multiqc !{params.multiqc_options} \
      --outdir multiqc \
      . 
  '''
}
