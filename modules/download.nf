process download {
  tag "${sra}"

  input:
  val(sra)

  output:
  tuple val(sra), file("sra/${sra}*fastq.gz"), emit: fastq

  shell:
  '''
    mkdir -p sra/

    sra=!{sra}

    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sra:0:6}/0${sra: -2}/!{sra}/!{sra}.fastq.gz

    if [ -f "!{sra}.fastq.gz" ] ; then mv !{sra}.fastq.gz sra/. ; fi
  '''
}