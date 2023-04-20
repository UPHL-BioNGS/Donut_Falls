process download {
  tag "Downloading subset15000"

  output:
  tuple val("subset15000"), file("nfcore_subset15000.fa.gz"), emit: fastq

  shell:
  '''
  wget -q https://github.com/nf-core/test-datasets/blob/23f5b889e4736798c8692e9b92810d9a3e37ee97/nanopore/subset15000.fq.gz?raw=true -O nfcore_subset15000.fa.gz
  '''
}