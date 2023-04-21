process download {
  tag "Downloading subset15000"
  cpus 1
  container 'staphb/gfastats:1.3.6'

  output:
  tuple val("subset15000"), file("nfcore_subset15000.fa.gz"), emit: fastq

  shell:
  '''
  wget -q https://github.com/nf-core/test-datasets/blob/23f5b889e4736798c8692e9b92810d9a3e37ee97/nanopore/subset15000.fq.gz?raw=true -O nfcore_subset15000.fa.gz
  '''
}

process great_dataset {
  tag "Downloading the great dataset"
  cpus 1
  container 'staphb/gfastats:1.3.6'

  output:
  tuple val("great_dataset"), file("reads.fastq.gz"), emit: fastq

  shell:
  '''
    wget -q https://bridges.monash.edu/ndownloader/files/23754659 -O great_dataset.tar.gz
    tar -xvf great_dataset.tar.gz
  '''
}