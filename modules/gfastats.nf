process gfastats {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'quay.io/biocontainers/gfastats:1.3.6--hd03093a_1'

  input:
  tuple val(sample), file(gfa)

  output:
  tuple val(sample), file("gfastats/${sample}_*.fasta"), optional: true, emit: fasta
  path "gfastats/*"

  shell:
  '''
  mkdir -p gfastats
  gfastats --version 

  gfastats  \
    !{gfa} \
    !{params.gfastats_options} \
    --threads !{task.cpus} \
    --tabular \
    --seq-report > gfastats/!{sample}_gfastats.txt

  while read line
  do
    header=$(echo $line | cut -f 2 -d ',')
    echo "header is $header"
    length=$(echo $line | cut -f 4 -d ',')
    echo "length is $length"
    circ=$(echo $line | cut -f 11 -d ',')
    echo "circ is $circ"
    if [ "$length" -ge 200 ]
    then
      if [[ "$circ" == "Y" ]]
      then
        echo ">$header length=$length circular=true" > gfastats/!{sample}_${header}_circular.fasta
        grep -w "^S" !{gfa} | grep -w $header | awk '{print $3}' >> gfastats/!{sample}_${header}_circular.fasta
      else
        echo ">$header length=$length circular=false" > gfastats/!{sample}_${header}_open.fasta
        grep -w "^S" !{gfa} | grep -w $header | awk '{print $3}' >> gfastats/!{sample}_${header}_open.fasta
      fi
    fi
  done < <(grep -v Header gfastats/!{sample}_gfastats.txt | tr "\\t" ",")
  '''
}
