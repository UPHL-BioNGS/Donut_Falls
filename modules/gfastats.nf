process gfastats {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/gfastats:1.3.6'

  input:
  tuple val(sample), file(gfa)

  output:
  tuple val(sample), file("gfastats/${sample}_*_{circular,open}.fasta"), optional: true, emit: fasta
  tuple val(sample), file("gfastats/${sample}.fasta"),                   optional: true, emit: assembly
  path "gfastats/*"
  path "gfastats/${sample}_gfastats_summary.csv",                                      emit: summary

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
    length=$(echo $line | cut -f 4 -d ',')
    circ=$(echo $line | cut -f 11 -d ',')

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
    echo ">!{sample}_${header} length=$length circular=$circ" >> gfastats/!{sample}.fasta
    grep -w "^S" !{gfa} | grep -w $header | awk '{print $3}' >> gfastats/!{sample}.fasta
  done < <(grep -v Header gfastats/!{sample}_gfastats.txt | tr "\\t" ",")

  head -n 1 gfastats/!{sample}_gfastats.txt | tr "\\t" "," | awk '{print "sample," $0 "circular" }' > gfastats/!{sample}_gfastats_summary.csv
  tail -n+2 gfastats/!{sample}_gfastats.txt | tr "\\t" "," | awk -v sample=!{sample} '{print sample "," $0 }' >> gfastats/!{sample}_gfastats_summary.csv
  '''
}
