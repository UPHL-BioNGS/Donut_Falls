#/usr/bin/bash

VERSION=0.0.20230221

USAGE="""
Admittedly, this is probably UPHL specific and not of use to others.
This script is meant to take a tab-delimited or comma-delimited file
1) concatenate fastq.gz files into one file
2) rename the file
3) find Illumina reads (if they exist)
4) create the sample sheet for donut falls

In order for this script to work, there needs to be file that can be used for a sample key. A header is not required, but the order matters.
lab accession, barcode, alternate id, illumina run

USAGE
uphl_sample_sheet.sh \\
\t-f <sample key> \\
\t-o <directory for combined and renamed files> \\
\t-p <path to fastq_pass directory> \\
\t-j <number of jobs for parallel>
"""

combined="combined"
wgs="/Volumes/IDGenomics_NAS/pulsenet_and_arln"
sample_key="samples.csv"
pass="fastq_pass"

while getopts 'o:f:i:j:p:w:hv' OPTION
do
  case "$OPTION" in
    o)
      combined=$OPTARG
    ;;
    f)
      sample_key=$OPTARG
    ;;
    p)
      pass=$OPTARG
    ;;
    j)
      jobs="-j $OPTARG"
    ;;
    h)
      echo -e "$USAGE"
      exit 0
    ;;
    v)
      echo "combine_rename version: $VERSION"
      exit 0
    ;;
    :)
      echo "Invalid option: $OPTARG requires an argument"
      echo -e "$USAGE"
      exit 1
    ;;
    \?)
      echo -e "$USAGE"
      exit 1
    ;;
  esac
done
shift "$(($OPTIND -1))"


mkdir -p $combined
if [ -d "$combined" ]
then
  echo "$(date): The combined and renamed files will be placed in $combined"
else
  echo "$(date): FATAL : output directory $combined could not be created. Set with -o <outdir>"
  exit 1
fi

if [ -d "$pass" ]
then
  echo "$(date): The corresponding directory to look for nanopore reads is $pass"
else
  echo "$(date): FATAL : output directory $pass was not found. Set with -p <dir>"
  exit 1
fi

if [ -f "$sample_key" ]
then
  echo "$(date): The file with relevant ids and run information is $sample_key"
else
  echo "$(date): FATAL : no sample key could be found. Set with -f <file>"
  exit 1
fi

if [ -z "$(which parallel)" ] ; then echo "$(date) : FATAL : parallel is not in PATH" ; exit 1 ; fi

echo "$(date): Concenating fastq.gz files and renaming"
cat $sample_key | parallel --colsep "," $jobs "cat $pass/{2}/*fastq.gz > $combined/{1}.fastq.gz"

echo "sample,fastq,fastq_1,fastq_2" > sample_sheet.csv

echo "$(date): getting R1 and R2 for each sample"
while read line
do
    # lab accession, barcode, alternate id, illumina run
    labid=$(echo $line | cut -f 1 -d , )
    brcde=$(echo $line | cut -f 2 -d , )
    altid=$(echo $line | cut -f 3 -d , )
    ilrun=$(echo $line | cut -f 4 -d , )
    R1=$(find $wgs/$ilrun -name "$labid*R1*fastq.gz" -o -name "$altid*R1*fastq.gz" | head -n 1 )
    R2=$(find $wgs/$ilrun -name "$labid*R2*fastq.gz" -o -name "$altid*R2*fastq.gz" | head -n 1 )
    if [ -z "$R1" ]
    then
      echo "Not found in $wgs/$ilrun"
      echo "Searching /Volumes/IDGenomics_NAS/pulsenet_and_arln/old_runs/2023/$ilrun"
      R1=$(find /Volumes/IDGenomics_NAS/pulsenet_and_arln/old_runs/2023/$ilrun -name "$labid*$ilrun*R1*fastq.gz" -o -name "$altid*$ilrun*R1*fastq.gz" | head -n 1 )
      R2=$(find /Volumes/IDGenomics_NAS/pulsenet_and_arln/old_runs/2023/$ilrun -name "$labid*$ilrun*R2*fastq.gz" -o -name "$altid*$ilrun*R2*fastq.gz" | head -n 1 )
    fi
    if [ -z "$R1" ]
    then
      echo "Searching /Volumes/IDGenomics_NAS/pulsenet_and_arln/old_runs/2022/$ilrun"
      R1=$(find /Volumes/IDGenomics_NAS/pulsenet_and_arln/old_runs/2022/$ilrun -name "$labid*$ilrun*R1*fastq.gz" -o -name "$altid*$ilrun*R1*fastq.gz" | head -n 1 )
      R2=$(find /Volumes/IDGenomics_NAS/pulsenet_and_arln/old_runs/2022/$ilrun -name "$labid*$ilrun*R2*fastq.gz" -o -name "$altid*$ilrun*R2*fastq.gz" | head -n 1 )
    fi
    if [ -z "$R1" ]
    then
      echo "Searching /Volumes/IDGenomics_NAS/pulsenet_and_arln/old_runs/2021/$ilrun"
      R1=$(find /Volumes/IDGenomics_NAS/pulsenet_and_arln/old_runs/2021/$ilrun -name "$labid*$ilrun*R1*fastq.gz" -o -name "$altid*$ilrun*R1*fastq.gz" | head -n 1 )
      R2=$(find /Volumes/IDGenomics_NAS/pulsenet_and_arln/old_runs/2021/$ilrun -name "$labid*$ilrun*R2*fastq.gz" -o -name "$altid*$ilrun*R2*fastq.gz" | head -n 1 )
    fi
    if [ -z "$R1" ]
    then
      echo "Could not find Illumina reads!"
      R1=""
      R2=""
    fi
    echo "$labid,$combined/$labid.fastq.gz,$R1,$R2" >> sample_sheet.csv
    echo "$(date): Information for $labid:"
    echo "$(date): Barcode is $brcde"
    echo "$(date): Combined nanopore fastq is $combined/$labid.fastq.gz"
    echo "$(date): Illumina reads are $R1 and $R2"
done <  $sample_key

echo "$(date): Everything is ready for Donut Falls!"
echo "$(date): Run with \"nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --assembler flye,unicycler --sample_sheet sample_sheet.csv\""

exit 0