#/usr/bin/bash

VERSION=0.0.20220429

USAGE="""
Admittedly, this is probably UPHL specific and not of use to others.

This script is meant to take a tab-delimited or comma-delimited file
1) use bgzip to compress fastq file to fastq.gz if not already compressed
2) concatenate fastq.gz files into one file
3) rename the file
4) copy the corresponding Illumina reads

In order for this script to work, there needs to be file that can be used for a sample key. A header is not required, but the order matters.

for csv : lab accession, barcode, alternate id, illumina run

OR

for tsv : lab accession\tbarcode\talternate id\tillumina run

USAGE
combine_rename.sh \\
\t-f <sample key> \\
\t-o <directory for combined and renamed files> \\
\t-i <directory for Illumina reads to use in polishing> \\
\t-w <directory where Illumina reads are> \\
\t-p <path to fastq_pass directory> \\
\t-j <number of jobs for parallel>
"""

combined="combined"
illumina="illumina"
wgs="/Volumes/IDGenomics_NAS/WGS_Serotyping"
sample_key="samples.tsv"
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
    i)
      illumina=$OPTARG
    ;;
    w)
      wgs=$OPTARG
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

mkdir -p $illumina
if [ -d "$illumina" ]
then
  echo "$(date): The corresponding paired-end Illumina files will be copied to $illumina"
else
  echo "$(date): FATAL : output directory $illumina could not be created. Set with -i <outdir>"
  exit 1
fi

if [ -d "$wgs" ]
then
  echo "$(date): The corresponding directory to look for paired-end Illumina files is $wgs"
else
  echo "$(date): FATAL : output directory $wgs was not found. Set with -w <dir>"
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

if [ -z "$(which bgzip)" ] ; then echo "$(date) : FATAL : bgzip is not in PATH" ; exit 1 ; fi
if [ -z "$(which parallel)" ] ; then echo "$(date) : FATAL : parallel is not in PATH" ; exit 1 ; fi

# sample key columns : lab accession, barcode, alternate id, illumina run
echo "$(date): Compressing any fastq files with bgzip"
ls $pass/*/fastq 2> /dev/null | parallel $jobs bgzip {}

echo "$(date): Concenating fastq.gz files and renaming"
cat $sample_key | tr "," "\t" | parallel --colsep "\t" $jobs "cat $pass/{2}/*fastq.gz > $combined/{1}.fastq.gz"

echo "$(date): Copying over Illumina reads"
echo "$(date): Copying R1"
cat $sample_key | tr "," "\t" | parallel --colsep "\t" $jobs "find $wgs/{4} -name \"{1}*R1*fastq.gz\" -o -name \"{3}*R1*fastq.gz\" | head -n 1 | parallel cp \{\} $illumina/{1}_R1.fastq.gz"
echo "$(date): Copying R2"
cat $sample_key | tr "," "\t" | parallel --colsep "\t" $jobs "find $wgs/{4} -name \"{1}*R2*fastq.gz\" -o -name \"{3}*R2*fastq.gz\" | head -n 1 | parallel cp \{\} $illumina/{1}_R2.fastq.gz"

echo "$(date): Everything is ready for Donut Falls!"
echo "$(date): Run with \"nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads combined --illumina illumina\""

exit 0
