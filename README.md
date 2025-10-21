# Donut Falls

Named after the beautiful [Donut Falls](https://en.wikipedia.org/wiki/Doughnut_Falls)

Location: 40.630°N 111.655°W, Elevation: 7,942 ft (2,421 m), [Hiking level: easy](https://www.alltrails.com/trail/us/utah/cecret-lake-trail)

|                                       |                                       |
|:-------------------------------------:|:-------------------------------------:|
| <img src="https://images.alltrails.com/eyJidWNrZXQiOiJhc3NldHMuYWxsdHJhaWxzLmNvbSIsImtleSI6InVwbG9hZHMvcGhvdG8vaW1hZ2UvNTIzNDg5MTUvODRiOGEzM2M3MjliMzUyYjk4YTJhYmY5Mjg1MWMzMDguanBnIiwiZWRpdHMiOnsidG9Gb3JtYXQiOiJqcGVnIiwicmVzaXplIjp7IndpZHRoIjoyMDQ4LCJoZWlnaHQiOjIwNDgsImZpdCI6Imluc2lkZSJ9LCJyb3RhdGUiOm51bGwsImpwZWciOnsidHJlbGxpc1F1YW50aXNhdGlvbiI6dHJ1ZSwib3ZlcnNob290RGVyaW5naW5nIjp0cnVlLCJvcHRpbWlzZVNjYW5zIjp0cnVlLCJxdWFudGlzYXRpb25UYWJsZSI6M319fQ==" width="500" /> | <img src="https://images.alltrails.com/eyJidWNrZXQiOiJhc3NldHMuYWxsdHJhaWxzLmNvbSIsImtleSI6InVwbG9hZHMvcGhvdG8vaW1hZ2UvNTg3MzUyOTgvZWM4YTQ5NTZlNDhiNmZmMmU4ZWEzMzE0NjhhOWIyYWYuanBnIiwiZWRpdHMiOnsidG9Gb3JtYXQiOiJqcGVnIiwicmVzaXplIjp7IndpZHRoIjoyMDQ4LCJoZWlnaHQiOjIwNDgsImZpdCI6Imluc2lkZSJ9LCJyb3RhdGUiOm51bGwsImpwZWciOnsidHJlbGxpc1F1YW50aXNhdGlvbiI6dHJ1ZSwib3ZlcnNob290RGVyaW5naW5nIjp0cnVlLCJvcHRpbWlzZVNjYW5zIjp0cnVlLCJxdWFudGlzYXRpb25UYWJsZSI6M319fQ==" width="500"/> |


([Image credit: User submitted photos at alltrails.com](https://www.alltrails.com/trail/us/utah/donut-falls-trail/photos))

More information about the trail leading up to this landmark can be found at [utah.com/hiking/donut-falls](https://utah.com/hiking/donut-falls)

Donut Falls is a [Nextflow](https://www.nextflow.io/) workflow developed by [@erinyoung](https://github.com/erinyoung) at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for long-read [nanopore](https://nanoporetech.com) sequencing of microbial isolates. Built to work on linux-based operating systems. Additional config options are needed for cloud batch usage.

Donut Falls is also included in the staphb toolkit [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit). 

We made a [wiki](https://github.com/UPHL-BioNGS/Donut_Falls/wiki). Please read it!

## Wiki table of contents:
- [Installation](https://github.com/UPHL-BioNGS/Donut_Falls/wiki/Installation)
- [Usage](https://github.com/UPHL-BioNGS/Donut_Falls/wiki/Usage)
  - [Using config files](https://github.com/UPHL-BioNGS/Donut_Falls/wiki/Usage#using-a-config-file)
  - [Parameters worth adjusting](https://github.com/UPHL-BioNGS/Donut_Falls/wiki/Usage#recommended-parameters-to-adjust)
  - [Examples](https://github.com/UPHL-BioNGS/Donut_Falls/wiki/Usage#examples)
  - [Using as a subworkflow](https://github.com/UPHL-BioNGS/Donut_Falls/wiki/Linking)
- [Workflow DAG](https://github.com/UPHL-BioNGS/Donut_Falls/wiki#basic-diagram-of-the-workflow-and-subworkflows)
- [Assessing Assembly Quality](https://github.com/UPHL-BioNGS/Donut_Falls/wiki/evaluation)
- [FAQ](https://github.com/UPHL-BioNGS/Donut_Falls/wiki/FAQ)

## Getting started

### Install dependencies
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Singularity](https://singularity.lbl.gov/install-linux) or [Docker](https://docs.docker.com/get-docker/)

## Quick start

```
nextflow run UPHL-BioNGS/Donut_Falls -profile <singularity or docker> --sample_sheet <sample_sheet.csv>
```

## Sample Sheets

Sample sheet is a csv file with the name of the sample and corresponding nanopore fastq.gz file on a single row with header `sample` and `fastq`. When Illumina fastq files are available for polishing or hybrid assembly, they are added to end of each row under column header `fastq_1` and `fastq_2`. 

Option 1 : just nanopore reads
```
sample,fastq
test,long_reads_low_depth.fastq.gz
```

Option 2 : nanopore reads and at least one sample has Illumina paired-end fastq files
```
sample,fastq,fastq_1,fastq_2
sample1,sample1.fastq.gz,sample1_R1.fastq.gz,sample1_R2.fastq.gz
sample2,sample2.fastq.gz,,
```

### Switching assemblers
There are currently several options for assembly
- [flye](https://github.com/fenderglass/Flye) (default)
- [raven](https://github.com/lbcb-sci/raven)
- [myloasm](https://github.com/bluenote-1577/myloasm)
- [unicycler](https://github.com/rrwick/Unicycler) (requires short reads for hybrid assembly)

These are specified with the `assembler` paramater. If Illumina reads are found, then flye, myloasm, and raven assemblies will be polished with those reads.

Note: more than one assembler can be chosen (i.e. `params.assembler = 'flye,raven'`). This will run the input files on each assembler listed. Listing an assembler more than once will not create additional assemblies with that tool (i.e. `params.assembler = 'flye,flye,flye'` will still only run the input files through flye once).

## Examples
```
# nanopore assembly with flye followed by polishing if illumina files are supplied
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --sample_sheet sample_sheet.csv

# or with docker and specifying the assembler
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --sample_sheet sample_sheet.csv --assembler flye

# hybrid assembly with unicycler where both nanopore and illumina files are required
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --sample_sheet sample_sheet.csv --assembler unicycler

# assembling with all asssemblers
# specifying the results to be stored in 'donut_falls_test_results' instead of 'donut_falls'
# using docker instead of singularity
nextflow run UPHL-BioNGS/Donut_Falls -profile docker --sample_sheet sample_sheet.csv --assembler unicycler_flye_raven_myloasm


# using some test files (requires internet connection)
nextflow run UPHL-BioNGS/Donut_Falls -profile docker,test

```

## Acknowledgement

Donut Falls would not be possible without
- [bandage](https://github.com/rrwick/Bandage) : visualize gfa files
- [bcftools](https://github.com/samtools/bcftools) : apply clair3 polishing
- [busco](https://gitlab.com/ezlab/busco) : assessment of assembly quality
- [bwa](https://github.com/lh3/bwa) : aligning reads for polypolish
- [circulocov](https://github.com/erinyoung/CirculoCov) : read depth per contig, aligning reads for clair3
- [clair3](https://github.com/HKU-BAL/Clair3) : long-read polishing
- [dnaapler](https://github.com/gbouras13/dnaapler) : rotation
- [fastp](https://github.com/OpenGene/fastp) : cleaning illumina reads
- [fastplong](https://github.com/OpenGene/fastplong) : cleaning nanopore reads
- [flye](https://github.com/fenderglass/Flye) : de novo assembly (default assembler)
- [gfastats](https://github.com/vgl-hub/gfastats) : assessment of assembly
- [mash](https://github.com/marbl/Mash) : determines distance between reads for QC check
- [multiqc](https://multiqc.info/) : amalgamation of results
- [myloasm](https://github.com/bluenote-1577/myloasm) : de novo assembly (params.assembler = 'myloasm')
- [polypolish](https://github.com/rrwick/Polypolish) : reduces sequencing artefacts through polishing with Illumina reads
- [pypolca](https://github.com/gbouras13/pypolca) : reduces sequencing artefacts through polishing with Illumina reads
- [rasusa](https://github.com/mbhall88/rasusa) : subsampling nanopore reads to 150X depth
- [raven](https://github.com/lbcb-sci/raven) : de novo assembly option (params.assembler = 'raven')
- [seqkit](https://github.com/shenwei356/seqkit) : fastq stats and fasta sorting
- [unicycler](https://github.com/rrwick/Unicycler) : hybrid assembly option (params.assembler = 'unicycler')
