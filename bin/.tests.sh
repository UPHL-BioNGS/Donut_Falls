nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads /home/eriny/sandbox/test_files/donut/combined \
    --outdir only_nanopore \
    -with-tower \
    -resume

nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads /home/eriny/sandbox/test_files/donut/combined \
    --illumina /home/eriny/sandbox/test_files/donut/combined/illumina \
    --outdir default \
    -with-tower \
    -resume

nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads /home/eriny/sandbox/test_files/donut/combined \
    --illumina /home/eriny/sandbox/test_files/donut/combined/illumina \
    --assembler raven \
    --outdir raven \
    -with-tower \
    -resume

nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads /home/eriny/sandbox/test_files/donut/combined \
    --illumina /home/eriny/sandbox/test_files/donut/combined/illumina \
    --assembler miniasm \
    --outdir miniasm \
    -with-tower \
    -resume

nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads /home/eriny/sandbox/test_files/donut/combined \
    --illumina /home/eriny/sandbox/test_files/donut/combined/illumina \
    --assembler unicycler \
    --outdir unicycler \
    -with-tower \
    -resume
