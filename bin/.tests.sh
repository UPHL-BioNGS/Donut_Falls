# nextflow run /home/eriny/sandbox/Donut_Falls -profile singularity --reads /home/eriny/sandbox/test_files/donut/combined --illumina /home/eriny/sandbox/test_files/donut/illumina -with-tower -resume

echo "$(date): testing only LR with flye"
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads     /home/eriny/sandbox/test_files/donut/combined \
    --outdir    only_nanopore \
    -with-tower \
    -resume

echo "$(date): testing defaults"
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads     /home/eriny/sandbox/test_files/donut/combined \
    --illumina  /home/eriny/sandbox/test_files/donut/illumina \
    --outdir    default \
    --sequencing_summary /home/eriny/sandbox/test_files/donut/sequencing_summary_FAS76150_35058c5c.txt \
    -with-tower \
    -resume

echo "$(date): testing raven"
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads     /home/eriny/sandbox/test_files/donut/combined \
    --illumina  /home/eriny/sandbox/test_files/donut/illumina \
    --assembler raven \
    --outdir    raven \
    -with-tower \
    -resume

echo "$(date): testing miniasm"
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads     /home/eriny/sandbox/test_files/donut/combined \
    --illumina  /home/eriny/sandbox/test_files/donut/illumina \
    --assembler miniasm \
    --outdir    miniasm \
    -with-tower \
    -resume

echo "$(date): testing unicycler"
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile    singularity \
    --reads     /home/eriny/sandbox/test_files/donut/combined \
    --illumina  /home/eriny/sandbox/test_files/donut/illumina \
    --assembler unicycler \
    --outdir    unicycler \
    -with-tower \
    -resume


echo "$(date): testing empty"
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads     shouldntexist \
    --illumina  wontexist \
    --outdir    nonexistent \
    --sequencing_summary doesntexit \
    -with-tower \
    -resume