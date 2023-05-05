# nextflow run /home/eriny/sandbox/Donut_Falls -profile singularity --reads  /home/eriny/sandbox/test_files/donut/combined -with-tower -resume 
# nextflow run /home/eriny/sandbox/Donut_Falls -profile singularity --sample_sheet /home/eriny/sandbox/test_files/donut/sample_sheet.csv -resume 

echo "$(date): testing only LR with flye" && \
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads  /home/eriny/sandbox/test_files/donut/combined \
    --outdir only_nanopore \
    -with-tower \
    -resume && \
echo "$(date): testing sample sheet" && \
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --sample_sheet /home/eriny/sandbox/test_files/donut/sample_sheet.csv \
    --outdir sample_sheet \
    --sequencing_summary /home/eriny/sandbox/test_files/donut/sequencing_summary_FAS76150_35058c5c.txt \
    -with-tower \
    -resume && \
echo "$(date): testing raven" && \
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --sample_sheet /home/eriny/sandbox/test_files/donut/sample_sheet.csv \
    --assembler raven \
    --outdir    raven \
    -with-tower \
    -resume && \
echo "$(date): testing miniasm" && \
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --sample_sheet /home/eriny/sandbox/test_files/donut/sample_sheet.csv \
    --assembler miniasm \
    --outdir    miniasm \
    -with-tower \
    -resume && \
echo "$(date): testing lr_unicycler" && \
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --sample_sheet /home/eriny/sandbox/test_files/donut/sample_sheet.csv \
    --assembler lr_unicycler \
    --outdir    lr_unicycler \
    -with-tower \
    -resume && \
echo "$(date): testing unicycler" && \
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile    singularity \
    --sample_sheet /home/eriny/sandbox/test_files/donut/sample_sheet.csv \
    --assembler unicycler \
    --outdir    unicycler \
    -with-tower \
    -resume && \
echo "$(date): testing masurca" && \
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile    singularity \
    --sample_sheet /home/eriny/sandbox/test_files/donut/sample_sheet.csv \
    --assembler masurca \
    --outdir    masurca \
    -with-tower \
    -resume && \
echo "$(date): testing empty" && \
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --reads     shouldntexist \
    --illumina  wontexist \
    --outdir    nonexistent \
    --sequencing_summary doesntexit \
    -with-tower \
    -resume


echo "$(date): testing trycycler" && \
nextflow run /home/eriny/sandbox/Donut_Falls \
    -profile singularity \
    --sample_sheet /home/eriny/sandbox/test_files/donut/sample_sheet.csv \
    --outdir    trycycler \
    --assembler trycycler \
    --trycycler_min_fasta 12 \
    -with-tower \
    -resume