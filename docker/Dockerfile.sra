# One image with everything the fused SRA task (QUANT_FROM_SRA) needs:
# sra-tools (prefetch + fasterq-dump) + salmon. The SRA counterpart of
# ./Dockerfile, which pairs gdc-client + samtools + salmon for QUANT_FROM_GDC.
#
# pigz is included so SRA_TO_FASTQ compresses with all cores; the process falls
# back to gzip if it is absent, but that is a lot slower on a 20 GB FASTQ.
#
# salmon must match the version used to build the index (see nextflow.config).
#
# Build in the same region as the Batch jobs, e.g. with Cloud Build:
#   gcloud builds submit docker/ --config=- --project=amaroq-495623 <<'EOF'
#   steps:
#     - name: gcr.io/cloud-builders/docker
#       args: ['build','-f','Dockerfile.sra','-t',
#              'us-east1-docker.pkg.dev/amaroq-495623/salmon-workflow/sra-salmon:1.0','.']
#   images: ['us-east1-docker.pkg.dev/amaroq-495623/salmon-workflow/sra-salmon:1.0']
#   EOF
# then set params.sra_fused_container to that tag.
FROM mambaorg/micromamba:1.5.10
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        sra-tools=3.1.1 salmon=1.10.3 pigz procps-ng \
    && micromamba clean -a -y
# Nextflow needs `bash` and `ps` in the image; put the conda env first on PATH
ENV PATH=/opt/conda/bin:$PATH
