# One image with everything the fused task (QUANT_FROM_GDC) needs: gdc-client + samtools + salmon.
# Versions mirror the conda env used on the VM (see README). gdc-client breaks on Python 3.14 -> keep 3.11.
#
# Build in the same region as the Batch jobs, e.g. with Cloud Build:
#   gcloud artifacts repositories create salmon-workflow --repository-format=docker --location=us-east1 --project=amaroq-495623
#   gcloud builds submit docker/ --tag us-east1-docker.pkg.dev/amaroq-495623/salmon-workflow/gdc-salmon:1.0 --project=amaroq-495623
FROM mambaorg/micromamba:1.5.10
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 gdc-client samtools=1.20 salmon=1.10.3 procps-ng \
    && micromamba clean -a -y
# Nextflow needs `bash` and `ps` in the image; put the conda env first on PATH
ENV PATH=/opt/conda/bin:$PATH
