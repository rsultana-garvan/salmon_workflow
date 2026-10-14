# Create local conda environment
```
conda create -n nextflow -c conda-forge -c bioconda nextflow=24.10 openjdk=17 python=3.11 gdc-client samtools=1.20 salmon=1.10.3
```

# Download the samplesheet
```
python3 fetch_gdc_metadata.py --project TCGA-LIHC -o samplesheet.csv
```

# Run 
```
nextflow run main.nf \
  --samplesheet samplesheet.csv \
  --limit_rows 2 \
  --gdc_token /data/controlled_data/gdc-user-token.2026-10-14T07_53_35.198Z.txt \
  --salmon_index /data/annotation/salmon_index \
  --gene_map /data/annotation/gencode.v49.genemap.tsv.gz \
  --num_bootstraps 0 \
  --outdir results
```
