#!/usr/bin/env nextflow
/*
 * GDC BAM -> FASTQ -> Salmon quantification, with cleanup of large intermediates.
 *
 * Steps per sample:
 *   1. DOWNLOAD_BAM   gdc-client download of the (genomic, STAR) BAM using a GDC token
 *   2. BAM_TO_FASTQ   samtools collate | samtools fastq  (unmapped reads included)
 *   3. SALMON_QUANT   salmon quant (lanes / multiple BAMs of one sample are combined)
 *   4. CLEANUP        delete the BAM and FASTQ files, only after Salmon succeeded
 *
 * Samplesheet (CSV, header required):
 *   sample,gdc_file_id[,case,tissue_type,...]
 * Create it with:  python3 fetch_gdc_metadata.py --project TCGA-LIHC -o samplesheet.csv
 *            or:    python3 fetch_gdc_metadata.py ids_or_old_samplesheet.csv -o samplesheet.csv
 * (adds a unique sample name and the tissue type, Tumor/Normal, from the GDC API).
 * Only the sample and gdc_file_id columns are used by the pipeline; extra columns are ignored.
 * Sample names must be unique per BAM (rows sharing a name are merged only with --merge_same_sample).
 *
 * Salmon-only mode: if you already have FASTQs, pass --fastq_samplesheet instead
 * (CSV with columns sample,fastq_1,fastq_2). Download, BAM->FASTQ and cleanup are skipped.
 *
 * SRA mode: pass --sra_samplesheet (CSV with columns sample,srr_run) to fetch the
 * reads from SRA instead of the GDC -- used for the CCLE/DepMap cell lines, whose
 * RNA-seq lives in BioProject PRJNA523380. Steps per run:
 *   1. SRA_TO_FASTQ   prefetch + fasterq-dump --split-3 (no token needed; public data)
 *   2. SALMON_QUANT   as for the GDC path
 * Build the samplesheet with:
 *   python3 fetch_sra_metadata.py --bioproject PRJNA523380 \
 *       --model-csv Model.csv --lineage Liver --paired-only -o sra_samplesheet.csv
 * Only `sample` and `srr_run` are used; the DepMap columns are carried for reference.
 * As with the GDC path, --fused true runs download + dump + Salmon in ONE task on
 * local disk (QUANT_FROM_SRA), which is what you want on Google Batch.
 *
 * All parameters and their defaults are defined in nextflow.config.
 */

nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    Usage -- pick exactly ONE input mode:

      GDC BAMs (TCGA):
        nextflow run main.nf -profile docker \\
          --samplesheet samplesheet.csv \\
          --gdc_token /path/to/gdc-user-token.txt \\
          --salmon_index /path/to/salmon_index_v49 \\
          --gene_map gencode.v49.annotation.gtf.gz \\
          --outdir results

      SRA runs (CCLE/DepMap cell lines; no token needed):
        nextflow run main.nf -profile docker \\
          --sra_samplesheet sra_samplesheet.csv \\
          --salmon_index /path/to/salmon_index_v49 \\
          --gene_map gencode.v49.annotation.gtf.gz \\
          --outdir results_ccle

      Existing FASTQs:
        nextflow run main.nf --fastq_samplesheet fastqs.csv \\
          --salmon_index ... --gene_map ... --outdir results

    Add --fused true (with --samplesheet or --sra_samplesheet) to do
    download + conversion + Salmon in one task on local disk; use it on
    Google Batch. Batching: --skip_rows N --limit_rows M.
    """.stripIndent()
}

process DOWNLOAD_BAM {
    tag "${sample}:${file_id}"
    label 'download'
    maxForks params.max_downloads

    input:
    tuple val(sample), val(file_id)
    path token

    output:
    tuple val(sample), val(file_id), path("${file_id}/*.bam"), emit: bam

    script:
    """
    gdc-client download \\
        --token-file ${token} \\
        --dir . \\
        --n-processes ${task.cpus} \\
        ${file_id}
    """
}

process BAM_TO_FASTQ {
    tag "${sample}:${file_id}"
    label 'samtools'
    publishDir "${params.outdir}/fastq_stats", pattern: '*.fastq_stats.txt', mode: 'copy'

    input:
    tuple val(sample), val(file_id), path(bam)

    output:
    tuple val(sample), path("${file_id}_R1.fastq.gz"), path("${file_id}_R2.fastq.gz"), emit: fastq
    path "${file_id}.fastq_stats.txt", emit: stats

    script:
    """
    # Name-collate so mates are adjacent, then write paired FASTQ.
    # samtools fastq excludes secondary/supplementary alignments by default (-F 0x900),
    # so every read is written once; unmapped reads are kept.
    samtools collate -@ ${task.cpus} -u -O ${bam} ${file_id}.collate_tmp | \\
        samtools fastq -@ ${task.cpus} -n -c 4 \\
            -1 ${file_id}_R1.fastq.gz \\
            -2 ${file_id}_R2.fastq.gz \\
            -0 /dev/null \\
            -s ${file_id}_singletons.fastq.gz \\
            -

    # An (almost) empty R1 means the BAM was single-end or pairing is broken
    R1_BYTES=\$(stat -c%s ${file_id}_R1.fastq.gz)
    if [ "\$R1_BYTES" -lt 100 ]; then
        echo "ERROR: R1 FASTQ is empty; BAM may be single-end or lack paired reads" >&2
        exit 1
    fi

    # Many singletons would indicate broken pairing in the BAM
    N_SINGLE=\$(( \$(gzip -dc ${file_id}_singletons.fastq.gz | wc -l) / 4 ))
    printf 'file_id\\tsingleton_reads\\n%s\\t%s\\n' "${file_id}" "\$N_SINGLE" > ${file_id}.fastq_stats.txt
    rm -f ${file_id}_singletons.fastq.gz
    """
}

/*
 * SRA run -> paired FASTQ. Public data, so no token is involved.
 *
 * sra-tools refuses to start without a settings file and, left alone, writes one
 * into $HOME -- which is read-only or shared in most containers. We write a
 * minimal one into the task directory instead. `report_instance_identity` lets
 * prefetch recognise it is running inside GCP and read SRA's cloud mirror
 * directly, which is both faster and avoids egress on Google Batch.
 */
process SRA_TO_FASTQ {
    tag "${sample}:${srr}"
    label 'sra'
    maxForks params.max_sra_downloads
    publishDir "${params.outdir}/fastq_stats", pattern: '*.fastq_stats.txt', mode: 'copy'

    input:
    tuple val(sample), val(srr)

    output:
    tuple val(sample), path("${srr}_R1.fastq.gz"), path("${srr}_R2.fastq.gz"), emit: fastq
    path "${srr}.fastq_stats.txt", emit: stats

    script:
    """
    export NCBI_SETTINGS="\$PWD/user-settings.mkfg"
    cat > "\$NCBI_SETTINGS" <<'MKFG'
/LIBS/GUID = "5c0d1f66-0a0b-4e1a-9a5a-2d1f6b5a7c31"
/libs/cloud/report_instance_identity = "true"
MKFG

    # --max-size u: prefetch's default 20 GB cap SKIPS larger runs with only a
    # warning, which would silently drop samples from the cohort.
    prefetch --max-size u --output-directory . ${srr}

    # --split-3: R1/R2 for properly paired spots, plus <srr>.fastq for unpaired.
    fasterq-dump --split-3 --threads ${task.cpus} \\
        --temp ./fasterq_tmp --outdir . ${srr}/${srr}.sra

    if [ ! -s ${srr}_1.fastq ] || [ ! -s ${srr}_2.fastq ]; then
        echo "ERROR: ${srr} produced no paired FASTQ; it is probably single-end." >&2
        echo "       Re-run fetch_sra_metadata.py with --paired-only to exclude it." >&2
        exit 1
    fi

    # Unpaired spots are the SRA analogue of BAM_TO_FASTQ's singletons.
    N_SINGLE=0
    if [ -s ${srr}.fastq ]; then
        N_SINGLE=\$(( \$(wc -l < ${srr}.fastq) / 4 ))
    fi
    printf 'srr_run\\tsingleton_reads\\n%s\\t%s\\n' "${srr}" "\$N_SINGLE" > ${srr}.fastq_stats.txt

    if command -v pigz >/dev/null 2>&1; then
        pigz -p ${task.cpus} -c ${srr}_1.fastq > ${srr}_R1.fastq.gz
        pigz -p ${task.cpus} -c ${srr}_2.fastq > ${srr}_R2.fastq.gz
    else
        gzip -c ${srr}_1.fastq > ${srr}_R1.fastq.gz
        gzip -c ${srr}_2.fastq > ${srr}_R2.fastq.gz
    fi

    # Drop the .sra and the uncompressed FASTQ now; only the .gz are outputs.
    rm -rf ${srr} ./fasterq_tmp ${srr}_1.fastq ${srr}_2.fastq ${srr}.fastq
    """
}

/*
 * Salmon indexes are directories. Nextflow can stage them directly,
 * but a tar.gz is also accepted.
 */
process UNTAR_INDEX {
    label 'light'

    input:
    path archive

    output:
    path 'salmon_index'

    script:
    """
    mkdir salmon_index
    tar xzf ${archive} --directory salmon_index
    """
}

process SALMON_QUANT {
    tag "${sample}"
    label 'salmon'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    // stageAs avoids name collisions between lanes / BAMs
    tuple val(sample), path(fastq_1, stageAs: 'r1_*/*'), path(fastq_2, stageAs: 'r2_*/*')
    path index
    path gene_map

    output:
    tuple val(sample), path("${sample}"), emit: quant_dir

    script:
    """
    # Unzip the gene map if needed (salmon --geneMap wants plain text)
    GENE_MAP="${gene_map}"
    if [[ "\$GENE_MAP" == *.gz ]]; then
        gunzip -c "\$GENE_MAP" > "${gene_map.baseName}"
        GENE_MAP="${gene_map.baseName}"
    fi

    salmon quant \\
        --no-version-check \\
        --index ${index} \\
        --libType A \\
        --mates1 ${fastq_1} \\
        --mates2 ${fastq_2} \\
        --threads ${task.cpus} \\
        --numBootstraps ${params.num_bootstraps} \\
        --seqBias \\
        --gcBias \\
        --dumpEq \\
        --geneMap "\$GENE_MAP" \\
        --output ${sample}
    """
}

/*
 * Delete the large intermediates of a sample. This process only runs once
 * SALMON_QUANT has finished successfully for that sample (it needs its output),
 * and it re-checks that quant.sf exists before deleting anything.
 *
 * Inputs are staged as symlinks, so we resolve them and delete the real files
 * in the work directories of the upstream tasks. This assumes a local or shared
 * filesystem (not an executor that copies inputs to remote storage).
 */
process CLEANUP {
    tag "${sample}"
    label 'light'

    input:
    tuple val(sample), path(bams, stageAs: 'bam_*/*'), path(fastq_1, stageAs: 'r1_*/*'), path(fastq_2, stageAs: 'r2_*/*'), path(quant_dir)

    output:
    val sample, emit: done

    script:
    """
    if [ ! -s ${quant_dir}/quant.sf ]; then
        echo "ERROR: ${quant_dir}/quant.sf missing or empty; refusing to delete inputs" >&2
        exit 1
    fi

    # Resolve every symlink to its real file BEFORE deleting anything (a dangling link cannot be resolved).
    for f in ${bams}; do
        target=\$(readlink -f "\$f")
        echo "Removing \$target"
        # also the BAM index files (.bam.bai / .bai / .csi) that gdc-client downloads next to the BAM; they are not staged
        rm -fv "\$target" "\${target}.bai" "\${target%.bam}.bai" "\${target}.csi"
    done

    for f in ${fastq_1} ${fastq_2}; do
        target=\$(readlink -f "\$f")
        echo "Removing \$target"
        rm -f "\$target"
    done
    """
}

/*
 * Fused mode (Google Batch): one task per sample does download -> BAM->FASTQ -> Salmon entirely on the
 * task's local disk, so no BAM/FASTQ is ever written to the bucket or moved between tasks, and the
 * large files are deleted as soon as they are no longer needed (BAM after the FASTQs exist, FASTQs
 * after Salmon succeeded). Only the Salmon output is published. One BAM = one sample here.
 */
process QUANT_FROM_GDC {
    tag "${sample}:${file_id}"
    label 'fused'
    publishDir "${params.outdir}", mode: 'copy', pattern: "${sample}"
    publishDir "${params.outdir}/fastq_stats", mode: 'copy', pattern: '*.fastq_stats.txt'

    input:
    tuple val(sample), val(file_id)
    path token
    path index
    path gene_map

    output:
    tuple val(sample), path("${sample}"), emit: quant_dir
    path "${file_id}.fastq_stats.txt", emit: stats

    script:
    """
    # 1. download the BAM
    gdc-client download \\
        --token-file ${token} \\
        --dir . \\
        --n-processes 4 \\
        ${file_id}
    BAM=\$(ls ${file_id}/*.bam)

    # 2. BAM -> paired FASTQ (same commands as BAM_TO_FASTQ), then drop the BAM and its index
    samtools collate -@ 4 -u -O "\$BAM" ${file_id}.collate_tmp | \\
        samtools fastq -@ 4 -n -c 4 \\
            -1 ${file_id}_R1.fastq.gz \\
            -2 ${file_id}_R2.fastq.gz \\
            -0 /dev/null \\
            -s ${file_id}_singletons.fastq.gz \\
            -
    rm -rf ${file_id} ${file_id}.collate_tmp*

    R1_BYTES=\$(stat -c%s ${file_id}_R1.fastq.gz)
    if [ "\$R1_BYTES" -lt 100 ]; then
        echo "ERROR: R1 FASTQ is empty; BAM may be single-end or lack paired reads" >&2
        exit 1
    fi
    N_SINGLE=\$(( \$(gzip -dc ${file_id}_singletons.fastq.gz | wc -l) / 4 ))
    printf 'file_id\\tsingleton_reads\\n%s\\t%s\\n' "${file_id}" "\$N_SINGLE" > ${file_id}.fastq_stats.txt
    rm -f ${file_id}_singletons.fastq.gz

    # 3. Salmon
    GENE_MAP="${gene_map}"
    if [[ "\$GENE_MAP" == *.gz ]]; then
        gunzip -c "\$GENE_MAP" > "${gene_map.baseName}"
        GENE_MAP="${gene_map.baseName}"
    fi

    salmon quant \\
        --no-version-check \\
        --index ${index} \\
        --libType A \\
        --mates1 ${file_id}_R1.fastq.gz \\
        --mates2 ${file_id}_R2.fastq.gz \\
        --threads ${task.cpus} \\
        --numBootstraps ${params.num_bootstraps} \\
        --seqBias \\
        --gcBias \\
        --dumpEq \\
        --geneMap "\$GENE_MAP" \\
        --output ${sample}

    # 4. only now (Salmon succeeded) remove the FASTQs, after checking the output looks complete
    if [ ! -s ${sample}/quant.sf ] || [ "\$(wc -l < ${sample}/quant.genes.sf)" -lt 1000 ]; then
        echo "ERROR: ${sample}/quant.sf or quant.genes.sf missing/too small" >&2
        exit 1
    fi
    rm -f ${file_id}_R1.fastq.gz ${file_id}_R2.fastq.gz
    """
}

/*
 * Fused SRA mode (Google Batch): the SRA analogue of QUANT_FROM_GDC. One task per
 * sample does prefetch -> fasterq-dump -> Salmon entirely on the task's local disk,
 * so no .sra or FASTQ is ever written to the bucket, and each is deleted as soon as
 * it is no longer needed. Only the Salmon output is published. One run = one sample
 * here (build the samplesheet with --one-run-per-line).
 */
process QUANT_FROM_SRA {
    tag "${sample}:${srr}"
    label 'sra_fused'
    publishDir "${params.outdir}", mode: 'copy', pattern: "${sample}"
    publishDir "${params.outdir}/fastq_stats", mode: 'copy', pattern: '*.fastq_stats.txt'

    input:
    tuple val(sample), val(srr)
    path index
    path gene_map

    output:
    tuple val(sample), path("${sample}"), emit: quant_dir
    path "${srr}.fastq_stats.txt", emit: stats

    script:
    """
    # 1. settings file in the task dir (see SRA_TO_FASTQ for why)
    export NCBI_SETTINGS="\$PWD/user-settings.mkfg"
    cat > "\$NCBI_SETTINGS" <<'MKFG'
/LIBS/GUID = "5c0d1f66-0a0b-4e1a-9a5a-2d1f6b5a7c31"
/libs/cloud/report_instance_identity = "true"
MKFG

    # 2. fetch and convert, then drop the .sra and the temp dir
    prefetch --max-size u --output-directory . ${srr}
    fasterq-dump --split-3 --threads ${task.cpus} \\
        --temp ./fasterq_tmp --outdir . ${srr}/${srr}.sra
    rm -rf ${srr} ./fasterq_tmp

    if [ ! -s ${srr}_1.fastq ] || [ ! -s ${srr}_2.fastq ]; then
        echo "ERROR: ${srr} produced no paired FASTQ; it is probably single-end." >&2
        exit 1
    fi
    N_SINGLE=0
    if [ -s ${srr}.fastq ]; then
        N_SINGLE=\$(( \$(wc -l < ${srr}.fastq) / 4 ))
    fi
    printf 'srr_run\\tsingleton_reads\\n%s\\t%s\\n' "${srr}" "\$N_SINGLE" > ${srr}.fastq_stats.txt
    rm -f ${srr}.fastq

    # 3. Salmon reads the uncompressed FASTQ directly -- no point gzipping a file
    #    that is deleted in the next step.
    GENE_MAP="${gene_map}"
    if [[ "\$GENE_MAP" == *.gz ]]; then
        gunzip -c "\$GENE_MAP" > "${gene_map.baseName}"
        GENE_MAP="${gene_map.baseName}"
    fi

    salmon quant \\
        --no-version-check \\
        --index ${index} \\
        --libType A \\
        --mates1 ${srr}_1.fastq \\
        --mates2 ${srr}_2.fastq \\
        --threads ${task.cpus} \\
        --numBootstraps ${params.num_bootstraps} \\
        --seqBias \\
        --gcBias \\
        --dumpEq \\
        --geneMap "\$GENE_MAP" \\
        --output ${sample}

    # 4. only now (Salmon succeeded) remove the FASTQs
    if [ ! -s ${sample}/quant.sf ] || [ "\$(wc -l < ${sample}/quant.genes.sf)" -lt 1000 ]; then
        echo "ERROR: ${sample}/quant.sf or quant.genes.sf missing/too small" >&2
        exit 1
    fi
    rm -f ${srr}_1.fastq ${srr}_2.fastq
    """
}

workflow {
    if (!params.salmon_index || !params.gene_map) {
        helpMessage()
        error "Missing required parameter(s): --salmon_index, --gene_map"
    }
    def chosen = ['--samplesheet': params.samplesheet,
                  '--fastq_samplesheet': params.fastq_samplesheet,
                  '--sra_samplesheet': params.sra_samplesheet].findAll { k, v -> v }.keySet()
    if (chosen.size() > 1) {
        error "Use exactly one input mode, got: ${chosen.join(', ')}"
    }
    if (chosen.isEmpty()) {
        helpMessage()
        error "Provide --samplesheet (GDC), --sra_samplesheet (SRA) or --fastq_samplesheet (existing FASTQs)"
    }

    if (params.sra_samplesheet) {
        // SRA mode: CSV columns sample,srr_run (extra columns ignored).
        // Build it with fetch_sra_metadata.py; see the header of this file.
        def sra_rows  = file(params.sra_samplesheet, checkIfExists: true).splitCsv(header: true)
        def sra_n_all = sra_rows.size()
        sra_rows = sra_rows.drop(params.skip_rows as int)
        if (params.limit_rows) { sra_rows = sra_rows.take(params.limit_rows as int) }
        log.info "SRA samplesheet: using ${sra_rows.size()} of ${sra_n_all} rows"

        def missing = sra_rows.findAll { !it.srr_run }
        if (missing) {
            error "${missing.size()} row(s) have no srr_run column. Expected CSV columns: sample,srr_run"
        }

        def sra_counts = sra_rows.countBy { it.sample }
        // Several runs of one cell line are normally lanes of one library, so
        // merging is usually right here -- but it is still opt-in, as for BAMs.
        def sra_dups = sra_counts.findAll { k, v -> v > 1 }.keySet()
        if (sra_dups && !params.merge_same_sample) {
            error "Sample name(s) used for more than one SRA run: ${sra_dups.join(', ')}. " +
                  "Pass --merge_same_sample true to merge them into one Salmon run, or re-run " +
                  "fetch_sra_metadata.py with --one-run-per-line."
        }

        sra_ch = Channel.fromList(sra_rows).map { row -> tuple(row.sample as String, row.srr_run as String) }

        if (params.fused) {
            if (sra_dups) { error "--fused needs one run per sample; duplicate sample name(s): ${sra_dups.join(', ')}" }
        }
        else {
            SRA_TO_FASTQ(sra_ch)
            reads_ch = SRA_TO_FASTQ.out.fastq
                .map { s, r1, r2 -> tuple(groupKey(s, sra_counts[s]), r1, r2) }
                .groupTuple()
                .map { key, r1, r2 -> tuple(key.getGroupTarget() as String, r1, r2) }
        }
    }
    else if (params.fastq_samplesheet) {
        // Salmon-only mode: skips download, BAM->FASTQ and cleanup.
        // CSV columns: sample,fastq_1,fastq_2 (use absolute paths; repeat a sample for several pairs)
        def fq_rows   = file(params.fastq_samplesheet, checkIfExists: true).splitCsv(header: true)
        def fq_n_all  = fq_rows.size()
        fq_rows = fq_rows.drop(params.skip_rows as int)
        if (params.limit_rows) { fq_rows = fq_rows.take(params.limit_rows as int) }
        log.info "FASTQ samplesheet: using ${fq_rows.size()} of ${fq_n_all} rows"
        def fq_counts = fq_rows.countBy { it.sample }

        def fq_dups = fq_counts.findAll { k, v -> v > 1 }.keySet()
        if (fq_dups && !params.merge_same_sample) {
            error "Sample name(s) used for more than one FASTQ pair: ${fq_dups.join(', ')}. " +
                  "Use unique sample names, or pass --merge_same_sample true to merge lanes of one library."
        }

        reads_ch = Channel.fromList(fq_rows)
            .map { row -> tuple(groupKey(row.sample as String, fq_counts[row.sample]),
                                file(row.fastq_1, checkIfExists: true),
                                file(row.fastq_2, checkIfExists: true)) }
            .groupTuple()
            .map { key, r1, r2 -> tuple(key.getGroupTarget() as String, r1, r2) }
    }
    else {
        if (!params.gdc_token) {
            helpMessage()
            error "Missing required parameter: --gdc_token"
        }

        // Read the samplesheet once to learn how many BAMs each sample has.
        // groupKey(sample, n) lets groupTuple emit a sample as soon as its n items
        // are ready, instead of waiting for every sample to finish.
        def rows   = file(params.samplesheet, checkIfExists: true).splitCsv(header: true)
        // Optionally use only part of the samplesheet (data rows, header not counted)
        def n_all  = rows.size()
        rows = rows.drop(params.skip_rows as int)
        if (params.limit_rows) { rows = rows.take(params.limit_rows as int) }
        log.info "Samplesheet: using ${rows.size()} of ${n_all} rows"
        def counts = rows.countBy { it.sample }

        // Rows sharing a sample name are merged into ONE Salmon run. That is right for lanes of
        // the same library, but wrong for e.g. tumor + normal BAMs of the same case.
        def dups = counts.findAll { k, v -> v > 1 }.keySet()
        if (dups && !params.merge_same_sample) {
            error "Sample name(s) used for more than one BAM: ${dups.join(', ')}. " +
                  "Generate the samplesheet with 'python3 fetch_gdc_metadata.py' (unique sample name + tissue type per BAM), " +
                  "or pass --merge_same_sample true if they really are the same library."
        }

        samples_ch = Channel.fromList(rows).map { row -> tuple(row.sample as String, row.gdc_file_id as String) }
        token      = file(params.gdc_token, checkIfExists: true)

        if (params.fused) {
            if (dups) { error "--fused needs one BAM per sample; duplicate sample name(s): ${dups.join(', ')}" }
        }
        else {
            DOWNLOAD_BAM(samples_ch, token)
            BAM_TO_FASTQ(DOWNLOAD_BAM.out.bam)

            reads_ch = BAM_TO_FASTQ.out.fastq
                .map { s, r1, r2 -> tuple(groupKey(s, counts[s]), r1, r2) }
                .groupTuple()
                .map { key, r1, r2 -> tuple(key.getGroupTarget() as String, r1, r2) }

            bams_ch = DOWNLOAD_BAM.out.bam
                .map { s, id, bam -> tuple(groupKey(s, counts[s]), bam) }
                .groupTuple()
                .map { key, bams -> tuple(key.getGroupTarget() as String, bams) }
        }
    }

    // Salmon index: directory or tar.gz
    index_path = file(params.salmon_index, checkIfExists: true)
    if (params.salmon_index.toString().endsWith('.tar.gz') || params.salmon_index.toString().endsWith('.tgz')) {
        index_ch = UNTAR_INDEX(index_path)
    } else {
        index_ch = Channel.value(index_path)
    }
    gene_map_ch = Channel.value(file(params.gene_map, checkIfExists: true))

    if (params.fused && params.sra_samplesheet) {
        QUANT_FROM_SRA(sra_ch, index_ch, gene_map_ch)
    }
    else if (params.fused && params.samplesheet) {
        QUANT_FROM_GDC(samples_ch, token, index_ch, gene_map_ch)
    }
    else {
        SALMON_QUANT(reads_ch, index_ch, gene_map_ch)

        // CLEANUP deletes BAMs + FASTQs from work/, so it only applies to the GDC
        // path. SRA_TO_FASTQ already removes the .sra and uncompressed FASTQ inside
        // its own task, and in --fastq_samplesheet mode the FASTQs are user inputs.
        if (params.samplesheet && params.cleanup) {
            CLEANUP( bams_ch.join(reads_ch).join(SALMON_QUANT.out.quant_dir) )
        }
    }
}
