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
 * All parameters and their defaults are defined in nextflow.config.
 */

nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    Usage:
      nextflow run main.nf -profile docker \\
        --samplesheet samplesheet.csv \\
        --gdc_token /path/to/gdc-user-token.txt \\
        --salmon_index /path/to/salmon_index_v49 \\
        --gene_map gencode.v49.annotation.gtf.gz \\
        --outdir results
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

    for f in ${bams} ${fastq_1} ${fastq_2}; do
        target=\$(readlink -f "\$f")
        echo "Removing \$target"
        rm -f "\$target"
    done

    # BAM index files (.bam.bai / .bai / .csi) are not staged, so remove them next to each real BAM
    for f in ${bams}; do
        target=\$(readlink -f "\$f")
        rm -fv "\${target}.bai" "\${target%.bam}.bai" "\${target}.csi"
    done
    """
}

workflow {
    if (!params.salmon_index || !params.gene_map) {
        helpMessage()
        error "Missing required parameter(s): --salmon_index, --gene_map"
    }
    if (params.samplesheet && params.fastq_samplesheet) {
        error "Use either --samplesheet (download from GDC) or --fastq_samplesheet (existing FASTQs), not both"
    }
    if (!params.samplesheet && !params.fastq_samplesheet) {
        helpMessage()
        error "Provide --samplesheet (GDC download mode) or --fastq_samplesheet (Salmon-only mode)"
    }

    if (params.fastq_samplesheet) {
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

    // Salmon index: directory or tar.gz
    index_path = file(params.salmon_index, checkIfExists: true)
    if (params.salmon_index.toString().endsWith('.tar.gz') || params.salmon_index.toString().endsWith('.tgz')) {
        index_ch = UNTAR_INDEX(index_path)
    } else {
        index_ch = Channel.value(index_path)
    }
    gene_map_ch = Channel.value(file(params.gene_map, checkIfExists: true))

    SALMON_QUANT(reads_ch, index_ch, gene_map_ch)

    if (!params.fastq_samplesheet && params.cleanup) {
        CLEANUP( bams_ch.join(reads_ch).join(SALMON_QUANT.out.quant_dir) )
    }
}
