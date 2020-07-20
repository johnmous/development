# Title: Run STARsolo for the Li (http://dx.doi.org/10.1016/j.stem.2017.03.007) raw data. Rename the barcodes to match the cell IDs used by Li
# See slides in repo: https://github.com/zorrodong/HECA/blob/master/scRNA-seq_pipeline_hg38/Single-cell_RNA-seq_pipeline.pptx
# Also see this repo, shich contains the pipeline from authors: https://github.com/zorrodong/HECA/tree/master/scRNA-seq_pipeline_hg38
# Author: I. Moustakas

OUTPUT_PATH = "/path/to/project-248-FGC_RNAseq/analysis/STARsolo_masked_no_read_QC"
INPUT_PATH = "/path/to/ovogrowth-hpc/project-248-FGC_RNAseq/data/LiLiData/sra/late_meiot"
SRAs = ["SRR4081197", "SRR4081198", "SRR4081199", "SRR4081200", "SRR4081201",
        "SRR4199308", "SRR4199309", "SRR4199310", "SRR4199311", "SRR4199312",
        "SRR4199313", "SRR4199314", "SRR4199315", "SRR4199316", "SRR4199317",
        "SRR4199318", "SRR4199319", "SRR4199320"]

rule all:
    input: [OUTPUT_PATH + "/{sra}/Log.out".format(sra = sra ) for sra in SRAs]

rule STARsolo:
    input:
        fq_r1 = INPUT_PATH + "/{sra}_pass_1.fastq.gz",
        fq_r2 = INPUT_PATH + "/{sra}_pass_2.fastq.gz"
    output: OUTPUT_PATH + "/{sra}/Log.out"
    params:
        genome = "/path/to/ioannis/cellRangerPipelines/refdata-cellranger-GRCh38-3.0.0_STAR_2.7.3a_masked/star",
        BCWhitelist = "/path/to/project-248-FGC_RNAseq/src/HECA/scRNA-seq_pipeline_hg38/barcode_96_8bp_noCells.csv",
        outputdir = OUTPUT_PATH + "/{sra}/",
        sthreads = 10,
        mem = "40G",
        time= "18:0:0"
    conda:
        "/path/to/project-248-FGC_RNAseq/src/STARsoloPipe/STARsolo.yml"
    shell:
        """
        STAR \
        --soloType CB_UMI_Simple \
        --soloCBstart 1 \
        --soloCBlen 8 \
        --soloUMIstart 9 \
        --soloUMIlen 8 \
        --soloBarcodeReadLength 0 \
        --runThreadN {params.sthreads} \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM CR CY UR UY CB UB \
        --readFilesCommand zcat \
        --genomeDir {params.genome} \
        --soloCBwhitelist /path/to/project-248-FGC_RNAseq/src/HECA/scRNA-seq_pipeline_hg38/barcode_96_8bp_noCells.csv \
        --outFileNamePrefix {params.outputdir} \
        --readFilesIn {input.fq_r1} \
        {input.fq_r2}
        """

rule barcode2CellId:
    input: OUTPUT_PATH + "/{sra}/Log.out"
    output: OUTPUT_PATH + "/{sra}/Solo.out/Gene/filtered.cellid/barcodes.tsv.gz"
    params:
        srrid="{sra}",
        srr2sexweekfile="/path/to/project-248-FGC_RNAseq/data/LiLiData/SRRid2sexWeek.csv",
        barcd2cellidfile="/path/to/project-248-FGC_RNAseq/src/HECA/scRNA-seq_pipeline_hg38/barcode_96_8bp.txt",
        experimentbarcdfile=OUTPUT_PATH + "/{sra}/Solo.out/Gene/filtered/barcodes.tsv",
        experimentfolder=OUTPUT_PATH + "/{sra}/Solo.out/Gene/filtered/",
        experimentcellidfile=OUTPUT_PATH + "/{sra}/Solo.out/Gene/filtered.cellid/barcodes.tsv.gz",
        experimentcellidfolder=OUTPUT_PATH + "/{sra}/Solo.out/Gene/filtered.cellid/",
        qthreads = 1,
        h_vmem = "2G"
    conda:
        "/path/to/project-248-FGC_RNAseq/src/STARsoloPipe/STARsolo.yml"
    shell:
        """
        mkdir -p {params.experimentcellidfolder} &&
        cp {params.experimentfolder}/features.tsv  {params.experimentcellidfolder} &&
        cp {params.experimentfolder}/matrix.mtx  {params.experimentcellidfolder} &&
        gzip {params.experimentcellidfolder}* &&
        python barcode2CellId.py \
        --srrid {params.srrid} \
        --srr2sexweekfile {params.srr2sexweekfile} \
        --barcd2cellidfile {params.barcd2cellidfile} \
        --experimentbarcdfile {params.experimentbarcdfile} \
        --experimentcellidfile {params.experimentcellidfile} | gzip
        """

rule barcode2CellIdraw:
    input: OUTPUT_PATH + "/{sra}/Log.out"
    output: OUTPUT_PATH + "/{sra}/Solo.out/Gene/raw.cellid/barcodes.tsv.gz"
    params:
        srrid="{sra}",
        srr2sexweekfile="/path/to/project-248-FGC_RNAseq/data/LiLiData/SRRid2sexWeek.csv",
        barcd2cellidfile="/path/to/project-248-FGC_RNAseq/src/HECA/scRNA-seq_pipeline_hg38/barcode_96_8bp.txt",
        experimentbarcdfile=OUTPUT_PATH + "/{sra}/Solo.out/Gene/raw/barcodes.tsv",
        experimentfolder=OUTPUT_PATH + "/{sra}/Solo.out/Gene/raw/",
        experimentcellidfile=OUTPUT_PATH + "/{sra}/Solo.out/Gene/filtered.cellid/barcodes.tsv.gz",
        experimentcellidfolder=OUTPUT_PATH + "/{sra}/Solo.out/Gene/raw.cellid/",
        qthreads = 1,
        h_vmem = "2G"
    conda:
        "/path/to/project-248-FGC_RNAseq/src/STARsoloPipe/STARsolo.yml"
    shell:
        """
        mkdir -p {params.experimentcellidfolder} &&
        cp {params.experimentfolder}/features.tsv  {params.experimentcellidfolder} &&
        cp {params.experimentfolder}/matrix.mtx  {params.experimentcellidfolder} &&
        gzip {params.experimentcellidfolder}* &&
        python barcode2CellId.py \
        --srrid {params.srrid} \
        --srr2sexweekfile {params.srr2sexweekfile} \
        --barcd2cellidfile {params.barcd2cellidfile} \
        --experimentbarcdfile {params.experimentbarcdfile} \
        --experimentcellidfile {output} | gzip
        """
