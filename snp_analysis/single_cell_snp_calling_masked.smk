# Align files with STARsolo
# Demultiplex the alignment files and output one alignment file per cell
# Deduplicate the individual cell files usung the UMI information
# Author: I. Moustakas

OUTPUT_PATH = "/path/to/analysis/SNP_calling_masked_no_read_QC"
INPUT_PATH = "/epath/to/analysis/STARsolo_masked_no_read_QC"

#SRAs, = glob_wildcards(INPUT_PATH + "/{sra}_pass_1.fastq.gz")

## SRAs selected to include only female embryos
SRAs = ["SRR4199320", "SRR4199311", "SRR4081201", "SRR4199312", "SRR4199313",
        "SRR4199314", "SRR4199316", "SRR4199317", "SRR4199318", "SRR4199315",
        "SRR4199319"]

rule all:
    input: [OUTPUT_PATH + "/demultiplexed_cell_ids_sorted_flagstat/{sra}/".format(sra= sra) for sra in SRAs],
           [OUTPUT_PATH + "/demultiplexed_cell_ids_dedup_flagstat/{sra}/".format(sra= sra) for sra in SRAs]

rule samtools:
    input: INPUT_PATH + "/{sra}/Aligned.sortedByCoord.out.bam"
    output: INPUT_PATH + "/{sra}/Aligned.sortedByCoord.out.sorted.bam"
    params:
        sthreads = 2,
        mem = "2000",
        time = "1:0:0"
    conda:
        "/path/to/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        samtools sort -@ 2 -o {output} {input} && \
        samtools index {output}
        """

rule demultiplex:
    input: INPUT_PATH + "/{sra}/Aligned.sortedByCoord.out.sorted.bam"
    output: OUTPUT_PATH + "/demultiplexed/{sra}/demultiplexing_done.txt"
    params:
        barcodes = "/path/to/src/HECA/scRNA-seq_pipeline_hg38/barcode_96_8bp_noCells.csv",
        output_dir = OUTPUT_PATH + "/demultiplexed/{sra}/",
        sthreads = 2,
        mem = "4000",
        time = "2:0:0"
    conda:
        "/path/to/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        python splitOnCellBarcode.py \
        --barcodes {params.barcodes} \
        --multiplexed_bam {input} \
        --output_dir {params.output_dir}
        """

checkpoint barcode2CellId:
    input: OUTPUT_PATH + "/demultiplexed/{sra}/demultiplexing_done.txt"
    output: directory(OUTPUT_PATH + "/demultiplexed_cell_ids/{sra}/")
    params:
        srrid="{sra}",
        srr2sexweekfile="/path/to//data/LiLiData/SRRid2sexWeek.csv",
        barcd2cellidfile="/path/to/src/HECA/scRNA-seq_pipeline_hg38/barcode_96_8bp.txt",
        path2barcodes=OUTPUT_PATH + "/demultiplexed/{sra}/",
        path2cellids=OUTPUT_PATH + "/demultiplexed_cell_ids/{sra}/",
        sthreads = 1,
        mem = "1000",
        time = "1:0:0"
    conda:
        "/path/to/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        python barcode2CellId.py \
        --srrid {params.srrid} \
        --srr2sexweekfile {params.srr2sexweekfile} \
        --barcd2cellidfile {params.barcd2cellidfile} \
        --path2barcodes {params.path2barcodes} \
        --path2cellids {params.path2cellids}
        """

checkpoint add_readgroup:
    input: directory(OUTPUT_PATH + "/demultiplexed_cell_ids/{sra}/")
    output: directory(OUTPUT_PATH + "/demultiplexed_cell_ids_readgroup/{sra}/")
    params:
        input_bam_path=OUTPUT_PATH + "/demultiplexed_cell_ids/{sra}",
        output_bam_path=OUTPUT_PATH + "/demultiplexed_cell_ids_readgroup/{sra}/",
        sthreads = 2,
        mem = "3000",
        time = "3:0:0"
    conda:
        "/path/to/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        for f in `ls {params.input_bam_path}/*bam`;
        do
            basename=$(echo $f | sed -e 's/\.bam$//' | sed -e 's/.*\///g')
            java -Xmx4G \
            -jar /exports/sasc/ioannis/gatk/gatk-4.1.7.0/gatk-package-4.1.7.0-local.jar \
            AddOrReplaceReadGroups \
            -I $f \
            -O {params.output_bam_path}$basename.bam \
            --RGLB $basename \
            --RGPL ILLUMINA \
            --RGPU $basename \
            --RGSM $basename
        done
        """

checkpoint sort:
    input: directory(OUTPUT_PATH + "/demultiplexed_cell_ids_readgroup/{sra}/")
    output: directory(OUTPUT_PATH + "/demultiplexed_cell_ids_sorted/{sra}/")
    params:
        input_bam_path=OUTPUT_PATH + "/demultiplexed_cell_ids_readgroup/{sra}/",
        output_bam_path=OUTPUT_PATH + "/demultiplexed_cell_ids_sorted/{sra}/",
        sthreads = 2,
        mem = "3000",
        time = "3:0:0"
    conda:
        "/path/to/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        for f in `ls {params.input_bam_path}/*bam`;
        do
            basename=$(echo $f | sed -e 's/\.bam$//' | sed -e 's/.*\///g')
            samtools sort -@ 2 $f -o {params.output_bam_path}$basename.sorted.bam && samtools index {params.output_bam_path}$basename.sorted.bam
        done
        """

checkpoint flagstats:
    input: directory(OUTPUT_PATH + "/demultiplexed_cell_ids_sorted/{sra}/")
    output: directory(OUTPUT_PATH + "/demultiplexed_cell_ids_sorted_flagstat/{sra}/")
    params:
        input_bam_path=OUTPUT_PATH + "/demultiplexed_cell_ids_sorted/{sra}",
        output_bam_path=OUTPUT_PATH + "/demultiplexed_cell_ids_sorted_flagstat/{sra}/",
        sthreads = 2,
        mem = "3000",
        time = "1:0:0"
    conda:
        "/path/to/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        for f in `ls {params.input_bam_path}/*bam`;
        do
            basename=$(echo $f | sed -e 's/\.bam$//' | sed -e 's/.*\///g')
            samtools flagstat -@ 2 $f > {params.output_bam_path}/$basename.flagstat
        done
        """



## https://umi-tools.readthedocs.io/en/latest/reference/dedup.html
checkpoint umitools:
    input: directory(OUTPUT_PATH + "/demultiplexed_cell_ids_sorted/{sra}/")
    output: directory(OUTPUT_PATH + "/demultiplexed_cell_ids_dedup/{sra}/")
    params:
        input_bam_path=OUTPUT_PATH + "/demultiplexed_cell_ids_sorted/{sra}",
        output_bam_path=OUTPUT_PATH + "/demultiplexed_cell_ids_dedup/{sra}/",
        log_path=OUTPUT_PATH + "/demultiplexed_cell_ids_dedup/{sra}/logs",
        sthreads = 1,
        mem = "7000",
        time = "6:0:0"
    conda:
        "/path/to/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        mkdir -p {params.log_path}
        for f in `ls {params.input_bam_path}/*.sorted.bam`;
        do
            basename=$(echo $f | sed -e 's/\.sorted.bam$//' | sed -e 's/.*\///g')
            umi_tools dedup \
            --stdin=$f \
            --log={params.log_path}/$basename.log_umi.txt \
            --extract-umi-method=tag \
            --umi-tag=UR \
            --method directional \
            > {params.output_bam_path}/$basename.dedup.bam
            samtools index {params.output_bam_path}/$basename.dedup.bam
        done
        """

checkpoint flagstats_dedupl:
    input: directory(OUTPUT_PATH + "/demultiplexed_cell_ids_dedup/{sra}/")
    output: directory(OUTPUT_PATH + "/demultiplexed_cell_ids_dedup_flagstat/{sra}/")
    params:
        input_bam_path=OUTPUT_PATH + "/demultiplexed_cell_ids_dedup/{sra}",
        output_bam_path=OUTPUT_PATH + "/demultiplexed_cell_ids_dedup_flagstat/{sra}/",
        sthreads = 2,
        mem = "3000",
        time = "1:0:0"
    conda:
        "/path/to/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        for f in `ls {params.input_bam_path}/*bam`;
        do
            basename=$(echo $f | sed -e 's/\.bam$//' | sed -e 's/.*\///g')
            samtools flagstat -@ 2 $f > {params.output_bam_path}/$basename.flagstat
        done
        """

# Call SNPs and do a first filtering on the results
rule bcftools:
    input:  directory(OUTPUT_PATH + "/demultiplexed_cell_ids_dedup/{sra}/")
    output: OUTPUT_PATH + "/vcf_files_chr_X/{sra}.vcf"
    params:
        input_dir = OUTPUT_PATH + "/demultiplexed_cell_ids_dedup/{sra}/",
        sthreads = 1,
        mem = "4000",
        time = "48:0:0"
    conda:
        "/path/to/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        bcftools mpileup -Ou \
        --max-depth 1000 \
        --regions "X" \
        --annotate FORMAT/AD,FORMAT/DP \
        -f /path/to/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
        `ls {params.input_dir}/*.dedup.bam` |\
        bcftools call -m \
        -Ou \
        --gvcf 5 |\
        bcftools filter \
        -Oz \
        --include 'FORMAT/DP[*]>8' \
        -o {output}
        """
