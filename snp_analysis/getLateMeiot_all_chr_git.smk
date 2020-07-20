# Title: run bcf tools per chromosome
# Author: I. Moustakas

OUTPUT_PATH = "/path/to/project-248-FGC_RNAseq/analysis/SNP_calling_masked_no_read_QC"
CHROMOSOMES = ["1", "2", "3", "4", "5",
               "6", "7", "8", "9", "10",
               "11", "12", "13", "14", "15",
               "16", "17", "18", "19", "20",
               "21", "22", "X"]



rule all:
    input: OUTPUT_PATH + "/bcftools/lateMeioticCells_all_chr_masked_intersect_imprinted.vcf"

rule bam_files:
    output: OUTPUT_PATH + "/bam_files.txt"
    params:
        path="/path/to/project-248-FGC_RNAseq/analysis/SNP_calling_masked_no_read_QC/demultiplexed_cell_ids_dedup",
        sthreads = 1,
        mem = "4000",
        time = "1:0:0"
    shell:
        """
        cellIDs=(F_18W_embryo1_sc37 F_18W_embryo1_sc66 F_18W_embryo1_sc69 F_12W_embryo1_sc61 F_18W_embryo2_sc31 F_18W_embryo2_sc57 F_18W_embryo2_sc60 F_18W_embryo2_sc65 F_18W_embryo2_sc70 F_20W_embryo1_sc31 F_20W_embryo1_sc36 F_20W_embryo1_sc39 F_20W_embryo1_sc44 F_20W_embryo1_sc47 F_20W_embryo1_sc49 F_20W_embryo1_sc51 F_20W_embryo1_sc52 F_20W_embryo1_sc54 F_20W_embryo1_sc55 F_20W_embryo1_sc56 F_20W_embryo1_sc58 F_20W_embryo1_sc61 F_20W_embryo1_sc62 F_20W_embryo1_sc63 F_20W_embryo1_sc64 F_20W_embryo1_sc65 F_20W_embryo1_sc76 F_20W_embryo2_sc12 F_20W_embryo2_sc46 F_20W_embryo2_sc49 F_20W_embryo2_sc50 F_20W_embryo2_sc52 F_20W_embryo2_sc53 F_20W_embryo2_sc54 F_20W_embryo2_sc55 F_20W_embryo2_sc58 F_20W_embryo2_sc59 F_20W_embryo2_sc70 F_20W_embryo2_sc71 F_20W_embryo2_sc72 F_20W_embryo2_sc73 F_20W_embryo2_sc76 F_20W_embryo2_sc77 F_23W_embryo1_sc40 F_23W_embryo1_sc41 F_23W_embryo1_sc43 F_23W_embryo1_sc44 F_23W_embryo1_sc45 F_23W_embryo1_sc47 F_23W_embryo1_sc48 F_23W_embryo1_sc49 F_23W_embryo1_sc50 F_23W_embryo1_sc51 F_23W_embryo1_sc52 F_23W_embryo1_sc53 F_23W_embryo1_sc54 F_23W_embryo1_sc55 F_23W_embryo1_sc56 F_23W_embryo1_sc57 F_23W_embryo1_sc58 F_23W_embryo1_sc59 F_23W_embryo1_sc60 F_23W_embryo1_sc61 F_23W_embryo1_sc62 F_23W_embryo1_sc63 F_23W_embryo1_sc64 F_23W_embryo1_sc65 F_23W_embryo1_sc66 F_23W_embryo1_sc67 F_23W_embryo1_sc68 F_23W_embryo1_sc70 F_23W_embryo1_sc72 F_23W_embryo1_sc74 F_23W_embryo2_sc4 F_23W_embryo2_sc39 F_23W_embryo2_sc40 F_23W_embryo2_sc41 F_23W_embryo2_sc43 F_23W_embryo2_sc44 F_23W_embryo2_sc45 F_23W_embryo2_sc47 F_23W_embryo2_sc48 F_23W_embryo2_sc49 F_23W_embryo2_sc51 F_23W_embryo2_sc52 F_23W_embryo2_sc53 F_23W_embryo2_sc56 F_23W_embryo2_sc57 F_23W_embryo2_sc59 F_23W_embryo2_sc60 F_23W_embryo2_sc61 F_23W_embryo2_sc62 F_23W_embryo2_sc67 F_24W_embryo1_sc39 F_24W_embryo1_sc50 F_24W_embryo1_sc54 F_24W_embryo1_sc55 F_24W_embryo1_sc56 F_24W_embryo1_sc64 F_24W_embryo1_sc66 F_24W_embryo1_sc70 F_24W_embryo2_sc27 F_24W_embryo2_sc39 F_24W_embryo2_sc40 F_24W_embryo2_sc41 F_24W_embryo2_sc44 F_24W_embryo2_sc48 F_24W_embryo2_sc49 F_24W_embryo2_sc50 F_24W_embryo2_sc54 F_24W_embryo2_sc55 F_24W_embryo2_sc63 F_24W_embryo2_sc64 F_24W_embryo2_sc65 F_24W_embryo2_sc66 F_24W_embryo2_sc67 F_26W_embryo1_sc41 F_11W_embryo1_sc33 F_11W_embryo1_sc39 F_18W_embryo1_sc19 F_18W_embryo1_sc29 F_18W_embryo1_sc31 F_12W_embryo1_sc43 F_12W_embryo1_sc44 F_12W_embryo1_sc45 F_12W_embryo1_sc50 F_12W_embryo1_sc53 F_12W_embryo1_sc60 F_12W_embryo1_sc66 F_12W_embryo1_sc68 F_12W_embryo1_sc83 F_12W_embryo1_sc90 F_12W_embryo1_sc94 F_18W_embryo2_sc17 F_20W_embryo1_sc26 F_20W_embryo2_sc25 F_20W_embryo2_sc28 F_20W_embryo2_sc29 F_20W_embryo2_sc45 F_20W_embryo2_sc56 F_20W_embryo2_sc63 F_23W_embryo1_sc1 F_23W_embryo1_sc3 F_23W_embryo1_sc6 F_23W_embryo1_sc7 F_23W_embryo1_sc8 F_23W_embryo1_sc9 F_23W_embryo1_sc10 F_23W_embryo1_sc12 F_23W_embryo1_sc18 F_23W_embryo1_sc21 F_23W_embryo1_sc22 F_23W_embryo1_sc24 F_23W_embryo1_sc25 F_23W_embryo1_sc27 F_23W_embryo1_sc29 F_23W_embryo1_sc32 F_23W_embryo1_sc33 F_23W_embryo1_sc36 F_23W_embryo2_sc2 F_23W_embryo2_sc3 F_23W_embryo2_sc7 F_23W_embryo2_sc10 F_23W_embryo2_sc13 F_23W_embryo2_sc17 F_23W_embryo2_sc20 F_23W_embryo2_sc22 F_23W_embryo2_sc25 F_23W_embryo2_sc28 F_24W_embryo2_sc10 F_24W_embryo2_sc43 F_26W_embryo1_sc48 F_26W_embryo1_sc78 F_11W_embryo1_sc4 F_11W_embryo1_sc5 F_11W_embryo1_sc6 F_11W_embryo1_sc9 F_11W_embryo1_sc12 F_11W_embryo1_sc13 F_11W_embryo1_sc20 F_11W_embryo1_sc22 F_11W_embryo1_sc38 F_18W_embryo1_sc11 F_18W_embryo2_sc2 F_18W_embryo2_sc4 F_18W_embryo2_sc7 F_18W_embryo2_sc21 F_18W_embryo2_sc30 F_20W_embryo1_sc20 F_20W_embryo2_sc8 F_20W_embryo2_sc9 F_20W_embryo2_sc39 F_20W_embryo2_sc47 F_20W_embryo2_sc74 F_23W_embryo1_sc16 F_23W_embryo1_sc17 F_23W_embryo1_sc26 F_23W_embryo1_sc37 F_23W_embryo1_sc38 F_23W_embryo1_sc39 F_23W_embryo1_sc42 F_23W_embryo2_sc6 F_23W_embryo2_sc9 F_23W_embryo2_sc24 F_23W_embryo2_sc26 F_23W_embryo2_sc29 F_23W_embryo2_sc35 F_23W_embryo2_sc54 F_24W_embryo1_sc4 F_24W_embryo1_sc31 F_24W_embryo1_sc40 F_24W_embryo1_sc68 F_24W_embryo2_sc7 F_24W_embryo2_sc12 F_24W_embryo2_sc14 F_24W_embryo2_sc15 F_24W_embryo2_sc20 F_24W_embryo2_sc25 F_24W_embryo2_sc28 F_24W_embryo2_sc32 F_24W_embryo2_sc33 F_24W_embryo2_sc45 F_24W_embryo2_sc56 F_26W_embryo1_sc14 F_26W_embryo1_sc65 F_26W_embryo1_sc68 F_26W_embryo1_sc75)

        touch bam_files.txt
        declare -a files
        for id in "${{cellIDs[@]}}"; do
            find {params.path} -name $id.dedup.bam -exec readlink -f {{}} \; >> {output}
        done
        """


rule bcftools:
    input: OUTPUT_PATH + "/bam_files.txt"
    output: OUTPUT_PATH + "/bcftools/lateMeioticCells_chr{chrom}_masked.vcf.gz"
    params:
        chrom = "{chrom}",
        sthreads = 2,
        mem = "4000",
        time = "20:0:0"
    conda:
        "/path/to/project-248-FGC_RNAseq/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        files=(`cat bam_files.txt`)

        bcftools mpileup -Ou \
        --max-depth 1000 \
        --regions {params.chrom} \
        --annotate FORMAT/AD,FORMAT/DP \
        -f /path/to/ioannis/cellRangerPipelines/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa \
        ${{files[@]}} |\
        bcftools call -m \
        -Ou \
        --variants-only |\
        bcftools filter \
        -Oz \
        --include 'COUNT(FORMAT/DP[*]>=10)>=5' \
        -o {output}
        """

rule concat:
    input: [OUTPUT_PATH + "/bcftools/lateMeioticCells_chr{chrom}_masked.vcf.gz".format(chrom=chrom) for chrom in CHROMOSOMES]
    output: OUTPUT_PATH + "/bcftools/lateMeioticCells_all_chr_masked.vcf.gz"
    params:
        sthreads = 1,
        mem = "4000",
        time = "1:0:0"
    conda:
        "/path/to/project-248-FGC_RNAseq/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        bcftools concat \
        -Oz \
        -o {output} \
        {input}
        """

rule intersect:
    input:
        vcf=OUTPUT_PATH + "/bcftools/lateMeioticCells_all_chr_masked.vcf.gz",
        bed="/path/to/project-248-FGC_RNAseq/data/imprinted_genes/imprinted_genes.bed"
    output: OUTPUT_PATH + "/bcftools/lateMeioticCells_all_chr_masked_intersect_imprinted.vcf"
    params:
        sthreads = 1,
        mem = "4000",
        time = "1:0:0"
    conda:
        "/path/to/project-248-FGC_RNAseq/src/SNP_calling/samtools_snpCall.yml"
    shell:
        """
        /path/to/ioannis/miniconda3/envs/samtools_snpCall/bin/bedtools intersect \
        -header \
        -a {input.vcf} \
        -b {input.bed} \
	-u \
        > {output}
        """
