# Workflow for aligning RNA sequence data to a sex specific reference with
# STAR, indexing BAM files with Samtools, and quantifying read counts with
# HTSeq.

configfile: "tcga_lihc_testing.config.json"

# Tools
STAR = "STAR"
SAMTOOLS = "samtools"
HTSEQ = "" # Path to HTSeq script htseq-count

# Reference genome files: XX with Y chromosome masked, XY with both X and Y
XX_REF = config["xx_GRCh38_ref_path"]
XX_REF_NAME = config["xx_GRCh38_ref_prefix"]
XX_GTF = config["xx_gtf_path"]
XX_STAR_INDEX = config["xx_star_index_dir"]
XY_REF = config["xy_GRCh38_ref_path"]
XY_REF_NAME = config["xy_GRCh38_ref_prefix"]
XY_GTF = config["xy_gtf_path"]
XY_STAR_INDEX = config["xy_star_index_dir"]

# Directories
FQ_DIR = config["tcga_lihc_rnaseq_stripped_fq_dir"] # path to directory with fastq files
AL_DIR = config["tcga_lihc_rnaseq_realignment_dir"] # path to directory for alignment files

# Samples
XX_SAMPLES = config["tcga_lihc_rnaseq_females"]
XY_SAMPLES = config["tcga_lihc_rnaseq_males"]
ALL_SAMPLES = config["tcga_lihc_rnaseq_samples"]

rule all:
	input:
    	#expand(FQ_DIR + "{sample}_rna_seq_1", FQ_DIR=FQ_DIR, sample=config["tcga_lihc_rnaseq_females"]),
    	#expand(FQ_DIR + "{sample}_rna_seq_2", FQ_DIR=FQ_DIR, sample=config["tcga_lihc_rnaseq_females"]),
        #expand(["{FQ_DIR}{sample}_rna_seq_1.fastq", "{FQ_DIR}{sample}_rna_seq_2.fastq"], FQ_DIR=FQ_DIR, sample=config["tcga_lihc_rnaseq_paired_end_samples"]"),
        #expand("{FQ_DIR}{sample}_rna_seq_1.fastq", FQ_DIR=FQ_DIR, sample=config["tcga_lihc_rnaseq_paired_end_samples"]),
		#expand("{FQ_DIR}{sample}_rna_seq_2.fastq", FQ_DIR=FQ_DIR, sample=config["tcga_lihc_rnaseq_paired_end_samples"]),
		expand(AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam", AL_DIR=AL_DIR, sample=config["tcga_lihc_rnaseq_females"]),
		#expand("{sample}_", sample=config["tcga_lihc_rnaseq_paired_end_samples"])
	output:
		#expand(AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam", AL_DIR=AL_DIR, sample=config["tcga_lihc_rnaseq_paired_end_samples"])
		#expand("{AL_DIR}{sample}_Aligned.sortedByCoord.out.bam", AL_DIR=AL_DIR, sample=config["tcga_lihc_rnaseq_paired_end_samples"]),
		#expand("{AL_DIR}{sample}_Aligned.sortedByCoord.out.bam.bai", AL_DIR=AL_DIR, sample=config["tcga_lihc_rnaseq_paired_end_samples"])
	params:
		expand("{AL_DIR}{sample}_", AL_DIR=AL_DIR, sample=config["tcga_lihc_rnaseq_females"])


#rule star_index:
#	input:
#		ref=REF,
#		star_index=star_index
#	output:
#	params:
#		threads=24
#   message: "Generating STAR genome index files for reference {REF} and GTF {GTF}"
#	shell:
#		"""
#		"{STAR} --runThreadN {params.threads} --runMode genomeGenerate "
#		"--genomeDir {star_index} --genomeFastaFiles {input.ref} --sjdbGTFfile {gtf}; "
#		"""

rule xx_align_reads:
	input:
		R1=FQ_DIR + "{sample}_1.fastq.gz",
		R2=FQ_DIR + "{sample}_2.fastq.gz",
		star_index=STAR_INDEX,
		gtf=GTF
	output: AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam"
	params:
		threads=24,
		name_prefix=AL_DIR + "{sample}_"
	message: "Mapping {sample} reads to {REF_NAME} with STAR. Threads: {params.threads}"
	shell:
		"""
		{STAR} --runThreadN {params.threads} --runMode alignReads
		--genomeDir {input.star_index} --outFileNamePrefix {params.name_prefix}
		--readFilesIn {input.R1} {input.R2} --outSAMtype BAM SortedByCoordinate
		--sjdbGTFfile {input.gtf}
		"""

rule xy_align_reads:
	input:
		R1= FQ_DIR + "{sample}_rna_seq_1.fastq",
		R2= FQ_DIR + "{sample}_rna_seq_2.fastq",
		star_index=STAR_INDEX,
		gtf=GTF
	output: AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam"
	params:
		threads=24,
		name_prefix=AL_DIR + "{sample}_"
	message: "Mapping {sample} reads to {REF_NAME} with STAR. Threads: {params.threads}"
	shell:
		"""
		{STAR} --runThreadN {params.threads} --runMode alignReads
		--genomeDir {input.star_index} --outFileNamePrefix {params.name_prefix}
		--readFilesIn {input.R1} {input.R2} --outSAMtype BAM SortedByCoordinate
		--sjdbGTFfile {input.gtf}
		"""

rule index_bam:
	input:
		BAM="{AL_DIR}{sample}_Aligned.sortedByCoord.out.bam"
	output: "{AL_DIR}{sample}_Aligned.sortedByCoord.out.bam.bai"
	message: "Indexing BAM file {BAM} with Samtools."
	params:
	shell:
		"""
		{SAMTOOLS} index {input.BAM}
		"""

rule quantify_readcounts:
   input:
        BAM="{}",
        GET="{GTF}"
    output:
        COUNTS="{SAMPLE}_raw_gene_counts.tsv"
    params:
        FORMAT="bam", # Format of the input data
        MODE="union", # Mode to handle reads overlapping more than one feature
        TYPE="exon", # Feature type (3rd column in GFF file) to be used, all features of other type are ignored
        ID_ATTRIBUTE="gene_id", # GFF attribute to be used as feature ID
        ORDER="pos", # For paired-end data, the alignment have to be sorted either by read name or by alignment position.
        STRANDED="no" # Whether the data is from a strand-specific assay
    message: "Quantifying read counts with HTSeq"
    shell:
        """
        python {HTSEQ} -m {params.MODE} -t {params.TYPE} -i
        {params.ID_ATTRIBUTE} -f {params.FORMAT} -r {params.ORDER}
        -s {params.STRANDED} {input.BAM} {input.GTF} > {output.COUNTS}
        """
