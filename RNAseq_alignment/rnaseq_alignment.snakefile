# Workflow for aligning RNA sequence data to a sex specific reference with
# STAR, indexing BAM files with Samtools, and quantifying read counts with
# HTSeq.

configfile: "tcga_lihc_testing.config.json"

# Tools
STAR = "STAR"
SAMTOOLS = "samtools"
HTSEQ = "htseq-count"

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
FQ_DIR = config["test_dir"] # path to directory with fastq files
AL_DIR = config["tcga_lihc_rnaseq_realignment_dir"] # path to directory for alignment files

# Samples
XX_SAMPLES = config["tcga_lihc_rnaseq_females"]
XY_SAMPLES= config["tcga_lihc_rnaseq_males"]
ALL_SAMPLES = config["tcga_lihc_rnaseq_samples"]

SAMPLES = ["OBG0044-1-010_S37_L005"]

#MAPPED = [AL_DIR + f + "_Aligned.sortedByCoord.out.bam" for f in SAMPLES]

# TODO:
# Fix output file name format

rule all:
	input:
		# Defining the files that snakemake will attempt to produce as an output.
		# If there is no rule defined to produce the file, or if the file already
		# exists, snakemake will throw "Nothing to be done"
		#MAPPED,
		expand(AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam", AL_DIR=AL_DIR, sample=SAMPLES),
		expand(AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam.bai", AL_DIR=AL_DIR, sample=SAMPLES),
		expand(AL_DIR + "{sample}_raw_gene_counts.tsv", AL_DIR=AL_DIR, sample=SAMPLES)
		#expand(AL_DIR + "{sample}_raw_gene_counts.tsv", AL_DIR=AL_DIR, sample=SAMPLES)

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

#rule xx_align_reads:
#	input:
#		R1 = lambda wildcards: FQ_DIR + config["tcga_lihc_rnaseq_stripped_fqs"][wildcards.sample][0],
#		R2 = lambda wildcards: FQ_DIR + config["tcga_lihc_rnaseq_stripped_fqs"][wildcards.sample][1],
#		star_index = XX_STAR_INDEX,
#		gtf = XX_GTF
#	output: AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam"
#	params:
#		threads=24,
#		read_files_command = "zcat",
#	message: "Mapping {sample} reads to {REF_NAME} with STAR. Threads: {params.threads}"
#	shell:
#		"""
#		{STAR} --runThreadN {params.threads}
#		--runMode alignReads
#		--genomeDir {input.star_index}
#		--readFilesCommand {params.read_files_command}
#		--outFileNamePrefix {input.AL_DIR}{wildcards.sample}_
#		--readFilesIn {input.R1} {input.R2}
#		--outSAMtype BAM SortedByCoordinate
#		--sjdbGTFfile {input.gtf}
#		"""

rule xy_align_reads:
	input:
		R1 = lambda wildcards: FQ_DIR + config["test_fqs"][wildcards.sample][0],
		R2 = lambda wildcards: FQ_DIR + config["test_fqs"][wildcards.sample][1],
		star_index = XY_STAR_INDEX,
		AL_DIR = AL_DIR,
		gtf = XY_GTF
	output: AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam"
	params:
		threads = 24,
		read_files_command = "zcat",
	message: "Mapping {wildcards.sample} reads to {XY_REF_NAME} with STAR. \
	Threads: {params.threads}"
	shell:
		"""
		{STAR} --runThreadN {params.threads} --runMode alignReads \
		--genomeDir {input.star_index} \
		--readFilesCommand {params.read_files_command} \
		--outFileNamePrefix {input.AL_DIR}{wildcards.sample}_ \
		--readFilesIn {input.R1} {input.R2} \
		--outSAMtype BAM SortedByCoordinate \
		--sjdbGTFfile {input.gtf}
		"""

rule index_bam:
	input:
		BAM = AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam"
	output: AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam.bai"
	message: "Indexing BAM file {input.BAM} with Samtools."
	params:
	shell:
		"""
		{SAMTOOLS} index {input.BAM}
		"""

rule xy_quantify_readcounts:
	input:
		BAM = AL_DIR + "{sample}_Aligned.sortedByCoord.out.bam",
		GTF = XY_GTF
	output:
		COUNTS = AL_DIR + "{sample}_raw_gene_counts.tsv"
	params:
		FORMAT = "bam", # Format of the input data
		MODE = "union", # Mode to handle reads overlapping more than one feature
		TYPE = "exon", # Feature type (3rd column in GFF file) to be used, all features of other type are ignored
		ID_ATTRIBUTE = "gene_id", # GFF attribute to be used as feature ID
		ORDER = "pos", # For paired-end data, the alignment have to be sorted either by read name or by alignment position.
		STRANDED = "no" # Whether the data is from a strand-specific assay
	message: "Quantifying read counts with HTSeq"
	shell:
		"""
		{SAMTOOLS} view -h | \
		{HTSEQ} -m {params.MODE} -t {params.TYPE} -i \
		{params.ID_ATTRIBUTE} -f {params.FORMAT} -r {params.ORDER} \
		-s {params.STRANDED} - {input.GTF} > {output.COUNTS}
		"""
