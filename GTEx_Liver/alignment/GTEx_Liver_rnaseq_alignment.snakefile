# Workflow for aligning RNA sequence data to a sex specific reference with
# HISAT2, sorting and indexing BAM files with Samtools, and quantifying
# read counts with Subread featureCounts.

configfile: "GTEx_v1.2.config.json"

# Tools
STAR = "STAR"
SAMTOOLS = "samtools"
HTSEQ = "htseq-count"
FEATURECOUNTS = "featureCounts"

# Reference genome files: XX with Y chromosome masked, XY with both X and Y
XX_REF = config["XX_GRCh38_ref_path"]
XX_REF_NAME = config["XX_GRCh38_ref_prefix"]
XX_GTF = config["XX_gtf_path"]
XX_STAR_INDEX = config["XX_star_index_dir"]
XX_HISAT2_INDEX = config["XX_hisat2_index"]
XX_WITH_VIRAL_REF = config["XX_GRCh38_ref_with_viral_genomes"]
XX_HISAT2_INDEX_WITH_VIRAL_REF = config["XX_HISAT2_index_GRCh38_ref_with_viral_genomes"]
XY_REF = config["XY_GRCh38_ref_path"]
XY_REF_NAME = config["XY_GRCh38_ref_prefix"]
XY_GTF = config["XY_gtf_path"]
XY_STAR_INDEX = config["XY_star_index_dir"]
XY_HISAT2_INDEX = config["XY_hisat2_index"]
XY_WITH_VIRAL_REF = config["XY_GRCh38_ref_with_viral_genomes"]
XY_HISAT2_INDEX_WITH_VIRAL_REF = config["XY_HISAT2_index_GRCh38_ref_with_viral_genomes"]


# Directories
FQ_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Liver/RNAseq/trimmed_fastqs/" # path to directory with fastq files
SAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Liver/RNAseq/GRCh38_SAM/" # path to directory for SAM alignment files
BAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Liver/RNAseq/GRCh38_BAM/" # path to directory for BAM alignment files
SORTED_BAM_AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Liver/RNAseq/GRCh38_sorted_BAM/" # path to directory for sorted BAM alignment files
HTSEQ_COUNTS_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Liver/RNAseq/GRCh38_HTSeq_gene_counts/" # path to directory for HTSeq gene count files
FEATURECOUNTS_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Liver/RNAseq/GRCh38_featurecounts/" # path to directory for FeatureCounts gene count files


# Samples
XX_SAMPLES = config["Liver_Female_RNA"]
XY_SAMPLES= config["Liver_Male_RNA"]
SAMPLES = config["Liver_RNA"]

#MAPPED = [AL_DIR + f + "_Aligned.sortedByCoord.out.bam" for f in SAMPLES]

rule all:
	input:
		# Defining the files that snakemake will attempt to produce as an output.
		# If there is no rule defined to produce the file, or if the file already
		# exists, snakemake will throw "Nothing to be done"
		#MAPPED,
		expand(SAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38.sam", SAM_AL_DIR=SAM_AL_DIR, sample=XX_SAMPLES),
		expand(BAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38.bam", BAM_AL_DIR=BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38_sortedbycoord.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38_sortedbycoord.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(SAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38.sam", SAM_AL_DIR=SAM_AL_DIR, sample=XY_SAMPLES),
		expand(BAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38.bam", BAM_AL_DIR=BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38_sortedbycoord.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		expand(SORTED_BAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38_sortedbycoord.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		#expand(SAM_AL_DIR + "{sample}_XX_Aligned.out.sam", SAM_AL_DIR=SAM_AL_DIR, sample=XX_SAMPLES),
		#expand(BAM_AL_DIR + "{sample}_XX_Aligned.out.bam", BAM_AL_DIR=BAM_AL_DIR, sample=XX_SAMPLES),
		#expand(SAM_AL_DIR + "{sample}_XY_Aligned.out.sam", SAM_AL_DIR=SAM_AL_DIR, sample=XY_SAMPLES),
		#expand(BAM_AL_DIR + "{sample}_XY_Aligned.out.bam", BAM_AL_DIR=BAM_AL_DIR, sample=XY_SAMPLES),
		#expand(SORTED_BAM_AL_DIR + "{sample}_XX_Aligned.sortedByCoord.out.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		#expand(SORTED_BAM_AL_DIR + "{sample}_XY_Aligned.sortedByCoord.out.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		#expand(SORTED_BAM_AL_DIR + "{sample}_XY_Aligned.sortedByCoord.out.bam", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XY_SAMPLES),
		#expand(SORTED_BAM_AL_DIR + "{sample}_XX_Aligned.sortedByCoord.out.bam.bai", SORTED_BAM_AL_DIR=SORTED_BAM_AL_DIR, sample=XX_SAMPLES),
		expand(FEATURECOUNTS_DIR + "{sample}_RNA_XX_HISAT2_transcript_featurecounts.tsv", HTSEQ_COUNTS_DIR=HTSEQ_COUNTS_DIR, sample=XX_SAMPLES),
		expand(FEATURECOUNTS_DIR + "{sample}_RNA_XY_HISAT2_transcript_featurecounts.tsv", HTSEQ_COUNTS_DIR=HTSEQ_COUNTS_DIR, sample=XY_SAMPLES),
		#expand(HTSEQ_COUNTS_DIR + "{sample}_XY_raw_gene_counts.tsv", HTSEQ_COUNTS_DIR=HTSEQ_COUNTS_DIR, sample=XY_SAMPLES)

#print (expand(SAM_AL_DIR + "{sample}_Aligned.out.sam", SAM_AL_DIR=SAM_AL_DIR, sample=XX_SAMPLES))

# rule hisat2_index:
# 	input:
#	REF = ""
# 	output:
# 	params:
# 	message:
# 	shell:
# 		"""
# 		hisat2-build -f {input.REF} {input.REF_PREFIX}
# 		"""

rule hisat2_xx_align_reads:
	input:
		fq1 = FQ_DIR + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		fq2 = FQ_DIR + "{sample}_trimmomatic_trimmed_paired_2.fastq"
	output:
		SAM = SAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38.sam"
	params:
		hisat2_index = XX_HISAT2_INDEX_WITH_VIRAL_REF,
		threads = 8
	message: "Mapping {wildcards.sample} reads to {params.hisat2_index} with HISAT2."
	shell:
		"""
		hisat2 -q --phred33 -p {params.threads} -x {params.hisat2_index} -s no -1 {input.fq1} -2 {input.fq2} -S {output.SAM}
		"""

rule hisat2_xy_align_reads:
	input:
		fq1 = FQ_DIR + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		fq2 = FQ_DIR + "{sample}_trimmomatic_trimmed_paired_2.fastq"
	output:
		SAM = SAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38.sam"
	params:
		hisat2_index = XY_HISAT2_INDEX_WITH_VIRAL_REF,
		threads = 8
	message: "Mapping {wildcards.sample} reads to {params.hisat2_index} with HISAT2."
	shell:
		"""
		hisat2 -q --phred33 -p {params.threads} -x {params.hisat2_index} -s no -1 {input.fq1} -2 {input.fq2} -S {output.SAM}
		"""

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

# rule star_xx_align_reads:
# 	input:
# 		R1 = lambda wildcards: FQ_DIR + config[wildcards.sample]["RNAseq_fastq_dir"][0] + "/" + config[wildcards.sample]["RNAseq_fastq_1"][0],
# 		R2 = lambda wildcards: FQ_DIR + config[wildcards.sample]["RNAseq_fastq_dir"][0] + "/" + config[wildcards.sample]["RNAseq_fastq_2"][0],
# 		star_index = XX_STAR_INDEX,
# 		SAM_AL_DIR = SAM_AL_DIR,
# 		gtf = XX_GTF
# 	output: SAM_AL_DIR + "{sample}_XX_Aligned.out.sam"
# 	params:
# 		threads = 24,
# 		read_files_command = "", # Use zcat, if input FASTQs are gzipped
# 		out_sam_type = "SAM", #  This ”unsorted” file can be directly used with downstream software such as HTseq, without the need of name sorting.
# 		quant_mode = "GeneCounts"
# 	message: "Mapping {wildcards.sample} reads to {XY_REF_NAME} with STAR. \
# 	Threads: {params.threads}"
# 	shell:
# 		"""
# 		{STAR} --runThreadN {params.threads} --runMode alignReads \
# 		--genomeDir {input.star_index} \
# 		--readFilesCommand {params.read_files_command} \
# 		--outFileNamePrefix {input.SAM_AL_DIR}{wildcards.sample}_ \
# 		--readFilesIn {input.R1} {input.R2} \
# 		--outSAMtype {params.out_sam_type} \
# 		--quantMode {params.quant_mode} \
# 		--sjdbGTFfile {input.gtf}
# 		"""
#
# rule star_xy_align_reads:
# 	input:
# 		R1 = lambda wildcards: FQ_DIR + config[wildcards.sample]["RNAseq_fastq_dir"][0] + "/" + config[wildcards.sample]["RNAseq_fastq_1"][0],
# 		R2 = lambda wildcards: FQ_DIR + config[wildcards.sample]["RNAseq_fastq_dir"][0] + "/" + config[wildcards.sample]["RNAseq_fastq_2"][0],
# 		star_index = XY_STAR_INDEX,
# 		SAM_AL_DIR = SAM_AL_DIR,
# 		gtf = XY_GTF
# 	output: SAM_AL_DIR + "{sample}_XY_Aligned.out.sam"
# 	params:
# 		threads = 24,
# 		read_files_command = "",  # Use zcat, if input FASTQs are gzipped
# 		out_sam_type = "SAM", #  This ”unsorted” file can be directly used with downstream software such as HTseq, without the need of name sorting.
# 		quant_mode = "GeneCounts"
# 	message: "Mapping {wildcards.sample} reads to {XY_REF_NAME} with STAR. \
# 	Threads: {params.threads}"
# 	shell:
# 		"""
# 		{STAR} --runThreadN {params.threads} --runMode alignReads \
# 		--genomeDir {input.star_index} \
# 		--readFilesCommand {params.read_files_command} \
# 		--outFileNamePrefix {input.SAM_AL_DIR}{wildcards.sample}_ \
# 		--readFilesIn {input.R1} {input.R2} \
# 		--outSAMtype {params.out_sam_type} \
# 		--quantMode {params.quant_mode} \
# 		--sjdbGTFfile {input.gtf}
# 		"""

rule xx_sam_to_bam:
	input:
		SAM = SAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38.sam"
	output:
		BAM = BAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38.bam"
	params:
	message: "Converting {input.SAM} to BAM, only outputting mapped reads."
	shell:
		"""
		{SAMTOOLS} view -b -F 4 {input.SAM} > {output.BAM}
		"""

rule xy_sam_to_bam:
	input:
		SAM = SAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38.sam"
	output:
		BAM = BAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38.bam"
	params:
	message: "Converting {input.SAM} to BAM, only outputting mapped reads."
	shell:
		"""
		{SAMTOOLS} view -b -F 4 {input.SAM} > {output.BAM}
		"""

rule xx_sort_bam:
	input:
		BAM = BAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38.bam"
	output:
		SORTED_BAM = SORTED_BAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38_sortedbycoord.bam"
	params:
	message: "Sorting BAM file {input.BAM}"
	shell:
		"""
		{SAMTOOLS} sort -O bam -o {output.SORTED_BAM} {input.BAM}
		"""

rule xy_sort_bam:
	input:
		BAM = BAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38.bam"
	output:
		SORTED_BAM = SORTED_BAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38_sortedbycoord.bam"
	params:
	message: "Sorting BAM file {input.BAM}"
	shell:
		"""
		{SAMTOOLS} sort -O bam -o {output.SORTED_BAM} {input.BAM}
		"""

rule xx_index_bam:
	input:
		BAM = SORTED_BAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38_sortedbycoord.bam"
	output: SORTED_BAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38_sortedbycoord.bam.bai"
	message: "Indexing BAM file {input.BAM} with Samtools."
	params:
	shell:
		"""
		{SAMTOOLS} index {input.BAM}
		"""

rule xy_index_bam:
	input:
		BAM = SORTED_BAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38_sortedbycoord.bam"
	output: SORTED_BAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38_sortedbycoord.bam.bai"
	message: "Indexing BAM file {input.BAM} with Samtools."
	params:
	shell:
		"""
		{SAMTOOLS} index {input.BAM}
		"""

# rule xx_htseq_quantify_readcounts:
# 	input:
# 		BAM = SORTED_BAM_AL_DIR + "{sample}_XX_HISAT2_aligned_sortedbycoord.bam",
# 		GTF = XX_GTF
# 	output:
# 		COUNTS = HTSEQ_COUNTS_DIR + "{sample}_XX_HISAT2_HTSEQ_raw_gene_counts.tsv"
# 	params:
# 		FORMAT = "sam", # Format of the input data
# 		MODE = "union", # Mode to handle reads overlapping more than one feature
# 		TYPE = "exon", # Feature type (3rd column in GFF file) to be used, all features of other type are ignored
# 		ID_ATTRIBUTE = "gene_id", # GFF attribute to be used as feature ID
# 		ORDER = "pos", # For paired-end data, the alignment have to be sorted either by read name or by alignment position.
# 		STRANDED = "no" # Whether the data is from a strand-specific assay
# 	message: "Quantifying read counts from SAM file {input.BAM} with HTSeq"
# 	shell:
# 		"""
# 		{HTSEQ} -m {params.MODE} -t {params.TYPE} -i \
# 		{params.ID_ATTRIBUTE} -f {params.FORMAT} -r {params.ORDER} \
# 		-s {params.STRANDED} {input.BAM} {input.GTF} > {output.COUNTS}
# 		"""

# rule xy_htseq_quantify_readcounts:
# 	input:
# 		BAM = SORTED_BAM_AL_DIR + "{sample}_XY_HISAT2_aligned_sortedbycoord.bam",
# 		GTF = XY_GTF
# 	output:
# 		COUNTS = HTSEQ_COUNTS_DIR + "{sample}_XY_HISAT2_HTSEQ_raw_gene_counts.tsv"
# 	params:
# 		FORMAT = "bam", # Format of the input data
# 		MODE = "union", # Mode to handle reads overlapping more than one feature
# 		TYPE = "exon", # Feature type (3rd column in GFF file) to be used, all features of other type are ignored
# 		ID_ATTRIBUTE = "gene_id", # GFF attribute to be used as feature ID
# 		ORDER = "pos", # For paired-end data, the alignment have to be sorted either by read name or by alignment position.
# 		STRANDED = "no" # Whether the data is from a strand-specific assay
# 	message: "Quantifying read counts from BAM file {input.BAM} with HTSeq"
# 	shell:
# 		"""
# 		{HTSEQ} -m {params.MODE} -t {params.TYPE} -r {params.ORDER} \
# 		-i {params.ID_ATTRIBUTE} -f {params.FORMAT} \
# 		-s {params.STRANDED} {input.BAM} {input.GTF} > {output.COUNTS}
# 		"""

rule xx_featurecounts:
	input:
		BAM = SORTED_BAM_AL_DIR + "{sample}_RNA_XX_HISAT2_GRCh38_sortedbycoord.bam",
		GTF = XY_GTF
	output:
		COUNTS = FEATURECOUNTS_DIR + "{sample}_RNA_XX_HISAT2_transcript_featurecounts.tsv"
	params:
		THREADS = 5
	message: "Quantifying read counts from BAM file {input.BAM} with Subread featureCounts"
	shell:
		"""
		{FEATURECOUNTS} -T {params.THREADS} --primary -p -s 0 -t exon -g gene_id -a {input.GTF} -o {output.COUNTS} {input.BAM}
		"""

rule xy_featurecounts:
	input:
		BAM = SORTED_BAM_AL_DIR + "{sample}_RNA_XY_HISAT2_GRCh38_sortedbycoord.bam",
		GTF = XY_GTF
	output:
		COUNTS = FEATURECOUNTS_DIR + "{sample}_RNA_XY_HISAT2_transcript_featurecounts.tsv"
	params:
		THREADS = 5
	message: "Quantifying read counts from BAM file {input.BAM} with Subread featureCounts"
	shell:
		"""
		{FEATURECOUNTS} -T {params.THREADS} --primary -p -s 0 -t exon -g gene_id -a {input.GTF} -o {output.COUNTS} {input.BAM}
		"""