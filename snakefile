import os

configfile: "eqtl.config.json"

# Tools
STAR = ""
SAMTOOLS = "samtools"
HTSEQ = ""

# Reference genome files
REF = config["GRCh38_ref_path"]
REF_NAME = config["GRCh38_ref_prefix"]
GTF = config["gtf_path"]
STAR_REF_DIR = config["star_ref_dir"]

# Directories
FQ_DIR = "/mnt/storage/CANCER_DOWNLOADS/LIHCFILES/data/rna/sequences_data/"
AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/LIHCFILES/RNAseq_realignment/"
FQ_DIR = config["tcga_rnaseq_fastq_dir"] # path to directory with fq files

# Getting unique fastq file and sample names
WC = glob_wildcards(os.path.join(FQ_DIR, "{sample}_{pair}.fastq"))
SAMPLES = set(WC.sample) # unique sample names
PAIR1, PAIR2 = set(WC.pair) # set of paired end specifiers

#rule all:
     #input:
         #expand("{S}.bam", S=SAMPLES)

rule star_index:
	input:
	output:
	params:
	shell:
		"""
	    "{STAR} --runThreadN {params.threads} --runMode genomeGenerate "
        "--genomeDir {star_ref_dir} --genomeFastaFiles {ref} --sjdbGTFfile {gtf}; "
		"""

rule map:
	input:
		R1="{fq_dir}/{wildcards.sample}_{pair}.fastq".format(pair=PAIR1),
	    R2="{fq_dir}/{wildcards.sample}_{pair}.fastq".format(pair=PAIR2),
		star_index=star_index,
		gtf=GTF,
		name_prefix:"{wildcards.sample}_"
	output: "{al_dir}/{wildcards.sample}_Aligned.sortedByCoord.out.bam"
	params:
		threads=24
	message: "Mapping reads to {REF_NAME} with STAR. Threads: {params.threads}"
	shell:
		"""
		{STAR} --runThreadN {params.threads} --runMode alignReads
		--genomeDir {input.star_index} --outFileNamePrefix {input.name_prefix}
		--readFilesIn {input.R1} {input.R2} --outSAMtype BAM SortedByCoordinate
		--sjdbGTFfile {input.gtf};
		"""


rule index_bam:
     input:
         BAM="{al_dir}/{wildcards.sample}_Aligned.sortedByCoord.out.bam"
     output: "{al_dir}/{wildcards.sample}_Aligned.sortedByCoord.out.bam.bai"
     message: "Indexing BAM with samtools"
	 params:
     shell:
		"""
        {SAMTOOLS} index {input.BAM}
        """

rule quantify:
	input:
		BAM:"{al_dir}/{wildcards.sample}_Aligned.sortedByCoord.out.bam",
		gtf=GTF
	output:
		count_output = "{al_dir}/{wildcards.sample}_raw_gene_counts.tsv"
	threads:
	message: "Counting reads with HTSeq"
	params:
	shell:
		"""
		{HTSEQ} -f bam -r pos -s no {input.BAM} {input.gtf}"
		"""
