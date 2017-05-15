import os

configfile: "eqtl.config.json"

# Tools
STAR = ""
SAMTOOLS = "samtools"
HTSEQ = ""

# Reference genome files
REF = config["GRCh38_ref_path"]
REF_NAME = config["GRCh38_ref_prefix"]
<<<<<<< HEAD
GTF = config["gtf_path"]
STAR_REF_DIR = config["star_ref_dir"]
=======
>>>>>>> master

# Directories
FQ_DIR = "/mnt/storage/CANCER_DOWNLOADS/LIHCFILES/data/rna/sequences_data/"
AL_DIR = "/mnt/storage/CANCER_DOWNLOADS/LIHCFILES/RNAseq_realignment/"
<<<<<<< HEAD
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
=======

FQ_DIR = config["tcga_rnaseq_fastq_dir"] # path to directory with fq files

WC = glob_wildcards(os.path.join(FQ_DIR, "{sample}_{pair}.fastq"))

SAMPLES = set(WC.sample) # unique sample names
PAIR1, PAIR2 = set(WC.pair) # set of paired end specifiers

rule all:
     input:
         expand("{S}.bam", S=SAMPLES)

# TODO:  RNA alignment whit STAR and HTSeq
#		 DNA alignment with bowtie2 and Samtools

rule map:
	input:
		R1="{fq_dir}/{{S}}_{pair}.fastq".format(fq_dir=tcga_rnaseq_fastq_dir, pair=PAIR1),
	    R2="{fq_dir}/{{S}}_{pair}.fastq".format(fq_dir=tcga_rnaseq_fastq_dir, pair=PAIR2),
		gtf=gtf_path
	output: "{al_dir}/{{S}}}_{ref_name}.bam".format(al_dir=tcga_rnaseq_alignment_dir, ref_name=REF_NAME)
	params:
		star_index = star_index_path,
		quantmode = quantmode,
>>>>>>> master
		threads=24
	message: "Mapping reads to {REF_NAME} with STAR. Threads: {params.threads}"
	shell:
		"""
<<<<<<< HEAD
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
=======
		{STAR} --runThreadN {params.threads} --runMode alignReads --genomeDir {input.star_index} --outFileNamePrefix {wildcards.sample}_ --readFilesIn {input.R1} {input.R2} --outSAMtype BAM SortedByCoordinate --quantMode {params.quant_param} --sjdbGTFfile {input.gtf};
		mv {wildcards.sample}_Aligned.out.bam {output}
		"""

rule index_bam:
     input:
         BAM="{al_dir}/{{S}}_.fastq".format(fq_dir=FQ_DIR, pair=PAIR1)
     output:
         "{S}.bam.bai"
     threads: 4
     log: "logs/bam.log"
     message: "Indexing BAM with samtools"
	 params:
	 	count_output = "_raw_gene_counts.tsv"
     shell:
		"""
        {samtools_path} index
		{output_folder}/{output_file}Aligned.sortedByCoord.out.bam;
        {read_count_method_path} -f bam -r pos -s no {output_folder}/{output_file}Aligned.sortedByCoord.out.bam {annotation_file_path} >
        {output_folder}/{output_file}{count_output} &> {log}
         """

rule prepare_reference_hg38:
	input:
		hg38_ref_path
	output:
		fai = hg38_ref_path + ".fai",
		amb = hg38_ref_path + ".amb",
		dict = hg38_ref_prefix + ".dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		# faidx
		shell("{params.samtools} faidx {input}")
		# .dict
		shell("{params.samtools} dict -o {output.dict} {input}")
		# bwa
		shell("{params.bwa} index {input}")
>>>>>>> master
