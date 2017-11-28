# Workflow for FastQC, MultiQC, and trimming of FASTQ read files.

configfile: "GTEx_Skin_v1.7.config.json"

# Tools
fastqc_path = "fastqc"
multiqc_path = "multiqc"
trimmomatic_path = "trimmomatic"

# Directory variables
fastq_directory = "/mnt/storage/public/dbgap-8834/skin_not_sun_exposed/"
adapter_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Skin/RNAseq/adapters/"
fastqc_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Skin/RNAseq/fastqc/"
trimmed_fastqs = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Skin/RNAseq/trimmed_fastqs/"
trimmed_fastqc_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Skin/RNAseq/fastqc/trimmed/"

rule all:
	input:
		expand(fastqc_directory + "{sample}_fq1_fastqc.html", fastqc_directory = fastqc_directory, sample=config["Skin_Not_Sun_Exposed_RNA"]),
		expand(fastqc_directory + "{sample}_fq2_fastqc.html", fastqc_directory = fastqc_directory, sample=config["Skin_Not_Sun_Exposed_RNA"]),
		(fastqc_directory + "multiqc_report.html"),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_1.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["Skin_Not_Sun_Exposed_RNA"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_1.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["Skin_Not_Sun_Exposed_RNA"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_2.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["Skin_Not_Sun_Exposed_RNA"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_2.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["Skin_Not_Sun_Exposed_RNA"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic.log", trimmed_fastqs=trimmed_fastqs, sample=config["Skin_Not_Sun_Exposed_RNA"]),
		expand(trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_1_fastqc.html", sample=config["Skin_Not_Sun_Exposed_RNA"]),
		expand(trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_2_fastqc.html", sample=config["Skin_Not_Sun_Exposed_RNA"]),
		(trimmed_fastqc_directory + "multiqc_report.html")

rule fastqc_analysis:
	input:
		fq1 = lambda wildcards: fastq_directory + wildcards.sample + "/" + wildcards.sample + "_fixed_1.fastq",
		fq2 = lambda wildcards: fastq_directory + wildcards.sample + "/" + wildcards.sample + "_fixed_2.fastq",
	output:
		fq1_zip = fastqc_directory + "{sample}_fq1_fastqc.zip",
		fq1_html = fastqc_directory + "{sample}_fq1_fastqc.html",
		fq2_zip = fastqc_directory + "{sample}_fq2_fastqc.zip",
		fq2_html = fastqc_directory + "{sample}_fq2_fastqc.html"
	params:
		fastqc_dir = fastqc_directory,
		fq1_prefix = lambda wildcards: wildcards.sample + "_fixed_1",
		fq2_prefix = lambda wildcards: wildcards.sample + "_fixed_2"
	shell:
		"""
		fastqc -o {params.fastqc_dir} {input.fq1};
		fastqc -o {params.fastqc_dir} {input.fq2};
		mv {params.fastqc_dir}{params.fq1_prefix}_fastqc.html {output.fq1_html};
		mv {params.fastqc_dir}{params.fq1_prefix}_fastqc.zip {output.fq1_zip};
		mv {params.fastqc_dir}{params.fq2_prefix}_fastqc.html {output.fq2_html};
		mv {params.fastqc_dir}{params.fq2_prefix}_fastqc.zip {output.fq2_zip}
		"""

rule multiqc:
	input:
	output:
		fastqc_directory + "multiqc_report.html"
	message: "Running MultiQC for FastQC reports located in {params.fastqc_dir}"
	params:
		fastqc_dir = fastqc_directory,
		output_dir = fastqc_directory
	shell:
		"""
		multiqc {params.fastqc_dir}*_fastqc.zip --outdir {params.output_dir} --interactive --verbose
		"""

rule trimmomatic:
	input:
		fq1 = lambda wildcards: fastq_directory + wildcards.sample + "/" + wildcards.sample + "_fixed_1.fastq",
		fq2 = lambda wildcards: fastq_directory + wildcards.sample + "/" + wildcards.sample + "_fixed_2.fastq",
		ADAPTER_FASTA = "/home/dlevy2/eQTL/Skin/QC/adapter_sequences.fa"
	output:
		paired_1 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		unpaired_1 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_1.fastq",
		paired_2 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_2.fastq",
		unpaired_2 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_2.fastq",
		logfile = trimmed_fastqs + "{sample}_trimmomatic.log"
	params:
		threads = 4,
		seed_mismatches = 2,
		palindrome_clip_threshold = 30,
		simple_clip_threshold = 10,
		leading = 3,
		trailing = 3,
		winsize = 4,
		winqual = 30,
		minlen = 50
	shell:
		"trimmomatic PE -threads {params.threads} -trimlog {output.logfile} {input.fq1} {input.fq2} {output.paired_1} {output.unpaired_1} {output.paired_2} {output.unpaired_2} ILLUMINACLIP:{input.ADAPTER_FASTA}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} LEADING:{params.leading} TRAILING:{params.trailing} SLIDINGWINDOW:{params.winsize}:{params.winqual} MINLEN:{params.minlen}"


rule fastqc_analysis_trimmomatic_trimmed_paired:
	input:
		fq1 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		fq2 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_2.fastq"
	output:
		html1 = trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_1_fastqc.html",
		html2 = trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_2_fastqc.html"
	params:
		fastqc_dir = trimmed_fastqc_directory
	shell:
		"""
		fastqc -o {params.fastqc_dir} {input.fq1} {input.fq2}
		"""

rule multiqc_trimmed_paired:
	input:
	output:
		trimmed_fastqc_directory + "multiqc_report.html",
	message: "Running MultiQC for post-trimming FastQC reports located in {params.fastqc_dir}"
	params:
		fastqc_dir = trimmed_fastqc_directory,
		output_dir = trimmed_fastqc_directory
	shell:
		"""
		multiqc {params.fastqc_dir}*trimmomatic*_fastqc.zip --outdir {params.output_dir} --interactive --verbose
		"""
