# Workflow for FastQC, MultiQC, adapter discovery, and trimming of FASTQ read files.

configfile: "GTEx_v1.2.config.json"

# Directory variables
fastq_directory = "/mnt/storage/public/dbgap-8834/liver/"
adapter_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Liver/RNAseq/adapters/"
fastqc_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Liver/RNAseq/fastqc/"
trimmed_fastqs = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Liver/RNAseq/trimmed_fastqs/"
trimmed_fastqc_directory = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/GTEx_Liver/RNAseq/fastqc/trimmed/"

# Tools
#xyalign_path = "/home/awolf4/XYalign/xyalign/xyalign.py"
bbmerge_sh_path = "bbmerge.sh"
bbduk_sh_path = "bbduk.sh"
fastqc_path = "fastqc"
trimmomatic_path = "trimmomatic"

rule all:
	input:
		# expand(fastqc_directory + "{sample}_fq1_fastqc.html", fastqc_directory = fastqc_directory, sample=config["Liver_RNA"]),
		# expand(fastqc_directory + "{sample}_fq2_fastqc.html", fastqc_directory = fastqc_directory, sample=config["Liver_RNA"]),
		# expand(adapter_directory + "{sample}.adapters.fa", adapter_directory = adapter_directory, sample=config["Liver_RNA"]),
		# (fastqc_directory + "multiqc_report.html")
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_1.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["Liver_RNA"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_1.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["Liver_RNA"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_2.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["Liver_RNA"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_2.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["Liver_RNA"]),
		expand(trimmed_fastqs + "{sample}_trimmomatic.log", trimmed_fastqs=trimmed_fastqs, sample=config["Liver_RNA"]),
		# expand(trimmed_fastqs + "{sample}_bbduk_trimmed_read1.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["Liver_RNA"]),
		# expand(trimmed_fastqs + "{sample}_bbduk_trimmed_read2.fastq", trimmed_fastqs=trimmed_fastqs, sample=config["Liver_RNA"])
		# expand(trimmed_fastqc_directory + "{sample}_bbduk_trimmed_read1_fastqc.html", sample=config["Liver_RNA"]),
		# expand(trimmed_fastqc_directory + "{sample}_bbduk_trimmed_read2_fastqc.html", sample=config["Liver_RNA"]),
		expand(trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_1_fastqc.html", sample=config["Liver_RNA"]),
		expand(trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_2_fastqc.html", sample=config["Liver_RNA"]),
		# (trimmed_fastqc_directory + "bbduk_multiqc_report.html"),
		# (trimmed_fastqc_directory + "bbduk_multiqc_data"),
		(trimmed_fastqc_directory + "trimmomatic_multiqc_report.html"),
		(trimmed_fastqc_directory + "trimmomatic_multiqc_data")

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
		fastqc = fastqc_path,
		fastqc_dir = fastqc_directory,
		fq1_prefix = lambda wildcards: wildcards.sample + "_fixed_1",
		fq2_prefix = lambda wildcards: wildcards.sample + "_fixed_2"
	shell:
		"""
		{params.fastqc} -o {params.fastqc_dir} {input.fq1};
		{params.fastqc} -o {params.fastqc_dir} {input.fq2};
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
		ADAPTER_FASTA = "/home/hnatri/eQTL/GTEx_liver/QC/adapter_sequences.fa"
	output:
		paired_1 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		unpaired_1 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_1.fastq",
		paired_2 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_2.fastq",
		unpaired_2 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_unpaired_2.fastq",
		logfile = trimmed_fastqs + "{sample}_trimmomatic.log"
	params:
		trimmomatic = trimmomatic_path,
		threads = 4,
		seed_mismatches = 2,
		palindrome_clip_threshold = 30,
		simple_clip_threshold = 10,
		leading = 3,
		trailing = 3,
		winsize = 4,
		winqual = 25,
		minlen = 50
	shell:
		"{params.trimmomatic} PE -threads {params.threads} -phred33 -trimlog {output.logfile} {input.fq1} {input.fq2} {output.paired_1} {output.unpaired_1} {output.paired_2} {output.unpaired_2} ILLUMINACLIP:{input.ADAPTER_FASTA}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} LEADING:{params.leading} TRAILING:{params.trailing} SLIDINGWINDOW:{params.winsize}:{params.winqual} MINLEN:{params.minlen}"

# rule adapter_discovery:
# 	input:
# 		fq1 = lambda wildcards: fastq_directory + config[wildcards.sample]["RNAseq_fastq_dir"][0] + "/" + config[wildcards.sample]["RNAseq_fastq_1"][0],
# 		fq2 = lambda wildcards: fastq_directory + config[wildcards.sample]["RNAseq_fastq_dir"][0] + "/" + config[wildcards.sample]["RNAseq_fastq_2"][0]
# 	output:
# 		adapter_directory + "{sample}.adapters.fa"
# 	params:
# 		bbmerge_sh = bbmerge_sh_path
# 	shell:
# 		"{params.bbmerge_sh} in1={input.fq1} in2={input.fq2} outa={output} reads=1m"
#
# rule trim_adapters_paired_bbduk:
# 	input:
# 		fq1 = lambda wildcards: fastq_directory + config[wildcards.sample]["RNAseq_fastq_dir"][0] + "/" + config[wildcards.sample]["RNAseq_fastq_1"][0],
# 		fq2 = lambda wildcards: fastq_directory + config[wildcards.sample]["RNAseq_fastq_dir"][0] + "/" + config[wildcards.sample]["RNAseq_fastq_2"][0],
# 		ADAPTER_FASTA = "/home/hnatri/eQTL/QC/adapter_sequences.fa"
# 	output:
# 		out_fq1 = trimmed_fastqs + "{sample}_bbduk_trimmed_read1.fastq",
# 		out_fq2 = trimmed_fastqs + "{sample}_bbduk_trimmed_read2.fastq"
# 	params:
# 		bbduk_sh = bbduk_sh_path
# 	shell:
# 		"{params.bbduk_sh} -Xmx1g in1={input.fq1} in2={input.fq2} out1={output.out_fq1} out2={output.out_fq2} ref={input.ADAPTER_FASTA} ktrim=r k=21 mink=11 hdist=2 tbo tbe qtrim=rl trimq=15 minlen=36"

# rule fastqc_analysis_bbduk_trimmed_paired:
# 	input:
# 		fq1 = trimmed_fastqs + "{sample}_bbduk_trimmed_read1.fastq",
# 		fq2 = trimmed_fastqs + "{sample}_bbduk_trimmed_read2.fastq"
# 	output:
# 		html1 = trimmed_fastqc_directory + "{sample}_bbduk_trimmed_read1_fastqc.html",
# 		html2 = trimmed_fastqc_directory + "{sample}_bbduk_trimmed_read2_fastqc.html"
# 	params:
# 		fastqc = fastqc_path,
# 		fastqc_dir = trimmed_fastqc_directory
# 	shell:
# 		"""
# 		{params.fastqc} -o {params.fastqc_dir} {input.fq1} {input.fq2};
# 		"""

rule fastqc_analysis_trimmomatic_trimmed_paired:
	input:
		fq1 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_1.fastq",
		fq2 = trimmed_fastqs + "{sample}_trimmomatic_trimmed_paired_2.fastq"
	output:
		html1 = trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_1_fastqc.html",
		html2 = trimmed_fastqc_directory + "{sample}_trimmomatic_trimmed_paired_2_fastqc.html"
	params:
		fastqc = fastqc_path,
		fastqc_dir = trimmed_fastqc_directory
	shell:
		"""
		{params.fastqc} -o {params.fastqc_dir} {input.fq1} {input.fq2}
		"""

rule multiqc_trimmed_paired:
	input:
	output:
		trimmed_fastqc_directory + "trimmomatic_multiqc_report.html",
		trimmed_fastqc_directory + "trimmomatic_multiqc_data",
	message: "Running MultiQC for post-trimming FastQC reports located in {params.fastqc_dir}"
	params:
		fastqc_dir = trimmed_fastqc_directory,
		output_dir = trimmed_fastqc_directory
	shell:
		"""
		multiqc {params.fastqc_dir}*trimmomatic*_fastqc.zip --outdir {params.output_dir} --interactive --verbose;
		mv {params.output_dir}multiqc_report.html {params.output_dir}trimmomatic_multiqc_report.html;
		mv /mnt/storage/CANCER_DOWNLOADS/PROCESSED/RNAseq/fastqc/trimmed/multiqc_data /mnt/storage/CANCER_DOWNLOADS/PROCESSED/RNAseq/fastqc/trimmed/trimmomatic_multiqc_data;
		"""
