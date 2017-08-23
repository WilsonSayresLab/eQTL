# Workflow for extracting FASTQ read files from BAM alignment files with
# XYAlign --STRIP_READS.

configfile: "EGA_v1.3.config.json"

# Tools
XYALIGN = "/home/hnatri/XYalign-noY_support/xyalign/xyalign.py" # Path to xyalign.py

# Reference genome files: XX with Y chromosome masked, XY with both X and Y
XY_REF = config["XY_hg19_ref_path"]
XY_REF_NAME = config["XY_hg19_ref_prefix"]
MASK_BED = "/home/hnatri/eQTL/QC/mask.bed"

# Directories
BAM_DIR = "/mnt/storage/euro/EGAD00001000131/"
RESULT_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/EGA_French/WES/xyalign_strip_reads/" # path to result directory
FASTQ_DIR = "/mnt/storage/CANCER_DOWNLOADS/PROCESSED/EGA_French/WES/stripped_fastqs/"

SAMPLES = config["EGAD1096_samples"]

rule all:
	input:
		expand(RESULT_DIR + "{sample}_WES_XYAlign_strip_reads/logfiles/{sample}_xyalign.log", RESULT_DIR=RESULT_DIR, sample=SAMPLES),
		#expand(FASTQ_DIR + "{sample}_WES_stripped_fq_1.fastq", FASTQ_DIR=FASTQ_DIR, sample=SAMPLES),
		#expand(FASTQ_DIR + "{sample}_WES_stripped_fq_2.fastq", FASTQ_DIR=FASTQ_DIR, sample=SAMPLES),

rule strip_reads:
	input:
		BAM = lambda wildcards: BAM_DIR + config[wildcards.sample]["Original_BAM"][0]
	output:
		LOG = RESULT_DIR + "{sample}_WES_XYAlign_strip_reads/logfiles/{sample}_xyalign.log"
	params:
		DIR = RESULT_DIR + "{sample}_WES_XYAlign_strip_reads/",
		SAMPLE_ID = "{sample}",
		cpus="8",
		xmx="4g",
		compression="3"
	message: "Stripping reads from {input.BAM} with XYAlign --STRIP_READS"
	shell:
		"set +u && source activate xyalign_environment && set -u && python {XYALIGN} --STRIP_READS --ref null --bam {input.BAM} --cpus {params.cpus} --xmx {params.xmx} --sample_id {params.SAMPLE_ID} --output_dir {params.DIR} --chromosomes ALL"

#TODO: rename and move stripped fastqs

#rule rename_fastqs:

#rule move_fastqs:
#   input:
#   output:
#   message:
#   shell:
#	   """
#	   mv STRIPPED_FQ FQ_DIR
#	   """
