# Workflow for extracting FASTQ read files from BAM alignment files with
# XYAlign --STRIP_READS.

configfile: "tcga_lihc_testing.config.json"

# Tools
XYALIGN="/home/hnatri/XYalign-noY_support/xyalign/xyalign.py"

# Directories
BAM_DIR=config["tcga_lihc_wgs_bam_dir"]
FQ_DIR=config["tcga_lihc_wgs_stripped_fq_dir"]

# Samples
SAMPLES = config["tcga_lihc_wgs_samples"]

rule all:
	input:
		expand("{BAM_DIR}{sample}_wgs_Illumina.bam", BAM_DIR=BAM_DIR, sample=config["tcga_lihc_wgs_samples"]),

rule strip_reads:
    input:
        BAM="{BAM_DIR}{sample}.bam"
    output:
    params:
        cpus="4",
        xmx="4g",
        compression="3"
    message: "Stripping reads from {BAM_DIR}{BAM}"
    shell:
        """
        python {XYALIGN} --STRIP_READS --ref null --bam {input.bam}
        --cpus {params.cpus} --xmx {params.xmx} --sample_id {sample}
        --output_dir {FQ_DIR} --chromosomes ALL
        --fastq-compression {params.compression}
        """

#TODO:
#rule move_fastqs:
#   input:
#   output:
#   message:
#   shell:
#       """
#       mv STRIPPED_FQ FQ_DIR
#       """
