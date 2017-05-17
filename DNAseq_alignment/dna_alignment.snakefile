# Workflow for aligning DNA sequence data to a sex specific reference with
# Bowtie2. Indexing BAM files with Samtools.

configfile: "tcga_lihc_testing.config.json"

# Tools
BOWTIE2 = "bowtie2"
SAMTOOLS = "samtools"

# Reference genome files: XX with Y chromosome masked, XY with both X and Y
XX_REF = config["xx_GRCh38_ref_path"]
XX_REF_NAME = config["xx_GRCh38_ref_prefix"]
XX_BOWTIE_INDEX = config["xx_bowtie_index"]
XY_REF = config["xy_GRCh38_ref_path"]
XY_REF_NAME = config["xy_GRCh38_ref_prefix"]
XY_BOWTIE_INDEX = config["xy_bowtie_index"]

# Directories
FQ_DIR=config["tcga_lihc_wgs_stripped_fq_dir"]
BAM_DIR=config["tcga_lihc_wgs_realignment_dir"]

# Samples
XX_SAMPLES = config["tcga_lihc_rnaseq_females"]
XY_SAMPLES = config["tcga_lihc_rnaseq_males"]
ALL_SAMPLES = config["tcga_lihc_rnaseq_samples"]

rule all:
    input:

rule xx_align_paired_reads:
    # Aligning paired FASTQs of XX samples to the XX reference genome. Writing
    # alignment statistics into a .txt file.
    input:
        R1="",
        R2=""
    output:
        SAM="",
        STATS=""
    params:
    message: "Aligning {SAMPLE} to {XX_REF_NAME} using Bowtie2."
    shell:
        """
        {BOWTIE2} -q -x {XX_BOWTIE_INDEX} -1 {R1} -2 {R2} -S {SAM} 2>
        {STATS}
        """

rule xy_align_paired_reads:
    # Aligning paired FASTQs of XY samples to the XY reference genome. Writing
    # alignment statistics into a .txt file.
    input:
    output:
    params:
    message: "Aligning {SAMPLE} to {XY_REF_NAME} with Bowtie2."
    shell:
        """
        """

rule convert_sam_to_bam:
    input:
        SAM=""
    output:
        BAM=""
    params:
    message: "Converting SAM file {SAM} to BAM with Samtools."
    shell:
        """
        """

rule sort_bam:
    input:
        BAM=""
    output:
        SORTED_BAM=""
    params:
    message: "Sorting BAM file {BAM} with Samtools."
    shell:
        """
        """

rule index_bam:
	input:
		BAM=""
	output: ".bam.bai"
	message: "Indexing BAM file {BAM} with Samtools"
	params:
	shell:
		"""
		{SAMTOOLS} index {input.BAM}
		"""
