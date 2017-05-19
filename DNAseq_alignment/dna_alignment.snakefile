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

XX_ALIGNMENTS = []
XY_ALIGNMENT = []

rule all:
    input:
		expand(AL_DIR + "{xx_alignment}.sam", AL_DIR=AL_DIR, alignment=XX_ALIGNMENTS),
		expand(AL_DIR + "{xy_alignment}.sam", AL_DIR=AL_DIR, alignment=XY_ALIGNMENTS),
		expand(AL_DIR + "{xx_alignment}_stats.txt", AL_DIR=AL_DIR, xx_alignment=XX_ALIGNMENTS),
		expand(AL_DIR + "{xy_alignment}_stats.txt" AL_DIR=AL_DIR, xy_alignment=XY_ALIGNMENTS),
		expand(AL_DIR + "{alignment}.bam", AL_DIR=AL_DIR, alignment=ALL_ALIGNMENTS),
		expand(AL_DIR + "{alignment}_sorted.bam", AL_DIR=AL_DIR, alignment=ALL_ALIGNMENTS),
		expand(AL_DIR + "{alignment}_sorted.bam.bai", AL_DIR=AL_DIR, alignment=ALL_ALIGNMENTS))
		
rule xx_align_paired_reads:
    # Aligning paired FASTQs of XX samples to the XX reference genome. Writing
    # alignment statistics into a .txt file.
    input:
        R1 = lambda wildcards: FQ_DIR + config["xx_test_fqs"][wildcards.sample][0],
        R2 = lambda wildcards: FQ_DIR + config["xx_test_fqs"][wildcards.sample][1]
    output:
        SAM = "{xx_alignment}.sam",
        STATS = "{xx_alignment}_stats.txt"
    params:
    message: "Aligning {wildcards.sample} to {XX_REF_NAME} using Bowtie2."
    shell:
        """
        {BOWTIE2} -q -x {XX_BOWTIE_INDEX} -1 {input.R1} -2 {input.2} -S {output.SAM} 2> \
        {output.STATS}
        """

rule xy_align_paired_reads:
    # Aligning paired FASTQs of XY samples to the XY reference genome. Writing
    # alignment statistics into a .txt file.
    input:
		R1 = lambda wildcards: FQ_DIR + config["xy_test_fqs"][wildcards.sample][0],
		R2 = lambda wildcards: FQ_DIR + config["xy_test_fqs"][wildcards.sample][1],
    output:
        SAM = AL_DIR + "{xy_alignment}.sam",
        STATS = AL_DIR + "{xy_alignment}_stats.txt"
    params:
    message: "Aligning {wildcards.sample} to {XY_REF_NAME} with Bowtie2."
    shell:
        """
        {BOWTIE2} -q -x {XY_BOWTIE_INDEX} -1 {input.R1} -2 {input.R2} -S {output.SAM} 2> \
        {output.STATS}
        """

rule convert_sam_to_bam:
    input:
        SAM = AL_DIR + "{alignment}.sam"
    output:
        BAM = AL_DIR + "{alignment}.bam"
    params:
    message: "Converting SAM file {input.SAM} to BAM with Samtools."
    shell:
        """
		{SAMTOOLS} view -bTS {input.REF} {input.SAM} > {output.BAM}
        """

rule sort_bam:
    input:
        UNSORTED_BAM = AL_DIR + "{alignment}.bam"
    output:
        SORTED_BAM = AL_DIR + "{alignment}_sorted.bam"
    params:
		cores=1
    message: "Sorting BAM file {input.UNSORTED_BAM} with Samtools."
    shell:
        """
		{SAMTOOLS} sort {input.UNSORTED_BAM} \
        -@ {params.cores} -o {output.SORTED_BAM}
        """

# TODO:
# Mark duplicates, for instance with samtools rmdup
# rule mark_duplicates:

rule index_bam:
	input:
		BAM = AL_DIR + "{alignment}_sorted.bam"
	output:
		BAI = AL_DIR + "{alignment}_sorted.bam.bai"
	message: "Indexing BAM file {BAM} with Samtools"
	params:
	shell:
		"""
		{SAMTOOLS} index {input.BAM}
		"""
