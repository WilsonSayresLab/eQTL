# Workflow for PCA, DE, and Go enrichment analysis with DEseq2,  identifying
# differentially expressed genes using edgeR, and runnin gene set enrichment
# analysis with GSEA.

configfile: "tcga_lihc_testing.config.json"

# Tools
GSEA = "gsea" # Path to GSEA.jar file

# Directories


# Data files
HTSEQ_COUNTS = "" # HTSeq count file for DESeq
COUNT_MATRIX = ""
# Use HUGO gene symbols across all GSEA data files
EXP_DATA = "" # res, gct, pcl, or txt. Contains features (genes or probes), samples, and an expression value for each feature in each sample.
PHEN_DATA = "" # cls. Contains phenotype labels and associates each sample with a phenotype.
GENE_SETS = "" # gmx or gmt. Contains one or more gene sets. For each gene set, gives the gene set name and list of features (genes or probes) in that gene set.
# CHIP_FILE = "" # Chip. Lists each probe on a DNA chip and its matching HUGO gene symbol. Optional for the gene set enrichment analysis.
GSETS = "" # Gene set names, comma delimited

rule all:
	input:

rule deseq:
	input:
		HTSEQ_COUNTS = HTSEQ_COUNTS
	output:
	params:
	message:
	script:
		"scripts/DESeq.R"

rule edger:
	input:
		COUNT_MATRIX = COUNT_MATRIX
	output:
		EDGER_RESULT = ""
	params:
	message: "Running differential expression analysis with edgeR."
	script:
		"scripts/edgeR.R"

rule gsea_gene_set_enrichment:
	input:
		EXP_DATA = EXP_DATA
		PHEN_DATA = PHEN_DATA
		GENE_SETS = GENE_SETS
		OUTPUT_DIR = OUTPUT_DIR
	output:
	params:
		xmx = "512m",
		collapse = "false" # If gene identifiers in the data set match those in the gene sets, set collapse to false. If dataset uses HG_U133A probe identifiers and gene sets use gene symbols, set to true.
		mode = "Max_probe"
		nperm = 1000
		rpt_label = "my_analysis"
		num = 100
		plot_top_x = 20
		set_max = 500
		set_min = 15
	message:
	shell:
		"""
		java -cp {GSEA} -Xmx{params.xmx} xtools.gsea.Gsea \
		-collapse {params.collapse} -mode {params.mode} -norm meandiv \
		-nperm {params.nperm} -permute phenotype \
		-rnd_type no_balance -scoring_scheme weighted \
		-rpt_label {params.rpt_label} -metric Signal2Noise -sort real \
		-order descending -include_only_symbols true -make_sets true \
		-median false -num {params.num} -plot_top_x {params.plot_top_x} \
		-rnd_seed timestamp -save_rnd_lists false -set_max {params.set_max} \
		-set_min {params.set_min} -zip_report false -gui false \
		-res {input.EXP_DATA} -cls {input.PHEN_DATA} \
		-gmx {input.GENE_SETS} -out {input.OUTPUT_DIR}
		"""

rule gsea_leading_edge:
	# The leading-edge subset in a gene set are those genes that appear in the
	# ranked list at or before the point at which the running sum reaches its
	# maximum deviation from zero. The leading-edge subset can be interpreted
	# as the core that accounts for the gene set’s enrichment signal. After
	# running the gene set enrichment analysis, you use the leading edge
	# analysis to examine the genes that are in the leading-edge subsets of
	# the enriched gene sets. A gene that is in many of the leading-edge
	# subsets is more likely to be of interest than a gene that is in only a
	# few of the leading-edge subsets.
	input:
		GSEA_RESULT_DIR = "" # Path to GSEA result directory
		GSETS = "" # Gene set names, comma delimited
	output:
	params:
		xmx = "512m",
	message:
	shell:
		"""
		java -cp {GSEA} -Xmx{params.xmx} xtools.gsea.LeadingEdgeTool \
		-dir {input.RESULT_DIR} -gsets {input.GSETS}
		"""
