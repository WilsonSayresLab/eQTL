# Workflow for eQTL detection using PEER.

configfile: "tcga_lihc_testing.config.json"

# Tools

# Directories
RESULT_DIR = ""

rule all:
	input:

rule peer:
	input:
		GT_MATRIX = "",
		EXP_MATRIX = "",
		COV_MATRIX = ""
	output:
		PEER_RESULT = ""
	params:
	message: "Running eQTL analysis with PEER."
	script:
		"scripts/PEER.R"

rule matrix_eqtl:
	input:
		GT_MATRIX = "",
		EXP_MATRIX = "",
		COV_MATRIX = "",
		SNP_LOC = "",
		GENE_LOC = "",
		CIS_OUTPUT_NAME = "",
		TRANS_OUTPUT_NAME = ""
	output:
		MATRIXEQTL_RESULT = ""
	params:
	message: "Running eQTL analysis with Matrix eQTL."
	script:
		"scripts/MatrixeQTL.R"
