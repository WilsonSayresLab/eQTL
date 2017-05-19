# Workflow for eQTL detection using PEER.

configfile: "tcga_lihc_testing.config.json"

# Tools
PEER = ""

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
	params:
	message:
	script:
		"scripts/PEER.R"
