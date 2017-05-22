# Workflow for identifying differentially expressed genes using edgeR.

configfile: "tcga_lihc_testing.config.json"

# Tools


# Directories


rule all:
	input:

rule edger:
	input:
		COUNT_MATRIX = ""
	output:
		EDGER_RESULT = ""
	params:
	message: "Running differential expression analysis with edgeR."
	script:
		"scripts/edgeR.R"

rule gsea:
	input:
	output:
	params:
	message:
	shell:
		"""
		"""
