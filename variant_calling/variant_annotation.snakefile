# Workflow for annotating variants with SnpEff (dbSNP) and SnpSift (dbNSFP).

configfile: "tcga_lihc_testing.config.json"

# Tools
SNPEFF = "snpeff" # Path to SnpEff.jar file
SNPSIFT = "snpsift" # Path to SnpSift.jar file

# Database files
DBNSFP = config["dbnsfp"] # Path to dbNSFP.txt.gz
DBNSFP_INDEX=config["dbnsfp"] # Path to dbNSFP.txt.gz.tbi

# Directories
VCF_DIR = config["tcga_lihc_vcf_dir"]

# VCFs
VCFS = config["tcha_lihc_vcfs"]

# VCF prefixes

# SnpEff/dbSNP annotated VCFS

# SnpSift/dbNSFP annotated VCFs

rule snpeff_annotate:
    input:
    output:
    params:
		xmx="4g",
		mode="ann", # Annotate variants / calculate effects. Default: ann (no command or 'ann')
		genome_version="GRCh37.75"
    message: "Annotating {input.VCF} with SnpEff/dbSNP."
    shell:
        """
		java -jar {SNPEFF} -Xmx{params.xmx} {params.genome_version}
		{input.VCF} > {output.annotated_VCF}
        """

rule snpsift_annotate:
    input:
    output:
    params:
    message: "Annotation {VCF} with SnpSift/dbNSFP."
    shell:
        """
		java -jar {SNPSIFT} dbnsfp -v -db {DBNSFP} {iput.VCF} > {output.annotated_VCF}
        """
