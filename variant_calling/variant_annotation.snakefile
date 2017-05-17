# Workflow for annotating variants with SnpEff (dbSNP) and SnpSift (dbNSFP).

configfile: "tcga_lihc_testing.config.json"

# Tools
SNPEFF = "snpeff"

# Database files
DBSNP = config["dbsnp_vcf"]
DBNSFP = config["dbnsfp_vcf"]

# Directories
VCF_DIR = config["tcga_lihc_vcf_dir"]

# VCFs
VCFS = config["tcha_lihc_vcfs"]

rule snpeff_annotate:
    input:
    output:
    params:
    message: "Annotating {VCF} "
    shell:
        """
        """

rule snpsift_annotate:
    input:
    output:
    params:
    message:
    shell:
        """
        """
