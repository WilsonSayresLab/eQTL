# Workflow for calling germline SNP and indel variants using multiple variant
# callers.
# For germline SNP and indel variants:
#   Samtools mpileup and VarScan2,
#   GATK UnifiedGenotyper,
#   and Freebayes.
# For somatic SNP variants:
#   Samtools mpileup and Varscan2,
#   MuTect,
#   and SomaticSniper.
# For somatic indel variants:
#   Samtools mpileup and Varscan2,
#   XYZ,
#   and ABC.
# Variants called with at least two of the callers will be accepted.

# TODO:
# Add indel callers: PINDEL, DINDEL (only for Illumina data), Platypus (also
# calls SNPs, somatic and germline), Somatopus (somatic SNPs and short indels)

configfile: "tcga_lihc_testing.config.json"

# Tools
SAMTOOLS = "samtools"
VARSCAN2 = "varscan2" # Path to VarScan2 .jar
GATK = "gatk" # Path to GATK .jar
FREEBAYES = "freebayes"
MUTECT = "mutect" # Path to MuTect .jar
SOMATIC_SNIPER = "somaticsniper"

# Reference file indexed with Samtools faidx
REF = config["GRCh38_ref_path"]

# Database and interval iles
DBSNP = config["dbsnp_vcf"]
COSMIC = config["cosmic_vcf"]
INTERVAL_BED = config["interval_bed"]

# Directories
BAM_DIR = config["tcga_lihc_wgs_realignment_dir"]
VCF_DIR = config["tcga_lihc_vcf_dir"]

# Samples
NORMAL = config["tcga_lihc_wgs_normal_samples"]
TUMOR = config["tcga_lihc_wgs_tumor_samples"]

# Position sorted BAM alignment files
BAMS = config["tcga_lihc_wes_samples_realigned_bams"]

rule all:
    input:

rule mpileup_normal:
    input:
    output:
    params:
        min_MQ=20,
        min_BQ=50
    message: "Generating an mpileup file for normal samples."
    shell:
        """
        {SAMTOOLS} mpileup --fasta-ref {REF}
        --bam-list {bam_list} --min-MQ {params.min_MQ}} --min-BQ
        {params.min_BQ} -C50 -d10000
        {bam_files_list} -o {mpileup_file_path}
        """

rule mpileup_tumor:
    input:
    output:
    params:
        min_MQ=20,
        min_BQ=50
    message: "Generating an mpileup file for tumor samples."
    shell:
        """
        {SAMTOOLS} mpileup --fasta-ref {REF}
        --bam-list {bam_list} --min-MQ {params.min_MQ}} --min-BQ
        {params.min_BQ} -C50 -d10000
        {bam_files_list} -o {mpileup_file_path}
        """

rule varscan2_germline:
    input:
        MPILEUP =
    output:
    params:
        min_coverage=8, # Minimum read depth at a position to make a call
        min_reads2=2, # Minimum supporting reads at a position to call variants
        min_avg_qual=15, # Minimum base quality at a position to count a read
        mun_var_freq=0.01, # Minimum variant allele frequency threshold
        min_freq_for_hom=0.75, # Minimum frequency to call homozygote
        p_value=99e-02, # Default p-value threshold for calling variants
        strand_filter=1, # Ignore variants with >90% support on one strand
        output_vcf=1, # If set to 1, outputs in VCF format
        variants=1 # Report only variant (SNP/indel) positions
    message: "Calling germline SNP and indel variants from a Samtools mpileup \
    file with VarScan2 mpileup2cns."
    shell:
        """
        java -jar VarScan.jar mpileup2cns [mpileup file] OPTIONS
        {mpileup}
        """

rule gatk_unifiedgenotyper:
    input:
    output:
    params:
    message: "Calling germline SNP and indel variants with GATK \
    UnifiedGenotyper."
    shell:
        """
        """

rule freebayes:
    input:
    output:
    params:
    message: "Galling germline SNP and indel variants with Freebayes."
    shell:
        """
        """

rule varscan2_somatic:
    input:
    output:
    params:
        output_snp:"", # Output file for SNP calls [default: output.snp]
        output_indel="", # Output file for indel calls [default: output.indel]
        min_coverage=8, # Minimum coverage in normal and tumor to call variant
        min_coverage_normal=8, # Minimum coverage in normal to call somatic
        min_coverage_tumor=8, # Minimum coverage in tumor to call somatic
        min_var_freq=0.10, # Minimum variant frequency to call a heterozygote
        min_freq_for_hom=0.75, # Minimum frequency to call homozygote
        normal_purity=1.00, # Estimated purity (non-tumor content) of normal sample
        tumor_purity=1.00, # Estimated purity (tumor content) of tumor sample
        p_value=0.99, # P-value threshold to call a heterozygote
        somatic_p_value=0.05, # P-value threshold to call a somatic site
        strand_filter=1, # If set to 1, removes variants with >90% strand bias
        validation=0 # If set to 1, outputs all compared positions even if non-variant
    message: "Calling somatic SNP and indel variants from a Samtools mpileup \
    file with VarScan2 mpileup2cns."
    shell:
        """
        java -jar VarScan.jar somatic normal.pileup tumor.pileup output.basename
        """

rule mutect:
    input:
    output:
        VCF="",
        CALL_STATS="",
        COVERAGE_WIG=""
    params:
    message: "Calling somatic SNP variants with MuTect.""
    shell:
        """
        java -Xmx2g -jar {MUTECT}
        --analysis_type MuTect
        --reference_sequence {REF}
        --cosmic {COSMIC_VCF}
        --dbsnp {DBSNP_VCF}
        --intervals {INTERVAL_BED}
        --input_file:normal {input.NORMAL_BAM}
        --input_file:tumor {input.TUMOR_BAM}
        --out {output.CALL_STATS}}
        --coverage_file {output.COVERAGE_WIG}
        --vcf
        """

rule somatic_sniper:
    input:
    output:
    params:
    message: "Calling somatic SNP variants with SomaticSniper."
    shell:
        """
        """
