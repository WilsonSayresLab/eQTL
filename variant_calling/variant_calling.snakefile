# Workflow for calling germline SNP and indel variants using multiple variant
# callers.
# For germline SNP and indel variants:
#   Samtools mpileup and VarScan2,
#   GATK HaplotypeCaller,
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
REF_FA_INDEX = config[""] # Reference FASTA index file (.fa.fai)

# Database and interval files
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
		BAM_LIST = "" # File with paths of indexed BAM files, one line per file. These files should be better tagged with read groups; if not, one BAM per sample.
	output:
	params:
		min_MQ=20,
		min_BQ=50
	message: "Generating an mpileup file for normal samples."
	shell:
		"""
		{SAMTOOLS} mpileup --fasta-ref {REF}
		--bam-list {input.BAM_LIST} --min-MQ {params.min_MQ}}
		-min-BQ {params.min_BQ} -C50 -d10000
		{input.BAM_FILE_LIST} -o {output.MPILEUP}
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
		java -jar {VARSCAN2} mpileup2cns {input.MPILEUP}
		--min-coverage {params.min_coverage}
		--min_reads2 {params.min_reads}
		--min_avg_qual {params.min_avg_qual}
		--mun_var_freq {params.min_var_freq}
		--min_freq_for_hom {params.min_freq_for_hom}
		--p_value {params.p_value}
		--strand_filter {params.strand_filter}
		--output_vcf {params.output_vcf}
		--variants {params.variants}

		"""

rule gatk_haplotypecaller:
	input:
	output:
	params:
		heterozygosity=0.001, # Heterozygosity value used to compute prior likelihoods for any locus
		heterozygosity_stdev=0.01, # Standard deviation of eterozygosity for SNP and indel calling
		indel_heterozygosity=1.25E-4, # Heterozygosity for indel calling
		maxReadsInRegionPerSample=10000, # Maximum reads in an active region
		min_base_quality_score=10, # Minimum base quality required to consider a base for calling
		minReadsPerAlignmentStart=10, # Minimum number of reads sharing the same alignment start for each genomic location in an active region
		sample_ploidy=2, # Ploidy per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy)
		standard_min_confidence_threshold_for_calling=10.0, # The minimum phred-scaled confidence threshold at which variants should be called
		max_alternate_alleles=6, # Maximum number of alternate alleles to genotype
		max_genotype_count=1024, # Maximum number of genotypes to consider at any site
		output_mode="EMIT_VARIANTS_ONLY" # Which type of calls we should output
	message: "Calling germline SNP and indel variants with GATK \
	HaplotypeCaller."
	shell:
		"""
		java -jar {GATK}}
		-R {REF}
		-T HaplotypeCaller
		-I sample1.bam [-I sample2.bam ...]
		--heterozygosity {params.heterozygosity}
		--heterozygosity_stdev {params.heterozygosity_stdev}
		--indel_heterozygosity {params.indel_heterozygosity}
		--maxReadsInRegionPerSample {params.maxReadsInRegionPerSample}
		--min_base_quality_score {params.min_base_quality_score}
		--minReadsPerAlignmentStart {params.minReadsPerAlignmentStart}
		--sample_ploidy {params.sample_ploidy}
		--standard_min_confidence_threshold_for_calling {params.standard_min_confidence_threshold_for_calling}
		--max_alternate_alleles {params.max_alternate_alleles}
		--max_genotype_count {params.max_genotype_count}
		--output_mode {params.output_mode}
		--dbsnp {DBSNP}
		[-L targets.interval_list]
		-o {OUTPUT}
		"""

rule freebayes:
	input: BAMs, REF
	output: VCF
	params:
		cores=36,
		ploidy=2,
		min_mapping_quality=1, # Exclude alignments from analysis if they have a mapping quality less than Q
		min_base_quality=0, # Exclude alleles from analysis if their supporting base quality is less than Q
		min_alternate_fraction=0.2, # Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position
		min_alternate_count=2, # Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position
		min_alternate_qsum=0, # Require at least this sum of quality of observations supporting an alternate allele within a single individual in order to evaluate the position
		use_best_n_alleles=4, # Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores. Set to 0 to use all.
		min_repeat_entropy=1, # The haplotype which is called has Shannon entropy less than --min-repeat-entropy, which is off by default but can be set to ~1 for optimal genotyping of indels in lower-complexity sequence
		phred=20 # Remove sites with estimated probability of not being polymorphic less than phred 20 (aka 0.01), or probability of polymorphism > 0.99
	message: "Galling germline SNP and indel variants with Freebayes."
	shell:
		"""
		freebayes-parallel <(fasta_generate_regions.py {REF_FA_INDEX} 100000)
		{params.cores}
		-f {input.REF}
		{input.BAMs}
		-p {params.ploidy}
		--min-mapping-quality {params.min_mapping_quality}
		--min-base-quality {params.min_base_quality}
		--min-alternate-fraction {params.min_alternate_fraction}
		--min-alternate-count {params.min_alternate_count}
		--min-alternate-qsum {params.min_alternate_qsum}
		--use-best-n-alleles {params.use_best_n_alleles}
		--min-repeat-entropy {params.min_repeat_entropy}
		| vcffilter -f "QUAL > {input.phred}" > {OUTPUT}
		"""

rule varscan2_somatic:
	input:
	output:
	params:
		output_snp:"", # Output file for SNP calls [default: output.snp]
		output_indel = "", # Output file for indel calls [default: output.indel]
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
		NORMAL_BAM = {},
		TUMOR_BAM = {},
		REF = REF
	output:
		VCF = "{}_mutect.vcf",
		CALL_STATS = "{}_mutect_callstats.txt",
		COVERAGE_WIG = "{}_mutect_coverage.wig"
	params:
		initial_tumor_lod=4.0, # Initial LOD threshold for calling tumor variant
		tumor_lod=6.3, # LOD threshold for calling tumor variant
		fraction_contamination=0, # Estimate of fraction (0-1) of physical contamination with other unrelated samples
		minimum_mutation_cell_fraction=0.00, # Minimum fraction of cells which are presumed to have a mutation, used to handle non-clonality and contamination
		normal_lod=2.2, # LOD threshold for calling normal non-germline
		normal_artifact_lod=1.0, # LOD threshold for calling normal non-variant
		strand_artifact_lod=2.0, # LOD threshold for calling strand bias
		minimum_normal_allele_fraction=0.00, # Minimum allele fraction to be considered in normal, useful for normal sample contaminated with tumor
		tumor_f_pretest=0.005, # For computational efficiency, reject sites with allelic fraction below this threshold
		min_qscore=10, # Threshold for minimum base quality score, default=5
		gap_events_threshold=3, # How many gapped events (ins/del) are allowed in proximity to this candidate
		heavily_clipped_read_fraction=0.30, # If this fraction or more of the bases in a read are soft/hard clipped, do not use this read for mutation calling
		max_alt_alleles_in_normal_count=2, # Threshold for maximum alternate allele counts in normal
		max_alt_allele_in_normal_fraction=0.03, # Threshold for maximum alternate allele fraction in normal

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
		--min_qscore
		--vcf
		"""

rule somatic_sniper:
	input:
		TUMOR_BAM = {},
		NORMAL_BAM = {},
		REF = REF
	output:
		VCF = "{}_somaticsniper.vcf"
	params:
		q=0, # Filtering reads with mapping quality less than INT
		Q=40, # Filtering somatic snv output with somatic quality less than INT
		L="-L", # Do not report LOH variants as determined by genotypes
		G="-G" # Do not report Gain of Referene variants as determined by genotypes
		p="-p" # Disable priors in the somatic calculation. Increases sensitivity for solid tumors.
		J="-J" # Use prior probabilities accounting for the somatic mutation rate
		s=0.01, # Prior probability of a somatic mutation (implies -J)
		T=0.850000, # Theta in maq consensus calling model (for -c/-g)
		N=2, # Number of haplotypes in the sample (for -c/-g)
		r=0.001000, # Prior of a diﬀerence between two haplotypes (for -c/-g) [0.001000]
		F="vcf" # Output format (vcf or classic)
	message: "Calling somatic SNP variants with SomaticSniper."
	shell:
		"""
		bam-somaticsniper -q {params.q} -Q {params.Q} {params.L} {params.G}
		{params.p} {params.J} -s {params.s} -T {params.T} -N {params.N}
		-r {params.r} -f {input.REF} {input.TUMOR_BAM} {input.NORMAL_BAM}
		-F {params.F} {output.VCF}
		"""
