# eQTL

This repository contains the code used in the eQTL project:

Author: Heini Natri (heini.natri@asu.edu)

Snakemake workflows:
  - Alignment
	- RNA-seq alignment with STAR
	- RNA-seq quantification with HTSeq
	- DNA-seq alignment with bowtie2
  - Variant calling and annotation
	- Samtools mpileup and VarScan2
	- GATK UnifiedGenotyper
	- Freebayes
	- MuTect
	- SomaticSniper
	- PINDEL/DINDEL/Platypus/Somatopus
	- SnpEff: dbSNP
	- SnpSift: dbNSFP
  - eQTL detection
	- PEER	
  - Downstream analyses
	- Differential expression analysis with edgeR
