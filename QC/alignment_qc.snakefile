# Workflow for yielding alignment quality control statistics.

configfile: "tcga_lihc_testing.config.json"

# Tools
PICARD = "" # Path to Picard.jar file

# Reference genome files
XX_REF = config["XX_GRCh38_ref_path"]
XY_REF = config["XY_GRCh38_ref_path"]
REF = config["XY_GRCh38_ref_path"]
REF_NAME = config["XY_GRCh38_ref_prefix"]

# Directories
QC_DIR = config["test_qc_dir"]

# Alignment files
XX_ALIGNMENTS = []
XY_ALIGNMENTS = []
ALL_ALIGNMENTS = []

rule all:
	input:
		(REF_DIR + "{REF_NAME}.dict", REF_DIR=REF_DIR, REF_NAME=REF_NAME),
		(REF_DIR + "{interval_name}.interval_list", REF_DIR=REF_DIR, interval_name=config["test_interval_bed_name"]),
		expand(QC_DIR + "{alignment}_PicardQYM.txt", QC_DIR=QC_DIR, alignment=ALL_ALIGNMENTS),
		expand(QC_DIR + "{alignment}_PicardASM.txt", QC_DIR=QC_DIR, alignment=ALL_ALIGNMENTS),
		expand(QC_DIR + "{alignment}_PicardISM.txt", QC_DIR=QC_DIR, alignment=ALL_ALIGNMENTS),
		expand(QC_DIR + "{alignment}_PicardHSM.txt", QC_DIR=QC_DIR, alignment=ALL_ALIGNMENTS),
		expand(QC_DIR + "{alignment}_PicardWGSM.txt", QC_DIR=QC_DIR, alignment=ALL_ALIGNMENTS)

rule create_reference_sequence_dictionary:
	input:
		REF = XX_REF
	output:
		REF_DICT = REF_DIR + "{REF_NAME}.dict"
	params:
	message: "Creating a reference sequence dictionary for {input.REF}"
	shell:
		"""
		java -jar {PICARD} CreateSequenceDictionary \
      	R={input.REF} \
      	O={output.REF_DICT}
		"""

rule create_picard_interval_list:
	input:
		BED = config["test_interval_bed"],
		REF_DICT = config["test_ref_dict"]
	output:
		INTERVAL_LIST = REF_DIR + "{interval_name}.interval_list"
	params:
	message: "Creating a Picard Interval List from a BED file {input.BED} and \
	reference sequence dictionary {input.REF_DICT}"
	shell:
		"""
		java -jar {PICARD} BedToIntervalList \
      	I={input.BED} \
      	O={output.INTERVAL_LIST} \
      	SD={input.REF_DICT}
		"""

rule collect_alignment_summary_metrics:
	input:
		BAM = AL_DIR + "{alignment}_sorted_wrgs_mdup.bam",
		REF = REF
	output:
	params:
	message:
	shell:
		"""
		java -Dsamjdk.buffer_size={params.buffer_size} -Xmx{params.xmx} \
		-XX:GCTimeLimit={params.gc_timelimit} \
        -XX:GCHeapFreeLimit={params.gc_heap_free_limit} \
		-jar {PICARD} CollectMultipleMetrics \
		I={input.BAM} \
		O={output.RESULT_FILE} \
		R={input.REF} AS=true PROGRAM=null \
		PROGRAM=CollectAlignmentSummaryMetrics FILE_EXTENSION=.txt
		"""

rule collect_quality_yield_metrics:
	input:
		BAM = AL_DIR + "{alignment}_sorted_wrgs_mdup.bam",
		REF = REF
	output:
		RESULT_FILE = QC_DIR + "{alignment}_PicardQYM"
	params:
		buffer_size = 131072,
		xmx = "4096m",
		gc_timelimit = 50,
		gc_heap_free_limit = 10
	message: "Collecting Picard's Quality Yield Metrics for {input.BAM}"
	shell:
		"""
		java -Dsamjdk.buffer_size={params.buffer_size} -Xmx{params.xmx} \
		-XX:GCTimeLimit={params.gc_timelimit} \
        -XX:GCHeapFreeLimit={params.gc_heap_free_limit} \
		-jar {PICARD} CollectMultipleMetrics \
        CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
        I={input.BAM} \
		O={output.RESULT_FILE} \
		R={input.REF} AS=true PROGRAM=null \
        PROGRAM=CollectQualityYieldMetrics FILE_EXTENSION=.txt
		"""

rule collect_insert_size_metrics:
	input:
	output:
	params:
		buffer_size = 131072,
		xmx = "4096m",
		gc_timelimit = 50,
		gc_heap_free_limit = 10,
		minimum_pct=0.0
	message: "Collecting Picard's Insert Size Metrics for {input.BAM}"
	shell:
		"""
		java -Dsamjdk.buffer_size={params.buffer_size} -Xmx{params.xmx} \
		-XX:GCTimeLimit={params.gc_timelimit} \
        -XX:GCHeapFreeLimit={params.gc_heap_free_limit} \
		-jar {PICARD} CollectInsertSizeMetrics \
        VALIDATION_STRINGENCY=SILENT \
		INPUT={input.BAM} \
        HISTOGRAM_FILE={output.HISTOGRAM_FILE} \
		OUTPUT={output.RESULT_FILE} \
        MINIMUM_PCT={params.minimum_pct}
		"""

rule collect_hs_metrics:
	input:
		BAM = "",
		INTERVAL_LIST = ""
		REF = REF
	output:
		RESULT_FILE = "",
		PER_TARGET_HS_METRICS = ""
	params:
		buffer_size = 131072,
		xmx = "4096m",
		gc_timelimit = 50,
		gc_heap_free_limit = 10,
		near_distance = 200
	message:
	shell:
		"""
		java -Dsamjdk.buffer_size={params.buffer_size} -Xmx{params.xmx} \
		-XX:GCTimeLimit={params.gc_timelimit} \
        -XX:GCHeapFreeLimit={params.gc_heap_free_limit} \
		-jar {PICARD} CollectHsMetrics \
        BAIT_INTERVALS={input.INTERVAL_LIST} \
		TARGET_INTERVALS={inout.INTERVAL_LIST} \
		NEAR_DISTANCE={params.near_distance} \
		INPUT={input.BAM} \
        OUTPUT={output.RESULT_FILE} \
		PER_TARGET_COVERAGE={output.PER_TARGET_HS_METRICS} \
        REFERENCE_SEQUENCE={input.REF}
		"""

rule collect_wgs_metrics:
	input:
		BAM = "",
		REF = REF
	output:
		RESULT_FILE =
	params:
		buffer_size = 131072,
		xmx = "4096m",
		gc_timelimit = 50,
		gc_heap_free_limit = 10,
		BQ = 20,
		MAQ = 20,
		coverage_cap = 250
	message:
	shell:
		"""
		java -Dsamjdk.buffer_size={params.buffer_size} -Xmx{params.xmx} \
		-XX:GCTimeLimit={params.gc_timelimit} \
        -XX:GCHeapFreeLimit={params.gc_heap_free_limit} \
		-jar {PICARD} CollectWgsMetrics \
        MINIMUM_MAPPING_QUALITY={params.MAQ} \
		MINIMUM_BASE_QUALITY={params.BQ} \
        COVERAGE_CAP={params.coverage_cap} \
		VALIDATION_STRINGENCY=SILENT R={input.REF} \
        I={input.BAM} OUTPUT={output.RESULT_FILE} &> \
        picard_wgsmetrics_stdout_stderr.txt
		"""
