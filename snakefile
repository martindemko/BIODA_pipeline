import math
import subprocess
from snakemake.utils import R
from snakemake.utils import report

#from pytools.persistent_dict import PersistentDict
#storage = PersistentDict("mystorage")

########################################
# DEFINITION OF DEFAULT INPUT SETTINGS
#
DEF_GENOME = "default_genome.fa"
DEF_ANALYSIS_TYPE  = "classic"
DEF_DATA_TYPE = "paired"
DEF_ORGANISM_TYPE = "mammals"
DEF_TARGET = "report_all.html"
DEF_INPUT_DIR = "."
DEF_OUTPUT_DIR = "."
DEF_ANNOTATION = "default_annotation.gff3"
DEF_GENOME_INDEX = ""
DEF_SAM_HEADER = "default_SAM_header"
DEF_SE_DATA = "data"
DEF_PE1_DATA = "data1"
DEF_PE2_DATA = "data2"
DEF_ADAPTERS = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/adapters_merge.fa"
DEF_STAR_PARAMS = ""
DEF_STAR_ALIGN_INTRON_MAX = "0"       # default value from STAR manual
DEF_STAR_ALIGN_MATES_GAP = "0"        # default value from STAR manual


########################################
# DEFINITION OF INPUT SETTINGS NAMES
#
CONF_GENOME_TERM = "genome"
CONF_ANALYSIS_TYPE_TERM = "analysis_type"
CONF_DATA_TYPE_TERM = "data_type"
CONF_ORGANISM_TYPE_TERM = "organism_type"
CONF_TARGET_TERM = "target"
CONF_OUTPUT_DIR_TERM = "output_dir"
CONF_INPUT_DIR_TERM = "input_dir"
CONF_ANNOTATION_TERM = "annotation"
CONF_GENOME_INDEX_TERM = "genome_index"
CONF_SAM_HEADER_TERM = "SAM_header"
CONF_SE_DATA_TERM = "data"
CONF_PE1_DATA_TERM = "data1"
CONF_PE2_DATA_TERM = "data2"
CONF_STAR_PARAMS_TERM = "STAR_params"
CONF_STAR_ALIGN_INTRON_MAX_TERM = "STAR_alignIntronMax"
CONF_STAR_ALIGN_MATES_GAP_TERM = "STAR_alignMatesGapMax"
CONF_ADAPTERS_TERM = "adapters"


###########################################
# DEFINITION OF INPUT SETTINGS VARIABLES
#
#.split() is used if it's neccesary to use more than one option
GENOME = config[CONF_GENOME_TERM] if (CONF_GENOME_TERM in config) else DEF_GENOME
ANALYSIS_TYPE  = config[CONF_ANALYSIS_TYPE_TERM].split() if (CONF_ANALYSIS_TYPE_TERM in config) else DEF_ANALYSIS_TYPE 
DATA_TYPE = config[CONF_DATA_TYPE_TERM].split() if (CONF_DATA_TYPE_TERM in config) else DEF_DATA_TYPE  
ORGANISM_TYPE = config[CONF_ORGANISM_TYPE_TERM] if (CONF_ORGANISM_TYPE_TERM in config) else DEF_ORGANISM_TYPE
TARGET = config[CONF_TARGET_TERM] if (CONF_TARGET_TERM in config) else DEF_TARGET
INPUT_DIR = config[CONF_INPUT_DIR_TERM] if (CONF_INPUT_DIR_TERM in config) else DEF_INPUT_DIR
OUTPUT_DIR = config[CONF_OUTPUT_DIR_TERM] if (CONF_OUTPUT_DIR_TERM in config) else DEF_OUTPUT_DIR
ANNOTATION = config[CONF_ANNOTATION_TERM] if (CONF_ANNOTATION_TERM in config) else DEF_ANNOTATION
GENOME_INDEX = config[CONF_GENOME_INDEX_TERM] if (CONF_GENOME_INDEX_TERM in config) else DEF_GENOME_INDEX
SAM_HEADER = config[CONF_SAM_HEADER_TERM] if (CONF_SAM_HEADER_TERM in config) else DEF_SAM_HEADER
DATA1 = config[CONF_PE1_DATA_TERM] if (CONF_PE1_DATA_TERM in config) else DEF_PE1_DATA
DATA2 = config[CONF_PE2_DATA_TERM] if (CONF_PE2_DATA_TERM in config) else DEF_PE2_DATA
ADAPTERS = config[CONF_ADAPTERS_TERM] if (CONF_ADAPTERS_TERM in config) else DEF_ADAPTERS
if (config[CONF_DATA_TYPE_TERM] == "single") :
  DATA = config[CONF_SE_DATA_TERM] if (CONF_SE_DATA_TERM in config) else DEF_SE_DATA
  SAMPLES = DATA.strip(".gz").strip(".fq|.fastq")
else:
  DATA = DATA1.split()+DATA2.split()
  SAMPLES = DATA1.strip(".gz").strip(".fq|.fastq").split()+DATA2.strip(".gz").strip(".fq|.fastq").split()

#if not os.path.exists(OUTPUT_DIR):
#  os.makedirs(OUTPUT_DIR)

# GENOME or ORGANISM dependent settings
#
if GENOME == "arabidopsis":
  GENOME_FILE = ""
  ORGANISM_TYPE = "plants"
  STAR_INTRON_MAX = "10000"   # 10k
  STAR_MATES_GAP = "10000"
elif GENOME == "human":
  GENOME_FILE = ""
  ORGANISM_TYPE = "mammals"
  STAR_INTRON_MAX = "100000"  # 100k
  STAR_MATES_GAP = "100000"
elif GENOME == "sifilis":
  GENOME_FILE = ""
  ORGANISM_TYPE = "bacteria"
  STAR_INTRON_MAX = "1"       # indicates no splicing
  STAR_MATES_GAP = "1000"
else:
  GENOME_FILE = GENOME
  ORGANISM_TYPE = "undefined"
  STAR_INTRON_MAX = DEF_STAR_ALIGN_INTRON_MAX
  STAR_MATES_GAP = DEF_STAR_ALIGN_MATES_GAP

STAR_INTRON_MAX = config[CONF_STAR_ALIGN_INTRON_MAX_TERM] if (CONF_STAR_ALIGN_INTRON_MAX_TERM in config) else STAR_INTRON_MAX
STAR_MATES_GAP = config[CONF_STAR_ALIGN_MATES_GAP_TERM] if (CONF_STAR_ALIGN_MATES_GAP_TERM in config) else STAR_MATES_GAP
STAR_PARAMS = config[CONF_STAR_PARAMS_TERM] if (CONF_STAR_PARAMS_TERM in config) else DEF_STAR_PARAMS
help = subprocess.Popen("grep -v '>' "+GENOME_FILE+" | wc -m",shell=True,stdout=subprocess.PIPE).communicate()[0]
STAR_GENOME_BASES_LOG = min(14,math.floor(math.log(float(int(help)),2)/2-1))

STRANDS = ["forward","reverse"]

####################################
# DEFINITION OF TOOLS
#
FASTQC = "fastqc"
ADAPTERS = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/adapters_merge.fa"
REAPER_SRC = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/reaper-15-065/src"
TRIMMOMATIC = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/trimmomatic-master/classes/trimmomatic.jar"
STAR = "STAR"
SAMTOOLS = "samtools"
RSEM_PATH = ""
BBMAP = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/bbmap"
FEATURE_COUNTS = "featureCounts"
UCSC_SCRIPTS = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools"
GTF_TO_BED12 = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/gtf2bed12.py"
PICARD = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/picard-tools-1.119"
PICARD_JAR = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/picard-tools-1.119/picard-1.119.jar"
PRESEQ = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/preseq_v2.0.2/preseq" #"/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/preseq_v2.0.1/preseq"
RSEQC = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/RSeQC-2.6.4/build/scripts-2.7"
DUPRADAR = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/dupRadar.R"

PHRED_FILTER = 5  # Trim the 3' end of read if four consequent bases have average base quality smaller than this value
LEN_FILTER = 20   # Filter sequences shorter than this value after quality trimming
RD_LENGTH = 75    # Trim all sequences to this maximal length - depends on experiment design and I like to trim them to specified length


####################################
# DEFINITION OF FINAL RULES
#
rule all:
  input:  "report_all.html"
  
  
rule final_report:
  input:  #expand("{dir}/1st_qc",dir=OUTPUT_DIR),
          #expand("{dir}/minion/{data}.minion.compare",data=SAMPLES,dir=OUTPUT_DIR) if config[CONF_DATA_TYPE_TERM] == "single" else expand("{dir}/minion/{data}.minion.compare",data=DATA2.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR),
          expand("{dir}/trimming",dir=OUTPUT_DIR),
          expand("{dir}/2nd_qc",dir=OUTPUT_DIR),
          expand("{dir}/aligned/Log.final.out",dir=OUTPUT_DIR),
          expand("{dir}/RSEM/RSEM_calc_expr.pdf",dir=OUTPUT_DIR),
          expand("{dir}/FeatureCounts/feature_counts.txt",dir=OUTPUT_DIR),
          expand("{dir}/mapped_QC/{data}.rRNA.intervalListBody.txt", dir=OUTPUT_DIR, data=ANNOTATION.rstrip(".gtf")),
          expand("{dir}/strandness.txt", dir=OUTPUT_DIR),
          expand("{dir}/preseq/Aligned.sortedByCoord.estimates.txt",dir=OUTPUT_DIR),
          expand("{dir}/mapped_QC/featureCounts.quantSeq.rev.biotype_counts.txt",dir=OUTPUT_DIR),
          expand("{dir}/mapped_QC/Aligned.sortedByCoord.dupRadar_dupMatrix.txt",dir=OUTPUT_DIR)
  output: "report_all.html"
  shell:  "touch {output}"


#rule step_zero_report:
#  input:  expand("{dir}/1st_qc/{data}_fastqc.html", dir=OUTPUT_DIR, data=SAMPLES),
#          expand("{dir}/minion/{data}.minion.compare",data=SAMPLES,dir=OUTPUT_DIR) if config[CONF_DATA_TYPE_TERM] == "single" else expand("{dir}/minion/{data}.minion.compare",data=DATA2.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR)
#  output: "report_step_zero.html"
#  shell:  "touch {output}"
  
  
####################################
# DEFINITION OF ZERO-STEP RULES (FIRST QC, ADAPTERS CHECK)
#
rule first_qc:
  input:  reads = expand("{dir}/{data}",data=DATA,dir=INPUT_DIR) if config[CONF_DATA_TYPE_TERM] == "single" else expand(["{dir}/{data1}","{dir}/{data2}"],data1=DATA1,data2=DATA2,dir=INPUT_DIR)
  output: html = expand("{dir}/1st_qc/{data}_fastqc.html", dir=OUTPUT_DIR, data=SAMPLES)
  log:    run = expand("{dir}/1st_qc/run_stats.log", dir=OUTPUT_DIR),
          report = expand("{dir}/report.html", dir=OUTPUT_DIR)
  params: extra = "--noextract --format fastq --nogroup",
          prefix = expand("{dir}/1st_qc", dir=OUTPUT_DIR)
  threads: 4
  run:  
          shell("""
          {FASTQC} -o {params.prefix}/ {params.extra} --threads {threads} {input.reads} > {log.run} 2>&1
          """)
        ### REPORT PART ###
          report("""
          =======================
          The title of the report
          =======================
  
          Write your report here, explaining your results. Don't fear to use math
          it will be rendered correctly in any browser using MathJAX,
          e.g. inline :math:`\sum_{{j \in E}} t_j \leq I`,
          or even properly separated:
  
          .. math::
  
              |cq_{{0ctrl}}^i - cq_{{nt}}^i| > 0.5
  
          Include your files using their keyword name and an underscore: reads_.
          Original html report from FastQC is here: html_.
  
          Access your global and local variables like within shell commands, e.g. INPUT_DIR={INPUT_DIR}.
          """, *log.report, **input, **output)
          
       
if config[CONF_DATA_TYPE_TERM] == "single":
  rule check_adapters:
    input:  reads = expand("{dir}/{data}",data=DATA,dir=INPUT_DIR)
    output: minion = expand("{dir}/minion/{data}.minion.fa",data=DATA.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR)
    log:    run = expand("{dir}/minion/{data}.log",data=DATA.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR),
            report = expand("{dir}/report.html", dir=OUTPUT_DIR)
    params: minion = "-show 3 -write-fasta"
    run:    
            shell("""
            {REAPER_SRC}/minion search-adapter -i {input.reads} {params.minion} {output.minion} > {log.run} 2>&1
            """)
          ### REPORT PART ###
            report("""
            =====================
            Minion part of report
            =====================
            
            Original output is minion_.
            """, *log.report, **input, **output)
else:
  rule check_adapters:
    input:  reads1 = expand("{dir}/{data}",data=DATA1,dir=INPUT_DIR),
            reads2 = expand("{dir}/{data}",data=DATA2,dir=INPUT_DIR)
    output: minion1 = expand("{dir}/minion/{data}.minion.fa",data=DATA1.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR),
            minion2 = expand("{dir}/minion/{data}.minion.fa",data=DATA2.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR)
    log:    run1 = expand("{dir}/minion/{data}.log",data=DATA1.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR),
            run2 = expand("{dir}/minion/{data}.log",data=DATA2.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR),
            report = expand("{dir}/report.html", dir=OUTPUT_DIR)
    threads:  1
    params: minion = "-show 3 -write-fasta"
    run:    
            shell("""
            {REAPER_SRC}/minion search-adapter -i {input.reads1} {params.minion} {output.minion1} > {log.run1} 2>&1
  	        {REAPER_SRC}/minion search-adapter -i {input.reads2} {params.minion} {output.minion2} > {log.run2} 2>&1
            """)
            ### REPORT PART ###
            report("""
            =====================
            Minion part of report
            =====================
            
            Original output is minion1_ and minion2_.
            """, *log.report, **input, **output)

       
if config[CONF_DATA_TYPE_TERM] == "single":
  rule check_adapters_swan:
    input:  minion = expand("{dir}/minion/{data}.minion.fa",data=DATA.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR)
    output: swan = expand("{dir}/minion/{data}.minion.compare",data=DATA.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR)
    shell:  """
  	        {REAPER_SRC}/swan -r {ADAPTERS} -q {input.minion} > {output.swan} 2>&1
            """
else:
  rule check_adapters_swan:
    input:  minion1 = expand("{dir}/minion/{data}.minion.fa",data=DATA1.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR),
            minion2 = expand("{dir}/minion/{data}.minion.fa",data=DATA2.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR)
    output: swan1 = expand("{dir}/minion/{data}.minion.compare",data=DATA1.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR),
            swan2 = expand("{dir}/minion/{data}.minion.compare",data=DATA2.strip(".gz").strip(".fq|.fastq"),dir=OUTPUT_DIR)
    threads:  1
    shell:  """
  	        {REAPER_SRC}/swan -r {ADAPTERS} -q {input.minion1} > {output.swan1} 2>&1
  	        {REAPER_SRC}/swan -r {ADAPTERS} -q {input.minion2} > {output.swan2} 2>&1
            """


#########################################################
# DEFINITION OF FIRST-STEP RULES (TRIMMING, SECOND QC)
#
if config[CONF_ANALYSIS_TYPE_TERM] == "quant":
  if config[CONF_DATA_TYPE_TERM] == "single":
    rule trimming_by_bbmap:
      input:  reads = expand("{dir}/{data}",data=DATA,dir=INPUT_DIR),
              adapters = expand("{data}",data=ADAPTERS)
      output: dir = expand("{dir}/trimming",dir=OUTPUT_DIR),
              clean = expand(["{dir}/trimming/{data}.clean.fastq.gz"],data=DATA.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR),
              waste = expand(["{dir}/trimming/{data}.trimmed.fastq.gz"],data=DATA.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR)
      log:    run = expand("{dir}/trimming/run_stats.log",dir=OUTPUT_DIR),
              stats = expand("{dir}/trimming/cont_stats.log",dir=OUTPUT_DIR),
              rpkm = expand("{dir}/trimming/rpkm_stats.log",dir=OUTPUT_DIR),
              kmers = expand("{dir}/trimming/kmers_stats.log",dir=OUTPUT_DIR)
      threads: 4
      params: extra = "literal=GGGGGGGGG,AAAAAAAAA k=13 useshortkmers=t mink=5 qtrim=rl trimq=10 minlength=20" # -Xmx10g for JVM is possible
      shell:  """
              {BBMAP}/bbduk.sh in={input.reads} out={output.clean} outm={output.waste} ref={input.adapters} stats={log.stats} rpkm={log.rpkm} dump={log.kmers} threads={threads} ktrim=r {params.extra} > {log.run} 2>&1
              """
  else:
    rule trimming_by_bbmap:
      input:  r1 = expand("{dir}/{data1}",data1=DATA1,dir=INPUT_DIR),
              r2 = expand("{dir}/{data2}",data2=DATA2,dir=INPUT_DIR),
              adapters = expand("{data}",data=ADAPTERS)
      output: dir = expand("{dir}/trimming",dir=OUTPUT_DIR),
              c1 = expand(["{dir}/trimming/{data1}.clean.fastq.gz"],data1=DATA1.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR),
              c2 = expand(["{dir}/trimming/{data2}.clean.fastq.gz"],data2=DATA2.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR),
              w1 = expand(["{dir}/trimming/{data1}.trimmed.fastq.gz"],data1=DATA1.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR),
              w2 = expand(["{dir}/trimming/{data2}.trimmed.fastq.gz"],data2=DATA2.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR),
              sing = expand(["{dir}/trimming/{data}.singleton.fastq.gz"],data=DATA1.strip(".fastq.gz|.fq.gz").strip("_1|_R1|_L1|1"),dir=OUTPUT_DIR)
      log:    run = expand("{dir}/trimming/run_stats.log",dir=OUTPUT_DIR),
              stats = expand("{dir}/trimming/cont_stats.log",dir=OUTPUT_DIR),
              rpkm = expand("{dir}/trimming/rpkm_stats.log",dir=OUTPUT_DIR),
              kmers = expand("{dir}/trimming/kmers_stats.log",dir=OUTPUT_DIR)
      threads: 4
      params: extra = "literal=GGGGGGGGG,AAAAAAAAA k=13 useshortkmers=t mink=5 qtrim=rl trimq=10 minlength=20"
      shell:  """
              {BBMAP}/bbduk.sh in={input.r1} in2={input.r2} out={output.c1} out2={output.c2} outm={output.w1} outm2={output.w2} outs={output.sing} ref={input.adapters} stats={log.stats} rpkm={log.rpkm} dump={log.kmers} threads={threads} ktrim=r {params.extra} > {log.run} 2>&1
              """
elif config[CONF_ANALYSIS_TYPE_TERM] == "classic":
  if config[CONF_DATA_TYPE_TERM] == "single":
    rule trimming_by_Trimmomatic:
      input:  reads = expand("{dir}/{data}",data=DATA,dir=INPUT_DIR)
      output: dir = expand("{dir}/trimming",dir=OUTPUT_DIR),
              clean = expand(["{dir}/trimming/{data}.clean.fastq.gz"],data=DATA.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR)
      log:    run = expand("{dir}/trimming/run_stats.log",dir=OUTPUT_DIR)
      threads: 4
      params: phred = "-phred33",
              leading = "3",
              trailing = "3",
              crop = expand("{par}",par=RD_LENGTH),
              minlen = expand("{par}",par=LEN_FILTER),
              slid_w_1 = "4",
              slid_w_2 = expand("{par}",par=PHRED_FILTER)
      wrapper:
              "file:/mnt/nfs/shared/999993-Bioda/scripts/martin/test/wrapper_trim_se"
#      shell:  """
#              java -jar {TRIMMOMATIC} SE -threads {threads} {params.phred} {input.reads} {output.clean} LEADING:{params.leading} TRAILING:{params.trailing} 
#              CROP:{params.crop} SLIDINGWINDOW:{params.slid_w_1}:{params.slid_w_2} MINLEN:{params.minlen} > {log.run} 2>&1
#              """
  else:
    rule trimming_by_Trimmomatic:
      input:  r1 = expand("{dir}/{data1}",data1=DATA1,dir=INPUT_DIR),
              r2 = expand("{dir}/{data2}",data2=DATA2,dir=INPUT_DIR)
      output: dir = expand("{dir}/trimming",dir=OUTPUT_DIR),
              r1p = expand(["{dir}/trimming/{data1}.clean.fastq.gz"],data1=DATA1.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR),
              r1u = expand(["{dir}/trimming/{data1}.unpaired.fastq.gz"],data1=DATA1.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR),
              r2p = expand(["{dir}/trimming/{data2}.clean.fastq.gz"],data2=DATA2.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR),
              r2u = expand(["{dir}/trimming/{data2}.unpaired.fastq.gz"],data2=DATA2.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR)
      log:    run = expand("{dir}/trimming/run_stats.log",dir=OUTPUT_DIR)
      threads: 4
      params: phred = "-phred33",
              leading = "3",
              trailing = "3",
              crop = expand("{par}",par=RD_LENGTH),
              minlen = expand("{par}",par=LEN_FILTER),
              slid_w_1 = "4",
              slid_w_2 = expand("{par}",par=PHRED_FILTER)
      shell:  """
              java -jar {TRIMMOMATIC} PE -threads {threads} {params.phred} {input.r1} {input.r2} {output.r1p} {output.r1u} {output.r2p} {output.r2u} LEADING:{params.leading} TRAILING:{params.trailing} CROP:{params.crop} SLIDINGWINDOW:{params.slid_w_1}:{params.slid_w_2} MINLEN:{params.minlen} > {log.run} 2>&1
              """
   
            
rule second_qc:
  input:  reads = expand("{dir}/trimming/{data}.clean.fastq.gz",data=DATA.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR) if config[CONF_DATA_TYPE_TERM] == "single" else expand(["{dir}/trimming/{data1}.clean.fastq.gz","{dir}/trimming/{data2}.clean.fastq.gz"],data1=DATA1.strip(".fastq.gz|.fq.gz"),data2=DATA2.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR)
  output: dir = "{OUTPUT_DIR}/2nd_qc"
  log:    run = "{OUTPUT_DIR}/2nd_qc/run_stats.log"
  threads: 4
  params: "--noextract --format fastq --nogroup"
  shell:  """
          {FASTQC} -o {output.dir}/ {params} --threads {threads} {input.reads} > {log.run} 2>&1
          """


#########################################################
# DEFINITION OF MAPPING RULES
#
if GENOME_INDEX == DEF_GENOME_INDEX :
  rule STAR_gen_index:
    input:  genome = expand("{data}",data=GENOME_FILE),
            ref = expand("{data}",data=ANNOTATION)
    output: dir = expand("{dir}/genome_index/{data}",dir=OUTPUT_DIR,data=GENOME_FILE.rstrip(".gz").rstrip(".fa|.fasta")),
            prefix = expand("{dir}/genome_index/",dir=OUTPUT_DIR)
    threads:  12
    params: extra = "",
            Nbases = expand("--genomeSAindexNbases {val}",val=STAR_GENOME_BASES_LOG)
    shell:  """
            mkdir -p {output.dir}
            {STAR} --runMode genomeGenerate --runThreadN {threads} --genomeDir {output.dir} --genomeFastaFiles {input.genome} --sjdbGTFfile {input.ref} --outFileNamePrefix {output.prefix} {params.Nbases} {params.extra}
            """


import os
rule STAR_alignment:
  input:  reads = expand("{dir}/trimming/{data}.clean.fastq.gz",data=DATA.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR) if config[CONF_DATA_TYPE_TERM] == "single" else expand(["{dir}/trimming/{data1}.clean.fastq.gz","{dir}/trimming/{data2}.clean.fastq.gz"],data1=DATA1.strip(".fastq.gz|.fq.gz"),data2=DATA2.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR),
          gen_index = expand("{data}",data=GENOME_INDEX) if GENOME_INDEX != DEF_GENOME_INDEX else expand("{dir}/genome_index/{data}",dir=OUTPUT_DIR,data=GENOME_FILE.rstrip(".gz").rstrip(".fa|.fasta")),
          ref = expand("{data}",data=ANNOTATION)
  output: final_log = expand("{dir}/aligned/Log.final.out",dir=OUTPUT_DIR),
          gen_bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR),
          trans_bam_tmp = temp(expand("{dir}/aligned/Aligned.toTranscriptome.out.bam",dir=OUTPUT_DIR)),
          trans_bam = expand("{dir}/aligned/Aligned.toTranscriptome.sortedByCoord.out.bam",dir=OUTPUT_DIR)
  threads: 12
  resources:  mem_mb=10000
  params: reads = expand("{dir}/trimming/{data}.clean.fastq.gz",data=DATA.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR) if config[CONF_DATA_TYPE_TERM] == "single" else expand("{dir}/trimming/{data1}.clean.fastq.gz {dir}/trimming/{data2}.clean.fastq.gz",data1=DATA1.strip(".fastq.gz|.fq.gz"),data2=DATA2.strip(".fastq.gz|.fq.gz"),dir=OUTPUT_DIR),
          dir = expand("{dir}/aligned/",dir=OUTPUT_DIR),
          read_command = "zcat",  #might be conditional to input
          sort_RAM = 10000000000, #could be also a config parameter
          sam_type = "--outSAMtype BAM SortedByCoordinate",
          quant = "--quantMode TranscriptomeSAM", # if config[CONF_ANALYSIS_TYPE_TERM] == "quant" else "",
          splicing = expand("--alignIntronMax {val}",val=STAR_INTRON_MAX),
          mates_gap = expand("--alignMatesGapMax {val}",val=STAR_MATES_GAP),
          extra = expand("{line}",line=STAR_PARAMS)
  shell:  """
          {STAR} --runMode alignReads --runThreadN {threads} --genomeDir {input.gen_index} --readFilesIn {params.reads} --readFilesCommand {params.read_command} --sjdbGTFfile {input.ref} --outFileNamePrefix {params.dir} {params.sam_type} --limitBAMsortRAM {params.sort_RAM} {params.splicing} {params.mates_gap} {params.quant} {params.extra}
          {SAMTOOLS} sort -@ {threads} {output.trans_bam_tmp} -o {output.trans_bam}
          """
          
rule STAR_index:
  input:  final_log = expand("{dir}/aligned/Log.final.out",dir=OUTPUT_DIR),
          gen_bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR),
          trans_bam = expand("{dir}/aligned/Aligned.toTranscriptome.sortedByCoord.out.bam",dir=OUTPUT_DIR)
  output: gen_bai = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam.bai",dir=OUTPUT_DIR),
          trans_bai = expand("{dir}/aligned/Aligned.toTranscriptome.sortedByCoord.out.bam.bai",dir=OUTPUT_DIR)
  threads:  12
  shell:  """
          {SAMTOOLS} index -@ {threads} {input.gen_bam}
          
          {SAMTOOLS} index -@ {threads} {input.trans_bam}
          """

## I want this to be optional
#rule SAM_to_BAM:
#  input:  "{OUTPUT_DIR}/aligned/{SAMPLE}.sam"
#  output: "{OUTPUT_DIR}/aligned/{SAMPLE}.bam"
#  threads:  4
#  shell:  """
#          {SAMTOOLS} view -@ {threads} -b {input} > {output}
#	        rm {input}
#          """


#######################################
# DEFINITION OF COUNTING RULES (RSEM)
#
rule RSEM_prep_ref:
  input:  ref = expand("{data}",data=ANNOTATION),
          genome = expand("{data}",data=GENOME_FILE)
  output: idx = expand("{dir}/RSEM/{data}.idx.fa",dir=OUTPUT_DIR,data=ANNOTATION.rstrip(".gz").rstrip(".gtf|.gff|.gff3"))
  log:    run = expand("{dir}/RSEM/RSEM_prep_ref.log",dir=OUTPUT_DIR)
  threads:  4
  params: rsem_ref = expand("{dir}/RSEM/{data}",dir=OUTPUT_DIR,data=ANNOTATION.rstrip(".gz").rstrip(".gtf|.gff|.gff3")),
          use_ref = "--gtf" if ANNOTATION.endswith(".gtf") else "--gff3"
  shell:  """
          {RSEM_PATH}rsem-prepare-reference --num-threads {threads} {params.use_ref} {input.ref} {input.genome} {params.rsem_ref} > {log.run} 2>&1
          """


rule RSEM_calc_expr:
  input:  bam = expand("{dir}/aligned/Aligned.toTranscriptome.sortedByCoord.out.bam",dir=OUTPUT_DIR),
          idx = expand("{dir}/RSEM/{data}.idx.fa",dir=OUTPUT_DIR,data=ANNOTATION.rstrip(".gz").rstrip(".gtf|.gff|.gff3"))
  output: dir = expand("{dir}/RSEM",dir=OUTPUT_DIR),
          help_bam = temp(expand("{dir}/RSEM/Aligned.toTranscriptome.converted_for_RSEM.bam",dir=OUTPUT_DIR)),
          pdf = expand("{dir}/RSEM/RSEM_calc_expr.pdf",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/RSEM/RSEM_calc_expr.log",dir=OUTPUT_DIR),
          convert = expand("{dir}/RSEM/RSEM_convert_input.log",dir=OUTPUT_DIR)
  threads: 4
  resources:  mem_mb=10000
  params: mem = 10000,
          paired = "" if config[CONF_DATA_TYPE_TERM] == "single" else "--paired-end",
          extra = "--estimate-rspd --calc-ci --no-bam-output --seed 12345 --forward-prob 0",
          bam = expand("{dir}/aligned/Aligned.toTranscriptome.sortedByCoord.out.bam",dir=OUTPUT_DIR),
          ref = expand("{dir}/RSEM/{data}",dir=OUTPUT_DIR,data=ANNOTATION.rstrip(".gz").rstrip(".gtf|.gff|.gff3")),
          help_bam = temp(expand("{dir}/RSEM/Aligned.toTranscriptome.converted_for_RSEM",dir=OUTPUT_DIR)),
          prefix = expand("{dir}/RSEM/RSEM_calc_expr",dir=OUTPUT_DIR)
  run:  
          shell("""
	        #cat <( {SAMTOOLS} view -H {input.bam} ) <( {SAMTOOLS} view -@ {threads} {input.bam} | awk '{{printf "%s", $0 " "; getline; print}}' | sort -S {params.mem}M -T ./ | tr ' ' '\\n' ) | {SAMTOOLS} view -@ {threads} -bS - > {output.help_bam}
	        {RSEM_PATH}convert-sam-for-rsem -p {threads} {input.bam} {params.help_bam} &> {log.convert}
          {RSEM_PATH}rsem-calculate-expression --alignments {params.paired} {params.extra} -p {threads} --ci-memory {params.mem} {output.help_bam} {params.ref} {params.prefix} >& {log.run}
	        {RSEM_PATH}rsem-plot-model {params.prefix} {output.pdf}
          """)


#rule RSEM_gen_data_matrix:


###############################################
# DEFINITION OF COUNTING RULES (FeatureCounts)
#
rule FeatureCounts:
  input:  ref = expand("{data}",data=ANNOTATION),
          bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR)
  output: expand("{dir}/FeatureCounts/feature_counts.txt",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/FeatureCounts/run_stats.log",dir=OUTPUT_DIR)
  threads:  4
  params: paired = "" if config[CONF_DATA_TYPE_TERM] == "single" else "-p"
          
  shell:  """
          {FEATURE_COUNTS} -t exon -g gene_id {params.paired} -T {threads} -a {input.ref} -o {output} {input.bam} > {log.run} 2>&1
          """

###############################################
# DEFINITION OF MAPPING QC RULES
#
rule mapping_qc_preparation:
  input:  ref = expand("{data}",data=ANNOTATION),
          bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR)
  output: tmp = temp(expand("{dir}/mapped_QC/temp_{data}",data=ANNOTATION, dir=OUTPUT_DIR)),
          bed12 = expand("{dir}/mapped_QC/{data}.genes.bed12", data=ANNOTATION.rstrip(".gtf"), dir=OUTPUT_DIR),
          tmp2 = temp(expand("{dir}/mapped_QC/temp_{data}.refFlat",data=ANNOTATION, dir=OUTPUT_DIR)),
          flat = expand("{dir}/mapped_QC/{data}.refFlat.txt", dir=OUTPUT_DIR, data=ANNOTATION.rstrip(".gtf")),
          head = expand("{dir}/mapped_QC/{data}.header", dir=OUTPUT_DIR, data=ANNOTATION.rstrip(".gtf")),
          list = expand("{dir}/mapped_QC/{data}.rRNA.intervalListBody.txt", dir=OUTPUT_DIR, data=ANNOTATION.rstrip(".gtf")),
          rrna = expand("{dir}/mapped_QC/{data}.rRNA.gtf", dir=OUTPUT_DIR, data=ANNOTATION.rstrip(".gtf"))
  shell:  """
          {UCSC_SCRIPTS}/gtfToGenePred -genePredExt -geneNameAsName2 {input.ref} {output.tmp}
          awk '{{print $2"\t"$4"\t"$5"\t"$1"\t0\t"$3"\t"$6"\t"$7"\t0\t"$8"\t"$9"\t"$10}}' {output.tmp} > {output.bed12}
          ####################
          #### TODO: resolve why single-exon transcripts starting from 1st position make troubles
          cat {output.bed12}  | awk '$7 != 0 || $10 != 1' > {output.tmp} && cp {output.tmp} {output.bed12}
          ####################
          {UCSC_SCRIPTS}/gtfToGenePred -genePredExt {input.ref} {output.tmp2}
          paste <(cut -f 12 {output.tmp2}) <(cut -f 1-10 {output.tmp2}) > {output.flat}
          {SAMTOOLS} view -H {input.bam} > {output.head}
          cat {output.head} > {output.list}
          grep 'gene_biotype \"rRNA' {input.ref} > {output.rrna} || echo $? > /dev/null
          cut -s -f 1,4,5,7,9 {output.rrna} >> {output.list}
          """

rule mapping_qc_picard:
  input:  bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR),
          bai = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam.bai",dir=OUTPUT_DIR),
          flat = expand("{dir}/mapped_QC/{data}.refFlat.txt", dir=OUTPUT_DIR, data=ANNOTATION.rstrip(".gtf")),
          list = expand("{dir}/mapped_QC/{data}.rRNA.intervalListBody.txt", dir=OUTPUT_DIR, data=ANNOTATION.rstrip(".gtf"))
  output: txt_fwd = expand("{dir}/picard/Aligned.sortedByCoord.output.RNA_Metrics.forward.txt", dir=OUTPUT_DIR),
          txt_rev = expand("{dir}/picard/Aligned.sortedByCoord.output.RNA_Metrics.reverse.txt", dir=OUTPUT_DIR)
  log:    run = expand("{dir}/picard/run_stats.log", dir=OUTPUT_DIR)
  params: pdf_for = expand("{dir}/picard/Aligned.sortedByCoord.npc.forward.pdf", dir=OUTPUT_DIR),
          pdf_rev = expand("{dir}/picard/Aligned.sortedByCoord.npc.reverse.pdf", dir=OUTPUT_DIR)
  threads:  4
  run:    
          shell("""
          java -jar {PICARD}/CollectRnaSeqMetrics.jar \
      		  I={input.bam} \
      		  O={output.txt_fwd} \
      		  REF_FLAT={input.flat} \
      		  STRAND=FIRST_READ_TRANSCRIPTION_STRAND \
      		  RIBOSOMAL_INTERVALS={input.list} \
      		  CHART={params.pdf_for} \
      		  VALIDATION_STRINGENCY=LENIENT > {log.run} 2>&1
      		
          java -jar {PICARD}/CollectRnaSeqMetrics.jar \
      		  I={input.bam} \
      		  O={output.txt_rev} \
      		  REF_FLAT={input.flat} \
      		  STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
      		  RIBOSOMAL_INTERVALS={input.list} \
      		  CHART={params.pdf_rev} \
      		  VALIDATION_STRINGENCY=LENIENT > {log.run} 2>&1
          """)
          
#### TODO: Ask Honza about warning - if it should stop the pipeline or not
rule mapping_qc_strandness:
  input:  txt_fwd = expand("{dir}/picard/Aligned.sortedByCoord.output.RNA_Metrics.forward.txt", dir=OUTPUT_DIR),
          txt_rev = expand("{dir}/picard/Aligned.sortedByCoord.output.RNA_Metrics.reverse.txt", dir=OUTPUT_DIR)
  output: expand("{dir}/strandness.txt", dir=OUTPUT_DIR)
  run:    
#          #with open(output[0],"w") as out:
#            with open(input.txt_fwd,"r") as file: 
#              for line in file: 
#                if not line.startswith("#") & re.match("\d",line):
#                  #FORWARD_PCT = float(line.split("\t")[17])
#                  #out.write(float(line.split("\t")[17]))
#                  storage.store("myvar", float(line.split("\t")[17]))

          R("""
          table <- read.table("{input.txt_fwd}", sep="\t", as.is=T, nrows = 2)  # TODO: nrows is temporary solution, should be done better
          fwd <- as.numeric(table[2,"V18"])
          table <- read.table("{input.txt_rev}", sep="\t", as.is=T, nrows = 2)  # TODO: nrows is temporary solution, should be done better
          rev <- as.numeric(table[2,"V18"])
          data <- "none"
          if(fwd > rev) {{ 
            if(fwd >= 0.75) data <- "forward" 
            else if(fwd < 0.55 && rev < 0.55) data <- "none" 
            else stop("Warning: Strandness is not obviously recognizable!") 
          }} else {{ 
            if(rev >= 0.75) data <- "reverse" 
            else if(fwd < 0.55 && rev < 0.55) data <- "none" 
            else stop("Warning: Strandness is not obviously recognizable!") 
          }}
          write(data, file = '{output}')
          """)
          

##### JUST FOR TESTING (DELETE FINALLY)
rule test:
  input:  expand("{dir}/strandness.txt", dir=OUTPUT_DIR)
  output:
  run:    
          with open(input[0],"r") as file:
            data = file.read().strip()
            #data = storage.fetch("myvar")
            print(data)
            if data == 'forward':
              shell(" echo {data} ")
              
            elif data == 'reverse':
              shell(" echo 'REVERSE' ")
            else:
              shell(" echo 'NONE' ")

## TODO: resolve automaticaly the insert gap size from mapper
rule mapping_qc_preseq:
  input:  bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR)
  output: extrap = expand("{dir}/preseq/Aligned.sortedByCoord.yield_estimates.txt",dir=OUTPUT_DIR),
          curve = expand("{dir}/preseq/Aligned.sortedByCoord.estimates.txt",dir=OUTPUT_DIR)
  params: paired = "" if config[CONF_DATA_TYPE_TERM] == "single" else "-pe",
          seglen = "1000000" #default is only 5000
  shell:  """
          {PRESEQ} lc_extrap -B {params.paired} -seg_len {params.seglen} -o {output.extrap} {input.bam}
          {PRESEQ} c_curve -B {params.paired} -seg_len {params.seglen} -o {output.curve} {input.bam}
          """


rule summary_biotypes:                              ################## CHECK WITH HONZA: -s 0/1/2 #####################
  input:  ref = expand("{data}",data=ANNOTATION),
          bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR)
  output: fwd = expand("{dir}/mapped_QC/featureCounts.quantSeq.fwd.biotype_counts.txt",dir=OUTPUT_DIR),
          rev = expand("{dir}/mapped_QC/featureCounts.quantSeq.rev.biotype_counts.txt",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/mapped_QC/featureCounts.run_stats.summary_biotypes.log",dir=OUTPUT_DIR)
  threads:  4
  params: fwd = expand("{dir}/mapped_QC/featureCounts.quantSeq.fwd.biotype",dir=OUTPUT_DIR),
          rev = expand("{dir}/mapped_QC/featureCounts.quantSeq.rev.biotype",dir=OUTPUT_DIR),
          paired = "" if config[CONF_DATA_TYPE_TERM] == "single" else "-p"
          
  shell:  """
          echo '###### FORWARD #######' > {log.run}
          {FEATURE_COUNTS} -t exon -g gene_biotype {params.paired} -s 1 -T {threads} -a {input.ref} -o {params.fwd} {input.bam} >> {log.run} 2>&1
          cut -f 1,7- {params.fwd} > {output.fwd}
          echo '###### REVERSE #######' >> {log.run}
          {FEATURE_COUNTS} -t exon -g gene_biotype {params.paired} -s 2 -T {threads} -a {input.ref} -o {params.rev} {input.bam} >> {log.run} 2>&1
          cut -f 1,7- {params.rev} > {output.rev}
          """

#### TODO: not working - needs to be resolved
rule RSeQC_read_distribution:
  input:  bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR),
          bed = expand("{dir}/mapped_QC/{data}.genes.bed12", data=ANNOTATION.rstrip(".gtf"), dir=OUTPUT_DIR)
  output: expand("{dir}/RSeQC/Aligned.sortedByCoord.read_distribution.txt",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/RSeQC/read_distribution.log",dir=OUTPUT_DIR)
  shell:  """
          {RSEQC}/read_distribution.py -i {input.bam} -r {input.bed} > {output} 2> {log.run}
          """
          
rule RSeQC_junction_saturation:
  input:  bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR),
          bed = expand("{dir}/mapped_QC/{data}.genes.bed12", data=ANNOTATION.rstrip(".gtf"), dir=OUTPUT_DIR)
  output: expand("{dir}/RSeQC/Aligned.sortedByCoord.junction_saturation.junctionSaturation_plot.r",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/RSeQC/junction_saturation.log",dir=OUTPUT_DIR)
  params: prefix = expand("{dir}/RSeQC/Aligned.sortedByCoord.junction_saturation",dir=OUTPUT_DIR)
  shell:  """
          {RSEQC}/junction_saturation.py -i {input.bam} -r {input.bed} -o {params.prefix} > {log.run} 2>&1
          """
       
rule RSeQC_junction_annotation:
  input:  bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR),
          bed = expand("{dir}/mapped_QC/{data}.genes.bed12", data=ANNOTATION.rstrip(".gtf"), dir=OUTPUT_DIR)
  output: expand("{dir}/RSeQC/Aligned.sortedByCoord.junction_annotation.junction.xls",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/RSeQC/junction_annotation.log",dir=OUTPUT_DIR)
  params: prefix = expand("{dir}/RSeQC/Aligned.sortedByCoord.junction_annotation",dir=OUTPUT_DIR)
  shell:  """
          {RSEQC}/junction_annotation.py -i {input.bam} -r {input.bed} -o {params.prefix} > {log.run} 2>&1
          """          
     
rule RSeQC_bam_stat:
  input:  bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR)
  output: expand("{dir}/RSeQC/Aligned.sortedByCoord.bam_stat.txt",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/RSeQC/bam_stat.log",dir=OUTPUT_DIR)
  shell:  """
          {RSEQC}/bam_stat.py -i {input.bam} > {output} 2> {log.run}
          """  
       
rule RSeQC_infer_experiment:
  input:  bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR),
          bed = expand("{dir}/mapped_QC/{data}.genes.bed12", data=ANNOTATION.rstrip(".gtf"), dir=OUTPUT_DIR)
  output: expand("{dir}/RSeQC/Aligned.sortedByCoord.infer_experiment.txt",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/RSeQC/infer_experiment.log",dir=OUTPUT_DIR)
  shell:  """
          {RSEQC}/infer_experiment.py -i {input.bam} -r {input.bed} > {output} 2> {log.run}
          """  
  
rule RSeQC_read_duplication:
  input:  bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR)
  output: expand("{dir}/RSeQC/Aligned.sortedByCoord.read_duplication.seq.DupRate.xls",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/RSeQC/read_duplication.log",dir=OUTPUT_DIR)
  params: prefix = expand("{dir}/RSeQC/Aligned.sortedByCoord.read_duplication",dir=OUTPUT_DIR)
  shell:  """
          {RSEQC}/read_duplication.py -i {input.bam} -o {params.prefix} > {log.run} 2>&1
          """
           
rule RSeQC_RPKM_saturation:
  input:  expand("{dir}/strandness.txt", dir=OUTPUT_DIR),
          bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR),
          bed = expand("{dir}/mapped_QC/{data}.genes.bed12", data=ANNOTATION.rstrip(".gtf"), dir=OUTPUT_DIR)
  output: expand("{dir}/RSeQC/Aligned.sortedByCoord.RPKM_saturation.rawCount.xls",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/RSeQC/RPKM_saturation.log",dir=OUTPUT_DIR)
  params: prefix = expand("{dir}/RSeQC/Aligned.sortedByCoord.RPKM_saturation",dir=OUTPUT_DIR)
  run:  
          with open(input[0],"r") as file:
            data = file.read().strip()
            print(data)
            if data == 'forward':
              if config[CONF_DATA_TYPE_TERM] == "single":
                shell(" {RSEQC}/RPKM_saturation.py -i {input.bam} -r {input.bed} -d '++,--' -o {params.prefix} > {log.run} 2>&1")
              else:
                shell(" {RSEQC}/RPKM_saturation.py -i {input.bam} -r {input.bed} -d '1++,1--,2+-,2-+' -o {params.prefix} > {log.run} 2>&1 ")
            elif data == 'reverse':
              if config[CONF_DATA_TYPE_TERM] == "single":
                shell(" {RSEQC}/RPKM_saturation.py -i {input.bam} -r {input.bed} -d '+-,-+' -o {params.prefix} > {log.run} 2>&1 ")
              else:
                shell(" {RSEQC}/RPKM_saturation.py -i {input.bam} -r {input.bed} -d '1+-,1-+,2++,2--' -o {params.prefix} > {log.run} 2>&1 ")
            else:
              shell(" {RSEQC}/RPKM_saturation.py -i {input.bam} -r {input.bed} -o {params.prefix} > {log.run} 2>&1 ")


###### TODO: resolve the used resources (memmory in this case)
rule picard_mark_duplicates:
  input:  bam = expand("{dir}/aligned/Aligned.sortedByCoord.out.bam",dir=OUTPUT_DIR)
  output: bam = expand("{dir}/mapped_QC/Aligned.sortedByCoord.markDups.bam",dir=OUTPUT_DIR),
          mtx = expand("{dir}/mapped_QC/Aligned.sortedByCoord.markDups_metrics.txt",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/mapped_QC/Aligned.sortedByCoord.markDups.log",dir=OUTPUT_DIR)
  threads:  4
  shell:  """
          java -Xmx2g -jar {PICARD}/MarkDuplicates.jar \
	        INPUT={input.bam} \
	        OUTPUT={output.bam} \
	        METRICS_FILE={output.mtx} \
	        REMOVE_DUPLICATES=false \
	        ASSUME_SORTED=true \
	        PROGRAM_RECORD_ID=null \
	        VALIDATION_STRINGENCY=LENIENT > {log.run} 2>&1
          """

rule dupradar_count_duplicates:
  input:  expand("{dir}/strandness.txt", dir=OUTPUT_DIR),
          bam = expand("{dir}/mapped_QC/Aligned.sortedByCoord.markDups.bam",dir=OUTPUT_DIR),
          ref = expand("{data}",data=ANNOTATION)
  output: expand("{dir}/mapped_QC/Aligned.sortedByCoord.dupRadar_dupMatrix.txt",dir=OUTPUT_DIR)
  log:    run = expand("{dir}/mapped_QC/Aligned.sortedByCoord.dupRadar.log",dir=OUTPUT_DIR)
  params: prefix = expand("{dir}/mapped_QC/Aligned.sortedByCoord.dupRadar",dir=OUTPUT_DIR),
          installation = ""
  threads:  4
  run:    
          with open(input[0],"r") as file:
            data = file.read().strip()
            print(data)
            if data == 'forward':
              if config[CONF_DATA_TYPE_TERM] == "single":
                shell(" {DUPRADAR} {input.bam} {input.ref} 'single' 1 {params.prefix} {threads} {params.installation} > {log.run} 2>&1 ")
              else:
                shell(" {DUPRADAR} {input.bam} {input.ref} 'paired' 1 {params.prefix} {threads} {params.installation} > {log.run} 2>&1 ")
            elif data == 'reverse':
              if config[CONF_DATA_TYPE_TERM] == "single":
                shell(" {DUPRADAR} {input.bam} {input.ref} 'single' 2 {params.prefix} {threads} {params.installation} > {log.run} 2>&1 ")
              else:
                shell(" {DUPRADAR} {input.bam} {input.ref} 'paired' 2 {params.prefix} {threads} {params.installation} > {log.run} 2>&1 ")
            else:
              if config[CONF_DATA_TYPE_TERM] == "single":
                shell(" {DUPRADAR} {input.bam} {input.ref} 'single' 0 {params.prefix} {threads} {params.installation} > {log.run} 2>&1 ")
              else:
                shell(" {DUPRADAR} {input.bam} {input.ref} 'paired' 0 {params.prefix} {threads} {params.installation} > {log.run} 2>&1 ")
          
          
