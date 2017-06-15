
from snakemake.shell import shell

TRIMMOMATIC = "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Tools/trimmomatic-master/classes/trimmomatic.jar"

shell(""" java -jar {TRIMMOMATIC} SE -threads {snakemake.threads} {snakemake.params.phred} {snakemake.input.reads} {snakemake.output.clean} LEADING:{snakemake.params.leading} TRAILING:{snakemake.params.trailing} CROP:{snakemake.params.crop} SLIDINGWINDOW:{snakemake.params.slid_w_1}:{snakemake.params.slid_w_2} MINLEN:{snakemake.params.minlen} > {snakemake.log.run} 2>&1 """)
