# BIODA_pipeline
This is the BIODA Group's prototype of [snakemake pipeline](http://snakemake.readthedocs.io/en/stable/index.html) for **RNA-Seq**, **VCF**, ... analysis.
...In order to use it, you should specify your own *config.yaml* file of the form:
```yaml
# use one of [single|paired]
data_type: paired
# define analysis type according to the input data and expected results, 
# currently one of [rna-quant|rna-classic|vcf]
analysis_type: rna-classic
# define the name (might be with path) of the folder where all data (for this analysis) will be stored
output_dir: name_of_output_dir
# define path to the folder where input data should be found
input_dir: /mnt/nfs/shared/999993-Bioda/scripts/martin/test/TS_data
# list of data files' prefixes without extensions and strandness (see the examples of input data files below)
sample: SRR609195_shorter,
        SRR609196_shorter,
        SRR609197_shorter
# if paired data are used it must contain trailing parts used in the names of files, 
# otherwise it must contain just empty string ""
strands:  "_1,_2"
# examples of data files we used for this config file
#data1:  SRR609195_shorter_1.fastq.gz,
#        SRR609196_shorter_1.fastq.gz,
#        SRR609197_shorter_1.fastq.gz
#data2:  SRR609195_shorter_2.fastq.gz,
#        SRR609196_shorter_2.fastq.gz,
#        SRR609197_shorter_2.fastq.gz
# genome might be defined either by one of shortcuts [human|mouse|trichinella|sifilis] in which case 
# the pre-defined genome, annotation and some settings will be used for analysis or by a genome file 
# with whole path to it you want to use instead of the pre-defined one (in that case you have to define 
# also an annotation - see below)
genome: /mnt/nfs/shared/999993-Bioda/scripts/martin/test/Trichinella_spiralis.Tspiralis1.27.dna.genome.fa.gz
# annotation is defined either along with the pre-defined genome as the pre-defined annotation or 
# by a list of annotation files with whole path to them which you want to use instead of the pre-defined one 
# (it must be defined by a file path if the genome is also defined by a file path)
annotation: /mnt/nfs/shared/999993-Bioda/scripts/martin/test/Trichinella_spiralis.Tspiralis1.27.gff3,
            /mnt/nfs/shared/999993-Bioda/scripts/martin/test/Trichinella_spiralis.Tspiralis1.27.new.gtf.gz
# if you know the adapters used in downstream analysis which you want to trim, specify a fasta file containing 
# such sequences, otherwise they won't be trimmed (we recommend to run pipeline up to the report_preprocessing 
# rule and check for any adapters identified by preprocessing in the results ...)
adapters: /mnt/nfs/shared/999993-Bioda/scripts/martin/test/adapters_merge.fa
```
