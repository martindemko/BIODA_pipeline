# BIODA_pipeline
This is the BIODA Group's prototype of [snakemake pipeline](http://snakemake.readthedocs.io/en/stable/index.html) for **RNA-Seq**, **WGS**, ... analysis.

## Configuration
In order to use it, you should specify your own *config.json* file. It has form of named configuration features (enclosed in double quotes) containing list of strings (values, enclosed by double quotes) separated by comma symbol. Each string belongs to different experiment and the number of all strings in a list defines also the total number of experiments and should be the same for all features (even if there would be an empty string in some cases). First we explain necessary configuration features on real example but remember that json doesn't allow commentary (lines started with symbol #) so all the next commentary lines are added just for better understanding and must be deleted:
```json
{
# define unique index value for each experiment
"index": ["S1", "S2", "S3"],

# define type of biological material (either DNA or RNA so far)
"material": ["RNA", "RNA", "RNA"],

# define type of NGS sequencing (either single or paired)
"data_type": ["single", "paired", "paired"],

# define analysis type according to the input data and expected results, 
# currently either classic or quant for RNA and
# .... for DNA
"analysis_type": ["classic", "classic", "quant"],

# define species for each experiment (for the future, by this name the default genome, 
# reference, genome index and some settings will be automatically set for pre-defined species)
"species":["trichinella", "trichinella", "trichinella"],

# define the name (or path) of the folder where data for particular experiment will be stored 
"output_dir": ["experiment_1", "experiment_2/classic-mode", "experiment_2/quant-mode"],

# define input data files' prefixes without extensions and strandness (see the examples 
# of input data files below)
"sample": ["SRR518926_shorter", "SRR609196_shorter", "SRR609197_shorter"],

# examples of input fastq data files we used for this config file (for single-end data 
# fill just RNA.fq.data1 and let RNA.fq.data2 contain an empty string)
# (supported file formats: .fq, .fastq, .fq.gz, .fastq.gz)
"RNA.fq.data1": ["/home/TS_data/SRR518926_shorter.fastq.gz", "/home/TS_data/SRR609196_shorter_R1.fastq.gz", "/home/TS_data/SRR609197_shorter_1.fastq.gz"],
"RNA.fq.data2": ["", "/home/TS_data/SRR609196_shorter_R2.fastq.gz", "/home/TS_data/SRR609197_shorter_2.fastq.gz"],

# define strandness of input fastq data files (for paired-end data use the exact 
# trailing parts used in the names of files - see example input data above - otherwise, 
# for single-end data use just an empty string "")
"strands": ["", "_R1,_R2", "_1,_2"],

# define genome files for each experiment even if they are the same 
# (supported file formats: .fa, .fasta, .fa.gz, .fasta.gz)
# (for the future, this feature will be optional)
"genome": ["/home/TS_data/Trichinella_spiralis.Tspiralis1.27.dna.genome.fa.gz", "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Trichinella_spiralis.Tspiralis1.27.dna.genome.fasta.gz", "/mnt/nfs/shared/999993-Bioda/scripts/martin/test/Trichinella_spiralis.Tspiralis1.27.dna.genome.fasta.gz"],

# define annotation files for each experiment even if they are the same
# (supported file formats: .gtf, .gff3, .gtf.gz, .gff3.gz)
# (for the future, this feature will be optional)
"annotation": ["/home/TS_data/Trichinella_spiralis.Tspiralis1.27.gtf", "/home/TS_data/Trichinella_spiralis.Tspiralis1.27.gff3", "/home/TS_data/Trichinella_spiralis.Tspiralis1.27.gff3.gz"],

# if you know the adapters used in upstream NGS which you want to trim, specify a fasta file 
# containing such sequences, otherwise they won't be trimmed and pipeline will inform you 
# (however, we recommend to run pipeline up to the preprocessing target and check for any 
# adapters identified by preprocessing in the results at "output_dir"/minion/.. and to define
# own fasta file for each experiment by hand if necessary and then to continue in pipeline)
"adapters": ["", "adapters_1.fa", "adapters_2.fa"],
```
