# m(Meta)AMPfinder v.1.1.0
To provide more opportunities for clinical solutions, we were thus motivated to design a homology-based gene prediction program (AMPfinder), an integrated tool stream that combines ORF prediction, and AMP classification to extract AMPs directly from genome or proteome sequence data.

1. AMPfinder is a simple, yet accurate, computational pipeline that processes either genomes and proteome sequences. The search for AMPs is based on alignment searching the existing antimicrobial peptide database and predicting on the feature model in amino acid sequence obtained from the translation of the original transcriptome sequence data. 

2. AMPfinder is a free software tool that combines ORF prediction and accurately classification of the AMPs from protein or nucleotide data, which performs translation of the input transcriptome data by using Prodigal, and selects short sequences containing ORF and signal peptide cleavage sites. Then, by default, DIAMOND and BLAST were used for homology detection and machine learning prediction model were used to search for potential AMPs, in which case all known or potential motifs will be revealed and classified. Due to the combination of various search methods, AMPfinder searcher allows to obtain the most complete repertoire of AMPs for one or more transcriptomes in a short amount of time. Therefore, AMPfinder seems to be the most suitable tool for rapid screening of potential AMPs.

<p align="center"><img src="images/workflow.png" alt="AMPfinder" width="600"></p>

##  Publications

> Jhong, J.H., et al., 
> dbAMP 2.0: updated resource for antimicrobial peptides with an enhanced scanning method for genomic and proteomic data. Nucleic Acids Res (2021).
> [10.1093/nar/gkab1080](https://doi.org/10.1093/nar/gkab1080)


## Install Development Version
### Install Dependencies
The following dependencies are required.
- python 3.8
- Biopython 1.79
- Fastp 0.23.0
- megahit 1.2.9
- Spades 3.15.5
- BLAST
- DIAMOND
- Prodigal 2.6.3
- modlamp 4.2.1
- scikit-learn 0.24.1
- pytest 3.0.0+

### Install AMPfinder

We have already integrate the environment in `conda_env.yml`. 
execute `conda env create --file conda_env.yml` to install packages required in a new created `AMPfinder` conda env.

Enter the enviornment with `conda activate AMPfinder` before further executions.

**Note**: Please ensure that all dependencies are installed before using AMPfinder

## Getting started
### Construct Database
We provide a JSON format file, **ampfinder_210729.json**, which contains AMP sequences from [dbAMP](https://awi.cuhk.edu.cn/~dbAMP/ "dbAMP"). You can directly use this file and run the following command to construct the database.

`macpfinder load --input_json ampfinder_210729.json`

Then you will see the following files in your directory of `db/`.
- _data/proteindb.fsa
- _db/protein_all.db.pin
- _db/protein_all.db.phr
- _db/protein_all.db.psq
- _db/protein_all.db.dmnd

### Running AMPfinder

We provided two test files: **example_protein.fasta** and **example_dna.fasta**. 
The following command will bring up AMPfinder's main help menu:

`python mampfinder main -h`

    AMPfinder - 1.1.0 - main
    
    optional arguments:
      -h, --help            show this help message and exit
      -1 SHORT1, --short1 SHORT1
                        FASTQ file (metagenomics reads required) of first short reads in each pair
      -2 SHORT2, --short2 SHORT2
                        FASTQ file (metagenomics reads required) of second short reads in each pair
      -i INPUT_SEQUENCE, --input_sequence INPUT_SEQUENCE
                            input file must be in FASTA (contig and peptide required) format!
      -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        path to the output directory (required)
      -t {read,contig,peptide}, --input_type {read,contig,peptide}
                            specify data input type (required)
      -n THREADS, --num_threads THREADS
                            number of threads (CPUs) to use in the BLAST search (default=16)
      --assembler {megahit,metaspades}
                            specify assembler tool (default = megahit)
      --alignment_tool {DIAMOND,BLAST}
                            specify alignment tool (default = BLAST)
      -v, --version         show mACPfinder software version number
      --debug               debug mode
