# AMPfinder V.0.0.3
To provide more opportunities for clinical solutions, we were thus motivated to design a homology-based gene prediction program (AMPfinder), an integrated tool stream that combines ORF prediction, and AMP classification to extract AMPs directly from genome or proteome sequence data.

1. AMPfinder is a simple, yet accurate, computational pipeline that processes either genomes and proteome sequences. The search for AMPs is based on alignment searching the existing antimicrobial peptide database and predicting on the feature model in amino acid sequence obtained from the translation of the original transcriptome sequence data. 

2. AMPfinder is a free software tool that combines ORF prediction and accurately classification of the AMPs from protein or nucleotide data, which performs translation of the input transcriptome data by using Prodigal, and selects short sequences containing ORF and signal peptide cleavage sites. Then, by default, DIAMOND and BLAST were used for homology detection and machine learning prediction model were used to search for potential AMPs, in which case all known or potential motifs will be revealed and classified. Due to the combination of various search methods, AMPfinder searcher allows to obtain the most complete repertoire of AMPs for one or more transcriptomes in a short amount of time. Therefore, AMPfinder seems to be the most suitable tool for rapid screening of potential AMPs.


## Install Development Version
### Install Dependencies
The following dependencies are required.
- python 3.6
- BLAST
- DIAMOND
- Prodigal 2.6.3
- Biopython1.73
- filetype 1.0.0+
- pytest 3.0.0+
- pandas 0.15.0+
- seaborn 0.8.1+
- pyfaidx 0.5.4.1+
- pyahocorasick 1.1.7+

### Install AMPfinder

`git clone https://github.com/BiOmicsLab/AMPfinder.git`

`cd AMPfinder`

`python AMPfinder main -h`

**Note**: Please ensure that all dependencies are installed before using AMPfinder

### Usage
#### Construct Database
We provide a JSON format file, **ampfinder_210729.json**, which contains AMP sequences from [dbAMP](https://awi.cuhk.edu.cn/~dbAMP/ "dbAMP"). You can directly use this file and run the following command to construct the database.

`python makedatabase.py database -i ampfinder_210729.json`

Then you will see the following files in your directory.
- protein_all.db.pin
- protein_all.db.phr
- protein_all.db.psq
- protein_all.db.dmnd
- database.fasta

#### Running AMPfinder

We provided two test files: **example_protein.fasta** and **example_dna.fasta**. 
The following command will bring up AMPfinder's main help menu:
`python ampfinder main -h`

    AMPfinder - 0.0.3 - main
    
    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT_SEQUENCE, --input_sequence INPUT_SEQUENCE
                            input file must be in FASTA (contig and protein) format! e.g myFile.fasta
      -o OUTPUT_FILE, --output_file OUTPUT_FILE
                            output folder and base filename
      -t {contig,protein}, --input_type {contig,protein}
                            specify data input type (default = contig)
      -a {DIAMOND,BLAST}, --alignment_tool {DIAMOND,BLAST}
                            specify alignment tool (default = BLAST)
      -n THREADS, --num_threads THREADS
                            number of threads (CPUs) to use in the BLAST search (default=48)
