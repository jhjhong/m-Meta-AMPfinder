from script.Database import Database
from script.Fastp import Fastp
from script.Megahit import Megahit
from script.Spades import Spades
from script.ORF import PyORF
from script.Blast import Blast
from script.Diamond import Diamond
from script.Predict import Predict
from script.settings import *
from script.Base import*

from Bio import SeqIO
import filetype
import os
import time
import gzip, zlib
import bz2


class AMP(AMPBase):
    def __init__(self, input_type='', short1=None, short2=None, input_sequence=None , threads=32, output_dir=None, data='na', assembler='megahit', aligner='blast', debug=False):
        # Change some arguments to full paths.
        if short1:
            self.short1 = os.path.abspath(short1)
        if short2:
            self.short2 = os.path.abspath(short2)
        else:
            self.short2 = None
        if input_sequence:
            self.input_sequence = os.path.abspath(input_sequence)

        self.input_type = input_type.lower()
        self.threads = threads
        self.output_dir = os.path.abspath(output_dir)
        self.assembler = assembler.lower()
        self.aligner = aligner.lower()

        self.debug = debug
        if self.debug:
	        logger.setLevel(10)
        
        # database folder
        self.db = path
        self.dp = data_path

        super(AMPBase, self).__init__()
    
    def make_output_directory(self):
        # Creates the output directory, if it doesn't already exist.
        if not os.path.exists(self.output_dir):
            try:
                os.makedirs(self.output_dir)
            except OSError:
                logger.error("mACPfinder was unable to make the output directory")
                exit()
            logger.info("Making output directory: {}".format(self.output_dir))
        elif os.listdir(self.output_dir):
            logger.warning("The output directory {} already exists, please change the parameter -o to another value to avoid overwriting.".format(self.output_dir))
            # exit()
        else:  # directory exists but is empty
            logger.info("The output directory already exists: {}".format(self.output_dir))

    def validate_inputs(self):
        # Checks that args are consistent
        #    - validate input file name and out file name
		# 	 - validation for mutually exclusive options e.g. protein sequence for contig input_type etc
        # Exits with an error message if not the case
        if self.input_type == "read":
            if not self.short1:
                if not self.input_sequence:
                    logger.error("FQ file is necessary for 'reads' command.")
                    exit()
                else:
                    self.short1 = self.input_sequence
            if not os.path.exists(self.short1):
                logger.error("input file does not exist: {}".format(self.short1))
                exit()
            else:
                # assign umcompressed short1 
                logger.info("Input short1: {} => {}".format(self.short1, filetype.guess(self.short1)))
                if self.check_filetype(self.short1) == None:
                    pass
                else:
                    self.short1 = self.check_filetype(self.short1)

                if self.short2:
                    if not os.path.exists(self.short2):
                        logger.error("input file does not exist: {}".format(self.short2))
                        exit()
                    elif self.short1 and self.short2 and self.short1 == self.short2:
                        logger.error("first and second read pair files cannot be the same file")
                        exit()
                    # assign umcompressed short2
                    logger.info("Input short2: {} => {}".format(self.short2, filetype.guess(self.short2)))
                    if self.check_filetype(self.short2) == None:
                        pass
                    else:
                        self.short2 = self.check_filetype(self.short2)
                else:
                    logger.info("Input short2: {}".format(self.short2))

        elif self.input_type in ["contig", "peptide"]:
            if not self.input_sequence:
                logger.error("FASTA File is necessary for 'contigs/peptides' command.")
            if not os.path.exists(self.input_sequence):
                logger.error("input file does not exist: {}".format(self.input_sequence))
                exit()
            else:
                logger.info("Input sequence: {} => {}".format(self.input_sequence, filetype.guess(self.input_sequence)))
                if self.check_filetype(self.input_sequence) == None:
                    pass
                else:
                    self.input_sequence = self.check_filetype(self.input_sequence)


        if self.threads > os.cpu_count():
            logger.error(
                "Argument num_threads illegal value, expected (>=1 and =<{}):  given `{}`)".format(os.cpu_count(),
                                                                                                   self.threads))
            exit()

    def check_filetype(self, file):
        kind = filetype.guess(file)
        if kind is None:
            # print("TEST THIS")
            # print("CHECK", self.is_fasta(file))
            if self.is_fasta(file) == False:
                logger.error("invalid fasta")
                exit()
            else:
                return
        else:
            if kind.extension in ["gz","bz2"]:
                # uncompressed input and use uncompressed file
                filename = os.path.basename(file).split(".")[0]
                umcompressed_file = os.path.join(self.output_dir, "{}.temp.uncompressed.fsa".format(filename))
                with open(umcompressed_file, "w") as file_out:
                    if kind.extension == "gz":
                        with gzip.open(file, "rt") as handle:
                            file_out.write(handle.read())
                    else:
                        with bz2.open(file, "rt") as handle:
                            file_out.write(handle.read())

                if self.is_fasta(umcompressed_file) == False:
                    logger.error("invalid fasta")
                    exit()
                return umcompressed_file
            else:
                logger.error("Sorry, no support for file format {}".format(kind.mime))
                exit()

    def is_fasta(self, file, extension=None):
        if extension is None:
            with open(file, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
# 有bug
                # print("RIGHT")
                self.check_record(fasta)
                return True
        elif extension in ["gz", "bz2"]:
            if extension == "gz":
                # print("GZ")
                # print(file)
                with gzip.open(file, "rt") as handle:
                    fasta = SeqIO.parse(handle, "fasta")
                    # print("check_record")
                    print(self.check_record(fasta))
            else:
                with bz2.open(file, "rt") as handle:
                    fasta = SeqIO.parse(handle, "fasta")
                    self.check_record(fasta)
            return True
        else:
            return False

    def check_record(self, fasta):
        # check each record in the file
 #       print("!!!")
 #       print(self.input_type)
 #       print(fasta)
        for record in fasta:
# for loop 進不來, hen 怪
            print("HERE")
            if any(record.id) == False or any(record.seq) == False:
                return False
            if self.input_type in ["read", "contig"]:
                print("READS")
                return self.is_dna(record.seq)
            if self.input_type == "peptide":
                return self.is_protein(record.seq)
    
    # @staticmethod
    def is_dna(sequence):
        #  dna codes
        nucleotide_dict = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0, 'U': 0,
                           #  other dna codes
                           'W': 0,  # W = A or T
                           'S': 0,  # S = C or G
                           'M': 0,  # M = A or C
                           'K': 0,  # K = G or T
                           'R': 0,  # R = A or G
                           'Y': 0,  # Y = C or T
                           'B': 0,  # B = C, G, or T
                           'D': 0,  # D = A, G, or T
                           'H': 0,  # H = A, C, or T
                           'V': 0  # V = A, C, or G
                           }

        for base in sequence:
            try:
                nucleotide_dict[base.upper()] += 1
            except Exception as e:
                logger.error("invalid nucleotide fasta due to: {}".format(e))
                return False
        logger.info("valid nucleotide fasta: {}".format(nucleotide_dict))
        print("dna")
        return True

    @staticmethod
    def is_protein(sequence):
        amino_acids_dict = {
            # common symbols between protein and dna codes
            'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N': 0, 'U': 0,
            # other amino acids
            'R': 0, 'D': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': 0,
            'L': 0, 'K': 0, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0,
            'W': 0, 'Y': 0, 'V': 0, 'X': 0, 'Z': 0, 'J': 0, 'B': 0
        }
        count = 0
        for amino_acid in sequence:
            try:
                amino_acids_dict[amino_acid.upper()] += 1
            except Exception as e:
                logger.error("invalid protein fasta due to: {}".format(e))
                return False

        for a in amino_acids_dict.keys():
            if a not in 'ATGCNU':
                count = count + amino_acids_dict[a]

        if count == 0:
            logger.error("invalid protein fasta: {}".format(amino_acids_dict))
            return False

        logger.info("valid protein fasta: {}".format(amino_acids_dict))
        return True


    def run(self):
        self.make_output_directory()
        self.validate_inputs()
        if self.input_type == "read":
            self.qc_inputs()
        if self.input_type == "read":
            self.run_assembly()
        if self.input_type in ["read", "contig"]:
            self.run_smorfs()
        if self.input_type in ["read", "contig", "peptide"]:
            self.run_blast()
        if self.input_type in ["read", "contig", "peptide"]:
            self.run_pred()
    
    def qc_inputs(self):
        # run fastp for quality control: remove N base.
        try:
            if self.short2:
                qc_obj = Fastp(short1=self.short1, short2=self.short2, output_dir=os.path.join(self.output_dir, "temp.qc"), num_threads=self.threads)
                qc_obj.run()
            else:
                qc_obj = Fastp(short1=self.short1, output_dir=os.path.join(self.output_dir, "temp.qc"), num_threads=self.threads)
                qc_obj.run()

        except Exception as e:
            logger.exception("failed to write orf file")
        else:
            pass
    
    def run_assembly(self):
        # Runs assembly.
        if self.assembler == "megahit":
            self.assembly_megahit()
        elif self.assembler == "metaspades":
            self.assembly_metaspades()
        else:
            exit()

    def assembly_megahit(self):
        # Process assembly paired reads.
        short1 = os.path.basename(self.short1)
        highqual_fsa_file_1 = os.path.join(self.output_dir, "temp.qc/{}.temp.highqual.fq".format(short1.split(".")[0]))
        if self.short2:
            short2 = os.path.basename(self.short2)
            highqual_fsa_file_2 = os.path.join(self.output_dir, "temp.qc/{}.temp.highqual.fq".format(short2.split(".")[0]))
        try:
            logger.info("Run Megahit assembler")
            if self.short2:
                if os.stat(highqual_fsa_file_1).st_size > 0 and os.stat(highqual_fsa_file_2).st_size > 0:
                    megahit_obj = Megahit(input1=highqual_fsa_file_1, input2=highqual_fsa_file_2, output_dir=self.output_dir, num_threads=self.threads)
                    megahit_obj.run()
                else:
                    logger.error("Each sequence quality are low!")
            else:
                if os.stat(highqual_fsa_file_1).st_size > 0 :
                    megahit_obj = Megahit(input1=highqual_fsa_file_1, output_dir=self.output_dir, num_threads=self.threads)
                    megahit_obj.run()
                else:
                    logger.error("Each sequence quality are low!")
        except Exception as e:
            logger.exception("failed to write orf file")
        else:
            pass
    
    def assembly_metaspades(self):
        # Process assembly paired reads.
        short1 = os.path.basename(self.short1)
        highqual_fsa_file_1 = os.path.join(self.output_dir, "temp.qc/{}.temp.highqual.fq".format(short1.split(".")[0]))
        if self.short2:
            short2 = os.path.basename(self.short2)
            highqual_fsa_file_2 = os.path.join(self.output_dir, "temp.qc/{}.temp.highqual.fq".format(short2.split(".")[0]))
        try:
            logger.info("Run Metaspades assembler")
            if self.short2:
                if os.stat(highqual_fsa_file_1).st_size > 0 and os.stat(highqual_fsa_file_2).st_size > 0:
                    metaspades_obj = Spades(input1=highqual_fsa_file_1, input2=highqual_fsa_file_2, output_dir=self.output_dir, num_threads=self.threads)
                    metaspades_obj.run()
                else:
                    logger.error("Each sequence quality are low!")
            else:
                if os.stat(highqual_fsa_file_1).st_size > 0 :
                    metaspades_obj = Spades(input1=highqual_fsa_file_1, output_dir=self.output_dir, num_threads=self.threads)
                    metaspades_obj.run()
                else:
                    logger.error("Each sequence quality are low!")
        except Exception as e:
            logger.exception("failed to write orf file")
        else:
            pass

    def run_smorfs(self):
        # run pyrodigal to predict genes.
        logger.info("Run Pyrodigal ORFfinder")
        input_file = os.path.join(self.output_dir, "final.contigs.fasta")

        self.clean = True
        self.low_quality = False

        try:
            if os.stat(input_file).st_size > 0:
                orf_obj = PyORF(input_file=input_file, output_dir=os.path.join(self.output_dir, "temp.smorfs"), num_threads=self.threads, clean=self.clean, low_quality=self.low_quality)
                orf_obj.run()

            else:
                logger.error("The contigs file are empty!")

        except Exception as e:
            logger.exception("failed to write orf file")
        else:
            pass


    def run_blast(self):
        # Process blast sequence(s).
        logger.info("Run alignment against AMPdb for known AMPs")
        if self.input_type == "peptide":
            input_file = self.input_sequence
        elif self.input_type in ["read", "contig"]:
            input_file = os.path.join(self.output_dir, "final.smorfs.fsa")

        try:
            if os.stat(input_file).st_size > 0:
                if self.aligner == "diamond":
                    diamond_obj = Diamond(input_file=input_file, output_dir=os.path.join(self.output_dir, "temp.alignment"), num_threads=self.threads)
                    diamond_obj.run()
                else:
                    blast_obj = Blast(input_file=input_file, output_dir=os.path.join(self.output_dir, "temp.alignment"), num_threads=self.threads)
                    blast_obj.run()
            else:
                self.write_stub_output_file()
        except Exception as e:
            logger.exception("failed to write orf file")
        else:
            pass

    def write_stub_output_file(self):
        # write empty output file if there are no open reading frames
        with open(os.path.join(self.output_dir, "out.txt"), 'w') as fout:
            fout.write(json.dumps({}))


    def run_pred(self):
        # run functional predict of AMPs.
        logger.info("Run functional predict of AMPs")
        input_file = os.path.join(self.output_dir, "final.smorfs.L100.fsa")

        try:
            if os.stat(input_file).st_size > 0:
                ampfinder_obj = Predict(input_file=input_file, output_dir=os.path.join(self.output_dir, "temp.prediction"), num_threads=self.threads)
                ampfinder_obj.run()

            else:
                logger.error("The contigs file are empty!")

        except Exception as e:
            logger.exception("failed to write orf file")
        else:
            pass