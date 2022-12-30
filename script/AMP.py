from script.Database import Database
from script.ORF import ORF
from script.Blast import Blast
from script.Diamond import Diamond
from script.settings import *
from script.Base import*

from Bio import SeqIO
import filetype
import os
import time
import gzip, zlib
import bz2


class AMP(AMPBase):
    def __init__(self, input_type='', short1=None, short2=None, input_sequence=None , threads=32, output_dir=None, data='na', aligner='blast', debug=False):
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
            logger.info("The output directory already exists and files may be reused or overwritten: {}".format(self.output_dir))
        else:  # directory exists but is empty
            logger.info("The output directory already exists: {}".format(self.output_dir))

    def validate_inputs(self):
        # Checks that args are consistent
        #    - validate input file name and out file name
		# 	 - validation for mutually exclusive options e.g. protein sequence for contig input_type etc
        # Exits with an error message if not the case
        if self.input_type == "reads":
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
                logger.info("Input short1: {} => {}".format(self.short1, filetype.guess(self.short1)))
                if self.check_filetype(self.short1) == False:
                    logger.error("invalid input file format")
                    exit()
                if self.short2:
                    if not os.path.exists(self.short2):
                        logger.error("input file does not exist: {}".format(self.short2))
                        exit()
                    elif self.short1 and self.short2 and self.short1 == self.short2:
                        logger.error("first and second read pair files cannot be the same file")
                        exit()

                    logger.info("Input short2: {} => {}".format(self.short2, filetype.guess(self.short2)))
                    if self.check_filetype(self.short2) == False:
                        logger.error("invalid input file format")
                        exit()
                else:
                    logger.info("Input short2: {}".format(self.short2))
                print("HE")
        
        elif self.input_type in ["contigs", "peptides"]:
            if not self.input_sequence:
                logger.error("FASTA File is necessary for 'contigs/peptides' command.")
            if not os.path.exists(self.input_sequence):
                logger.error("input file does not exist: {}".format(self.input_sequence))
                exit()
            else:
                logger.info("{} => {}".format(self.input_sequence, filetype.guess(self.input_sequence)))


        if self.threads > os.cpu_count():
            logger.error(
                "Argument num_threads illegal value, expected (>=1 and =<{}):  given `{}`)".format(os.cpu_count(),
                                                                                                   self.threads))
            exit()

    def check_filetype(self, file):
        kind = filetype.guess(file)
        if kind is None:
            if self.is_fasta(file) == False:
                logger.error("invalid fasta")
                exit()
            else:
                return True
        # else:
        #     if kind.extension in ["gz","bz2"]:
        #         if self.is_fasta(kind.extension) == False:
        #             logger.error("invalid fasta")
        #             exit()
        #         # uncompressed input and use uncompressed file
        #         filename = os.path.basename(self.input_sequence)
        #         umcompressed_file = os.path.join(self.output_dir, "{}.temp.uncompressed.fsa".format(filename))
        #         with open(umcompressed_file, "w") as file_out:
        #             if kind.extension == "gz":
        #                 with gzip.open(self.input_sequence, "rt") as handle:
        #                     file_out.write(handle.read())
        #             else:
        #                 with bz2.open(self.input_sequence, "rt") as handle:
        #                     file_out.write(handle.read())

        #         self.input_sequence = umcompressed_file
        #         self.umcompressed_file = umcompressed_file
        #     else:
        #         logger.error("Sorry, no support for file format {}".format(kind.mime))
        #         exit()

    def is_fasta(self, file, extension=None):
        if extension is None:
            with open(file, "r") as handle:
                fasta = SeqIO.parse(handle, "fasta")
                self.check_record(fasta)
                return True
        elif extension in ["gz", "bz2"]:
            if extension == "gz":
                with gzip.open(file, "rt") as handle:
                    fasta = SeqIO.parse(handle, "fasta")
                    self.check_record(fasta)
            else:
                with bz2.open(file, "rt") as handle:
                    fasta = SeqIO.parse(handle, "fasta")
                    self.check_record(fasta)
            return True
        else:
            return False

    def check_record(self, fasta):
        # check each record in the file
        for record in fasta:
            if any(record.id) == False or any(record.seq) == False:
                return False
            if self.input_type == "contig":
                return self.is_dna(record.seq)
            if self.input_type == "protein":
                return self.is_protein(record.seq)
    
    # @staticmethod
    def is_dna(self, sequence):
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
    def is_protein(self, sequence):
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
        # print("TEST1")
        # print(self.output_file)
        self.make_output_directory()
        self.validate_inputs()
        # self.run_blast()

    def run_blast(self):
        """Runs blast."""
        if self.input_type == "protein":
            self.process_protein()
        elif self.input_type == "contig":
            self.process_contig()
        else:
            exit()

    def process_metagenomics(self):
        """Process metagenomics short reads."""
        short1 = os.path.basename(self.short1)
        short2 = os.path.basename(self.short2)
        print("HERE")
        output = self.output_file
        orf_obj = ORF(input_file=self.input_sequence)
        orf_obj.contig_to_orf()
        contig_fsa_file = os.path.join(self.working_directory, "{}.temp.contig.fsa".format(file_name))
        try:
            if os.stat(contig_fsa_file).st_size > 0:
                if self.aligner == "diamond":
                    diamond_obj = Diamond(input_file=contig_fsa_file, output_file = output,num_threads=self.threads)
                    diamond_obj.run()
                else:
                    blast_obj = Blast(input_file=contig_fsa_file, output_file = output, num_threads=self.threads)
                    blast_obj.run()
            else:
                self.write_stub_output_file()
        except Exception as e:
            logger.exception("failed to write orf file")
        else:
            pass


    def process_protein(self):
        """Process protein sequence(s)."""
        file_name = os.path.basename(self.input_sequence)
        output = self.output_file
        if self.aligner == "diamond":
            diamond_obj = Diamond(self.input_sequence, output_file = output, num_threads=self.threads)
            diamond_obj.run()
        else:
            blast_obj = Blast(input_file= file_name,  output_file = output, num_threads=self.threads)
            blast_obj.run()


    def process_contig(self):
        """Process nuclotide sequence(s)."""
        file_name = os.path.basename(self.input_sequence)
        output = self.output_file
        orf_obj = ORF(input_file=self.input_sequence)
        orf_obj.contig_to_orf()
        contig_fsa_file = os.path.join(self.working_directory, "{}.temp.contig.fsa".format(file_name))
        try:
            if os.stat(contig_fsa_file).st_size > 0:
                if self.aligner == "diamond":
                    diamond_obj = Diamond(input_file=contig_fsa_file, output_file = output,num_threads=self.threads)
                    diamond_obj.run()
                else:
                    blast_obj = Blast(input_file=contig_fsa_file, output_file = output, num_threads=self.threads)
                    blast_obj.run()
            else:
                self.write_stub_output_file()
        except Exception as e:
            logger.exception("failed to write orf file")
        else:
            pass




