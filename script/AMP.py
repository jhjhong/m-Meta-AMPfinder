import filetype
import os
import time
from script.Database import Database
from script.ORF import ORF
from script.Blast import Blast
from script.Diamond import Diamond
from script.settings import *
from Bio import SeqIO
from script.Base import*


class AMP(AMPBase):
    def __init__(self, input_type='contig', input_sequence=None , threads=32, output_file=None, data='na', aligner='blast'):
        o_f_path, o_f_name = os.path.split(os.path.abspath(output_file))

        self.input_type = input_type.lower()
        self.input_sequence = os.path.abspath(input_sequence)
        self.threads = threads
        self.output_file = os.path.abspath(output_file)
        self.data = data
        self.aligner = aligner.lower()

        self.db = path
        self.dp = data_path

        self.working_directory = o_f_path

        super(AMPBase, self).__init__()

    def validate_inputs(self):

        if not os.path.exists(self.input_sequence):
            logger.error("input file does not exist: {}".format(self.input_sequence))
            exit()

        if self.output_file == self.input_sequence and self.clean:
            logger.error("output path same as input, must specify "
                         "different path when cleaning to prevent "
                         "accidental deletion of input files")
            exit()

        logger.info("{} => {}".format(self.input_sequence, filetype.guess(self.input_sequence)))
        kind = filetype.guess(self.input_sequence)

        if kind is None:
            if self.is_fasta() == False:
                logger.error("invalid fasta")
                exit()
        else:
            logger.error(kind.extension)
            logger.error(kind.mime)
            logger.warning("Sorry, no support for this format.")
            exit()
        if self.threads > os.cpu_count():
            logger.error(
                "Argument num_threads illegal value, expected (>=1 and =<{}):  given `{}`)".format(os.cpu_count(),
                                                                                                   self.threads))
            exit()

    def is_fasta(self):
        """Checks for valid fasta format."""
        with open(self.input_sequence, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            # check each record in the file
            for record in fasta:
                if any(record.id) == False or any(record.seq) == False:
                    return False
                if self.input_type == "contig":
                    return self.is_dna(record.seq)
                if self.input_type == "protein":
                    return self.is_protein(record.seq)
            return True

    @staticmethod
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

        self.validate_inputs()
        self.run_blast()

    def run_blast(self):
        """Runs blast."""
        if self.input_type == "protein":
            self.process_protein()
        elif self.input_type == "contig":
            self.process_contig()
        else:
            exit()


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




