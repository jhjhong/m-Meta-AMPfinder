import math
import filetype
import joblib
import numpy as np
import os
import time
import shutil
from glob import glob
from Bio.Blast import NCBIXML
from Database import Database
from ORF import ORF
from Blast import Blast
from Diamond import Diamond
from settings import *
from Bio import SeqIO
from Base import*

def mycopyfile(srcfile, dstpath):  # 复制函数

    if not os.path.isfile(srcfile):

        print("%s not exist!" % (srcfile))

    else:

        fpath, fname = os.path.split(srcfile)  # 分离文件名和路径

        if not os.path.exists(dstpath):
            os.makedirs(dstpath)  # 创建路径

        shutil.copy(srcfile, dstpath + fname)  # 复制文件

        print("copy %s -> %s" % (srcfile, dstpath + fname))


class AMP(AMPBase):
    """Class to predict resistome(s) from protein or nucleotide data based on CARD detection models."""

    def __init__(self, input_type='contig', input_sequence=None , threads=32, output_file=None, \
                 clean=True, data='na', aligner='blast',tax_class = None):
        """Creates RGI object for resistome(s) prediction."""

        o_f_path, o_f_name = os.path.split(os.path.abspath(output_file))

        self.input_type = input_type.lower()
        self.input_sequence = os.path.abspath(input_sequence)
        self.threads = threads
        self.output_file = os.path.abspath(output_file)
        self.clean = clean
        self.data = data
        self.aligner = aligner.lower()
        self.tax_class = tax_class


        self.db = path
        self.dp = data_path

        self.working_directory = o_f_path

        super(AMPBase, self).__init__()

    def __repr__(self):
        """Returns RGI class full object."""
        return "RGI({}".format(self.__dict__)

    @classmethod
    def from_string(cls, cmd_string):
        """Creates RGI object from string."""
        input_type, input_sequence, threads, num_sequences, output_file, aligner, database = cmd_string.split('@')
        return cls(input_type, input_sequence, threads, num_sequences, output_file, aligner, database)

    @classmethod
    def from_args(cls, *initial_data, **kwargs):
        """Creates RGI object from args."""
        for dictionary in initial_data:
            for key in dictionary:
                if key in ['input_type', 'loose', 'clean', 'aligner']:
                    setattr(cls, key, dictionary[key].lower())
                setattr(cls, key, dictionary[key])

        for key in kwargs:
            if key in ['input_type', 'loose', 'clean', 'aligner']:
                setattr(cls, key, kwargs[key].lower())
            setattr(cls, key, kwargs[key])

        return cls()

    def validate_inputs(self):
        """Validate inputs.

            - validate input file name and out file name
            - validation for mutually exclusive options e.g. protein sequence for contig input_type etc
        """
        if not os.path.exists(self.input_sequence):
            logger.error("input file does not exist: {}".format(self.input_sequence))
            exit()

        # otherwise you blow up your input when deleting intermediate files
        if self.output_file == self.input_sequence and self.clean:
            logger.error("output path same as input, must specify "
                         "different path when cleaning to prevent "
                         "accidental deletion of input files")
            exit()


        kind = filetype.guess(self.input_sequence)

        if kind is None:
            if self.is_fasta() == False:
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
                return False
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


        return True


    def run(self):
        """Runs AMP."""

        t0 = time.time()
        self.validate_inputs()
        self.run_blast()
        self.makejsonfile()
        self.getpredictionresult()
        # logger.info("Output......")
        # self.out()

    def getpredictionresult(self):
        if self.input_type == "protein" :

            print("11111111111111111111111111111112131324")
            para1 = time.strftime("%Y%m%d", time.localtime())
            para2 = time.strftime("%H%M%S", time.localtime())
            src_dir = self.input_sequence
            dst_dir = '/home/AMPfinder/public_html/test/temp/' + para1 + '/'
            print(src_dir)
            print(dst_dir)
            if not os.path.isdir(dst_dir):
                os.makedirs(dst_dir)

            if not os.path.isdir("/home/AMPfinder/public_html/test/weka-3-8-1/weka_temp/"+para1):
                os.makedirs("/home/AMPfinder/public_html/test/weka-3-8-1/weka_temp/"+para1)

            src_file_list = glob(src_dir)

            for srcfile in src_file_list:
                mycopyfile(srcfile, dst_dir)


            os.rename(dst_dir + os.path.basename(src_dir) , dst_dir +  "/" + para2 + ".fasta")
            print(dst_dir + os.path.basename(src_dir))
            print(dst_dir  + para2 + ".fasta")

            cmd = "python /home/AMPfinder/public_html/test/arff.py " + para1 + " " + para2 + " < " + dst_dir + para2 +".fasta"
            print(cmd)
            os.system(cmd)

            parameter1 = para1 + "/" + para2
            parameter2 = "all"

            cmd2 = "php /home/AMPfinder/public_html/test/mk.php " + parameter1 + " " + parameter2
            os.system(cmd2)

            o_f_path, o_f_name = os.path.split(os.path.abspath(self.output_file))
            src_dir = '/home/AMPfinder/public_html/test/temp/' + para1 + '/' + para2 +".anti_finish.txt"
            dst_dir = o_f_path + "/"
            src_file_list = glob(src_dir)
            for srcfile in src_file_list:
                mycopyfile(srcfile, dst_dir)

            os.rename(dst_dir + para2 +".anti_finish.txt" , dst_dir + "/" + "anti_finish.txt")

            src_dir = '/home/AMPfinder/public_html/test/temp/' + para1 + '/' + para2 + ".nonanti_finish.txt"
            dst_dir = o_f_path + "/"
            src_file_list = glob(src_dir)
            for srcfile in src_file_list:
                mycopyfile(srcfile, dst_dir)

            os.rename(dst_dir + para2 + ".nonanti_finish.txt", dst_dir + "/" + "nonanti_finish.txt")

        if self.input_type == "contig":
            print("222222222222222222222222222222211213132422")
            print("enter contig")
            para1 = time.strftime("%Y%m%d", time.localtime())
            para2 = time.strftime("%H%M%S", time.localtime())
            o_f_path, o_f_name = os.path.split(os.path.abspath(self.input_sequence))
            src_dir = o_f_path + "/" + o_f_name +".temp.contig.fsa"
            dst_dir = '/home/AMPfinder/public_html/test/temp/' + para1 + '/'
            print(src_dir)
            print(dst_dir)
            if not os.path.isdir(dst_dir):
                os.makedirs(dst_dir)

            if not os.path.isdir("/home/AMPfinder/public_html/test/weka-3-8-1/weka_temp/" + para1):
                os.makedirs("/home/AMPfinder/public_html/test/weka-3-8-1/weka_temp/" + para1)

            src_file_list = glob(src_dir)

            for srcfile in src_file_list:
                mycopyfile(srcfile, dst_dir)

            os.rename(dst_dir + os.path.basename(src_dir), dst_dir + "/" + para2 + ".fasta")
            print(dst_dir + os.path.basename(src_dir))
            print(dst_dir + para2 + ".fasta")

            cmd = "python /home/AMPfinder/public_html/test/arff.py " + para1 + " " + para2 + " < " + dst_dir + para2 + ".fasta"
            print(cmd)
            os.system(cmd)

            parameter1 = para1 + "/" + para2
            parameter2 = "all"

            cmd2 = "php /home/AMPfinder/public_html/test/mk.php " + parameter1 + " " + parameter2
            os.system(cmd2)

            o_f_path, o_f_name = os.path.split(os.path.abspath(self.output_file))
            src_dir = '/home/AMPfinder/public_html/test/temp/' + para1 + '/' + para2 + ".anti_finish.txt"
            dst_dir = o_f_path + "/"
            src_file_list = glob(src_dir)
            for srcfile in src_file_list:
                mycopyfile(srcfile, dst_dir)

            os.rename(dst_dir + para2 + ".anti_finish.txt", dst_dir + "/" + "anti_finish.txt")

            src_dir = '/home/AMPfinder/public_html/test/temp/' + para1 + '/' + para2 + ".nonanti_finish.txt"
            dst_dir = o_f_path + "/"
            src_file_list = glob(src_dir)
            for srcfile in src_file_list:
                mycopyfile(srcfile, dst_dir)

            os.rename(dst_dir + para2 + ".nonanti_finish.txt", dst_dir + "/" + "nonanti_finish.txt")




    def makejsonfile(self):
        list = []
        obj={

        }
        i = 0
        with open(self.output_file, 'r') as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        obj["queryID"] = blast_record.query
                        obj["subjectID"] = alignment.hit_id
                        obj["query"] = hsp.query
                        obj["subject"] = hsp.sbjct
                        obj["identities"] = hsp.identities
                        obj["qstart"] = hsp.query_start
                        obj["qend"] = hsp.query_end
                        obj["sstart"] = hsp.sbjct_start
                        obj["send"] = hsp.sbjct_end
                        obj["alignlength"] = hsp.align_length
                        obj["bitscore"] = hsp.bits
                        obj["evalue"] = hsp.expect
                        obj["querylength"] = blast_record.query_letters
                        obj["identity"] = (hsp.identities / hsp.align_length) * 100
                        obj["hspcov"] = ((hsp.query_end - hsp.query_start+1)/blast_record.query_letters) * 100

                        evalue = obj["evalue"]
                        if(obj["evalue"]==0):
                            evalue ==1000
                        else:
                            evalue = (-1) * math.log10(evalue)

                        x = np.c_[obj["identities"],evalue,obj["bitscore"],obj["hspcov"]]
                        model = joblib.load("/home/AMPfinder/public_html/pythoncode/myAMP_advanced1/BlastModel_DIAMOND_final_1223.pkl")
                        y = model.predict(x)

                        if obj["identity"]==100:
                            obj["level"] = "perfact"

                        if obj["identity"]!=100 and  y==1:
                            obj["level"] = "strict"

                        if obj["identity"]!=100 and  y==0:
                            obj["level"] = "loose"


                        list.append(obj)
                        obj = {

                        }

            print(list)
            o_f_path, o_f_name = os.path.split(os.path.abspath(self.output_file))
            print(o_f_path+'/record.json')
            with open(o_f_path+'/record.json', 'w') as f:
                json.dump(list, f,sort_keys=True, indent=4, separators=(',', ': '))







    def clean_directory(self, directory, basename_output_file):
        """Cleans files in directory."""
        logger.info(directory)
        files = glob.glob(os.path.join(directory, "*"))
        for f in files:
            if os.path.basename(self.input_sequence) + ".temp" in f and os.path.isfile(f):
                self.remove_file(f)
            if os.path.basename(self.input_sequence) + ".fai" in f and os.path.isfile(f):
                self.remove_file(f)
        # if os.path.basename(f)[:3] == "tmp" in f and os.path.isfile(f) and ".temp." in f:
        #	self.remove_file(f)
        # if ".temp.directory" in f and os.path.isdir(f):
        #	logger.info("Removed directory: {}".format(f))
        #	shutil.rmtree(f)

    def remove_file(self, f):
        """Removes file."""
        if os.path.exists(f):
            try:
                logger.info("Removed file: {}".format(f))
                os.remove(f)
            except Exception as e:
                raise e
        else:
            logger.warning("Missing file: {}".format(f))



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
        tax_class = self.tax_class
        if self.aligner == "diamond":
            diamond_obj = Diamond(self.input_sequence, output_file = output, tax_class = tax_class , num_threads=self.threads)
            diamond_obj.run()
        else:
            blast_obj = Blast(input_file= file_name,  output_file = output, tax_class = tax_class , num_threads=self.threads)
            blast_obj.run()


    def process_contig(self):
        """Process nuclotide sequence(s)."""
        file_name = os.path.basename(self.input_sequence)
        output = self.output_file
        tax_class = self.tax_class
        print("OKOKOKOKOKOK")
        orf_obj = ORF(input_file=self.input_sequence)
        orf_obj.contig_to_orf()
        contig_fsa_file = os.path.join(self.working_directory, "{}.temp.contig.fsa".format(file_name))

        try:
            if os.stat(contig_fsa_file).st_size > 0:
                if self.aligner == "diamond":
                    diamond_obj = Diamond(input_file=contig_fsa_file, output_file = output, tax_class = tax_class,num_threads=self.threads)
                    diamond_obj.run()
                else:
                    print(contig_fsa_file)
                    blast_obj = Blast(input_file=contig_fsa_file, output_file = output, tax_class = tax_class, num_threads=self.threads)
                    blast_obj.run()
            else:
                self.write_stub_output_file()
        except Exception as e:
            pass
        else:
            pass




    def output(self):
        pass






