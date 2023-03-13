
import logging
import tempfile, time, fileinput, math, multiprocessing, shutil
from Bio import SeqIO
from script.settings import *

class PyORF(object):
    # find open reading frames using Pyrodigal.
    def __init__(self, input_file, clean=True, output_dir=None, low_quality=False, num_threads=16):
        self.input_file = input_file
        self.clean = clean
        self.output_dir = output_dir
        self.low_quality = low_quality

    def __repr__(self):
        return "ORF({}".format(self.__dict__)

    def min_max_sequence_length(self):
    	sequences = []
    	for record in SeqIO.parse(self.input_file, "fasta"):
    		sequences.append(len(record.seq))
    	return min(sequences), max(sequences), len(sequences)

    def filter_smorfs(self):
        long_sequences = []

        cmd = "sed 's/*//g' {input_file} > {output_file}" \
        .format(
                input_file = os.path.join(self.output_dir, "temp.contig.fsa"),
                output_file = os.path.join(self.output_dir, "temp.contig.fsa1")
            )
        # print(cmd)
        os.system(cmd)

        input_file = os.path.join(self.output_dir, "temp.contig.fsa1")
        output_file = os.path.join(self.output_dir, "temp.contig.L100.fsa")
        output_handle = open(output_file, "w")

        try:
            for record in SeqIO.parse(input_file, "fasta") :
                if len(record.seq) >= 10 and len(record.seq) <= 100 :
                    # Add this record to our list
                    long_sequences.append(record)

            SeqIO.write(long_sequences, output_file, "fasta")
            logger.info("filter done")
            
        except Exception as e:
            logger.exception("failed to write orf file")
        else:
            pass

    def run(self):
        logger.info("predict ORF")
        minimum_sequence_length, maximum_sequence_length, number_of_sequences = self.min_max_sequence_length()

        quality = "-n -p single"
        if self.low_quality == True or minimum_sequence_length < 20000:
            quality = "-p meta"

        filename = os.path.basename(self.input_file)
        trans_file = "temp.contig.fsa"
        output_file = "temp.draft"
        nuc_file = "temp.contigToORF.fsa"
        potential_genes = "temp.potentialGenes"

        if not os.path.exists(self.output_dir):
            logger.info("Making temp.orf folder {}".format(self.output_dir))
            os.makedirs(self.output_dir)

        stdout = "2> /dev/null"

        cmd = "pyrodigal -m -a {trans_file} -i {input_file} -o {output_file} -d {nuc_file} -s {potential_genes} {quality} {stdout}" \
        .format(
                trans_file=os.path.join(self.output_dir, trans_file),
                input_file=self.input_file,
                output_file=os.path.join(self.output_dir, output_file),
                quality=quality,
                stdout=stdout,
                nuc_file=os.path.join(self.output_dir, nuc_file),
                potential_genes=os.path.join(self.output_dir, potential_genes)
            )

        logger.debug(cmd)
        os.system(cmd)

        # filter 10 <= sequence length <= 100 for next classifier prediction
        self.filter_smorfs()


        if self.clean == True and os.path.exists(os.path.join(self.output_dir, output_file)):
            os.remove("{output_file}".format(output_file=os.path.join(self.output_dir, output_file)))

        # move final file
        original = os.path.join(self.output_dir, "temp.contig.fsa")
        original_L100 = os.path.join(self.output_dir, "temp.contig.L100.fsa")
        target = os.path.join(os.path.abspath(os.path.join(self.output_dir, os.path.pardir)), "final.smorfs.fsa")
        target_L100 = os.path.join(os.path.abspath(os.path.join(self.output_dir, os.path.pardir)), "final.smorfs.L100.fsa")
        os.system("ln -fs {final} {soft_link}".format(final = original, soft_link = target))
        os.system("ln -fs {final} {soft_link}".format(final = original_L100, soft_link = target_L100))