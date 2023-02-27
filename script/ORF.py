
import logging
import tempfile, time, fileinput, math, multiprocessing, shutil
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

    # def min_max_sequence_length(self):
    # 	sequences = []
    # 	for record in SeqIO.parse(self.input_file, "fasta"):
    # 		sequences.append(len(record.seq))
    # 	return min(sequences), max(sequences), len(sequences)

    def run(self):
        logger.info("predict ORF")

        # check quality if it need meta mode (low_quality).
        records = list(SeqIO.parse(self.input_file, "fasta"))
        sequences = [ str(record.seq) for record in records ]
        minimum_sequence_length = min((len(seq) for seq in sequences), default=0)
        # minimum_sequence_length, maximum_sequence_length, number_of_sequences = self.min_max_sequence_length()

        quality = "-n -p single"
        if self.low_quality == True or minimum_sequence_length < 20000:
            quality = "-p meta"

        filename = os.path.basename(self.input_file)
        trans_file="{}.temp.contig.fsa".format(filename)
        output_file="{}.temp.draft".format(filename)
        nuc_file="{}.temp.contigToORF.fsa".format(filename)
        potential_genes="{}.temp.potentialGenes".format(filename)

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
                potential_genes=os.path.join(self.output_dir, potential_genes),
            )

        print(cmd)
        os.system(cmd)

        if self.clean == True and os.path.exists(os.path.join(self.output_dir, output_file)):
            os.remove("{output_file}".format(output_file=os.path.join(self.output_dir, output_file)))





