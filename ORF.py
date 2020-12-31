
import logging
import tempfile, time, fileinput, math, multiprocessing, shutil
from settings import os, SeqIO

class ORF(object):
	"""Class to find open reading frames from nucleotide sequence."""
	def __init__(self,input_file, clean=True, working_directory=None, low_quality=False):
		"""Creates ORF object for finding open reading frames."""
		self.input_file = input_file
		self.clean = clean
		self.working_directory = working_directory
		self.low_quality = low_quality

	def __repr__(self):
		"""Returns ORF class full object."""
		return "ORF({}".format(self.__dict__)

	def contig_to_orf(self):
		self.orf_prodigal()

	def min_max_sequence_length(self):
		"""Returns minimum and maximun sequence length in multi-fasta inputs"""
		sequences = []
		for record in SeqIO.parse(self.input_file, "fasta"):
			sequences.append(len(record.seq))
		return min(sequences), max(sequences), len(sequences)

	def orf_prodigal(self):
		"""Runs PRODIGAL to find open reading frames."""
		quality = "-n -p single"

		minimum_sequence_length, maximum_sequence_length, number_of_sequences = self.min_max_sequence_length()

		if self.low_quality == True or minimum_sequence_length < 20000:
			quality = "-p meta"

		filename = os.path.basename(self.input_file)

		stdout = "2> /dev/null"
		print("enter orf")
		print(filename)
		o_f_path, o_f_name = os.path.split(os.path.abspath(self.input_file))
		print(o_f_path)
		cmd = "prodigal -q -m -a {trans_file} -i {input_file} -o  {output_file} -d {nuc_file} -s {potential_genes} {quality} {stdout}" \
		.format(
				trans_file= o_f_path + "/{}.temp.contig.fsa".format(filename),
				input_file=self.input_file,
				output_file=o_f_path + "/{}.temp.draft".format(filename),
				quality=quality,
				stdout=stdout,
				nuc_file=o_f_path + "/{}.temp.contigToORF.fsa".format(filename),
				potential_genes= o_f_path + "/{}.temp.potentialGenes".format(filename)
			)

		# logger.debug(cmd)
		os.system(cmd)

		print("OK")

		# format the contig file headers to remove space
		#format_fasta_headers(working_directory+"/"+filename+".contig.fsa")

		if self.clean == True:
			os.remove(o_f_path + "/{}.temp.draft".format(filename))





