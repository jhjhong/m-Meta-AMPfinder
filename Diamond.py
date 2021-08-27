from settings import *

class Diamond(object):
<<<<<<< HEAD
	def __init__(self,input_file, output_file=None, program = 'blastp', num_threads=32):
=======
	"""Class to create Diamond object and align for protein and translated DNA searches."""
	def __init__(self,input_file, output_file=None, tax_class = None, program = 'blastp', num_threads=32):
		"""Creates Diamond object for running DIAMOND algorithm."""
>>>>>>> eacbca3da2cb9e7278c934f260583853a33684a4
		self.input_file = input_file
		if output_file == None:
			f_path, f_name = os.path.split(input_file)
			self.output_file = os.path.join(f_path,"{}.blastRes.xml".format(f_name))
		else:
			self.output_file = output_file
		self.db = PATH

		self.program = program
<<<<<<< HEAD
		self.num_threads = num_threads
		self.index_chunks = 1
		self.block_size = 1
		self.outfmt = 6

	def __repr__(self):
		return "Diamond({}".format(self.__dict__)

	def run(self):
		db_name = "protein_all.db"
=======
		self.tax_class = tax_class
		self.num_threads = num_threads
		self.index_chunks = 1
		self.block_size = 1
		self.outfmt = 5

	def __repr__(self):
		"""Returns Diamond class full object."""
		return "Diamond({}".format(self.__dict__)

	def run(self):
		"""Runs DIAMOND algorithm.""" 
		# TODO:: for versions of diamond > 0.8.36 use --xml-blord-format
		db_name = "protein_" + self.tax_class + ".db"
		print(db_name)
		print(os.getcwd())
>>>>>>> eacbca3da2cb9e7278c934f260583853a33684a4
		cmd = ('diamond {program} --in {in_ref} --db {db} \
				   --query {input} --outfmt {outfmt} --out {output_file}  \
				   --threads {num_threads}  --index-chunks {index_chunks} \
				   --block-size {block_size}  \
<<<<<<< HEAD
				   --salltitles  --quiet --more-sensitive' \
					.format(
						program=self.program,
						in_ref=os.path.join(self.db,"proteindb.fsa"),
						db=os.path.join(self.db,db_name),
=======
				   ' \
					.format(
						program=self.program,
						in_ref=os.path.join(self.db,"proteindb.fsa"),
						db=os.path.join(self.db+"/pythoncode/myAMP_advanced1",db_name),
>>>>>>> eacbca3da2cb9e7278c934f260583853a33684a4
						input=self.input_file,
						output_file=self.output_file,
						num_threads=self.num_threads,
						index_chunks=self.index_chunks,
						block_size=self.block_size,
						outfmt=self.outfmt
					)
				)

<<<<<<< HEAD
		os.system(cmd)
=======
		# logger.debug(cmd)
		os.system(cmd)
		print(cmd)
>>>>>>> eacbca3da2cb9e7278c934f260583853a33684a4
	
