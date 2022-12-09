from script.settings import *

class Diamond(object):
	def __init__(self,input_file, output_file=None, program = 'blastp', num_threads=32):
		self.input_file = input_file
		if output_file == None:
			f_path, f_name = os.path.split(input_file)
			self.output_file = os.path.join(f_path,"{}.blastRes.xml".format(f_name))
		else:
			self.output_file = output_file
		self.db = PATH

		self.program = program
		self.num_threads = num_threads
		self.index_chunks = 1
		self.block_size = 1
		self.outfmt = 6

	def __repr__(self):
		return "Diamond({}".format(self.__dict__)

	def run(self):
		db_name = "protein_all.db"
		cmd = ('diamond {program} --in {in_ref} --db {db} \
				   --query {input} --outfmt {outfmt} --out {output_file}  \
				   --threads {num_threads}  --index-chunks {index_chunks} \
				   --block-size {block_size}  \
				   --salltitles  --quiet --more-sensitive' \
					.format(
						program=self.program,
						in_ref=os.path.join(self.db,"proteindb.fsa"),
						db=os.path.join(self.db,db_name),
						input=self.input_file,
						output_file=self.output_file,
						num_threads=self.num_threads,
						index_chunks=self.index_chunks,
						block_size=self.block_size,
						outfmt=self.outfmt
					)
				)

		os.system(cmd)
	
