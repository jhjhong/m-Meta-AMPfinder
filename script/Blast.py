from settings import *


class Blast(object):

    def __init__(self, input_file, output_file=None, program='blastp', num_threads=32):
        self.input_file = input_file
        if output_file == None:
            f_path, f_name = os.path.split(input_file)
            self.output_file = os.path.join(f_path, "{}.blastRes.xml".format(f_name))
        else:
            self.output_file = output_file
        self.db = PATH

        self.program = program
        self.num_threads = num_threads
        self.outfmt = 6

    def __repr__(self):
        return "Blast({}".format(self.__dict__)

    def run(self):
        db_name = "protein_all.db"
        os.system('{program} -query {input} -db {path} \
					-num_threads {num_threads} -outfmt {outfmt} -out {output_file}' \
            .format(
            program=self.program,
            num_threads=self.num_threads,
            outfmt=self.outfmt,
            input=self.input_file,
            path=os.path.join(self.db, db_name),
            output_file=self.output_file
        )
        )