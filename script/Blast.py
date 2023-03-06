from script.settings import *


class Blast(object):

    def __init__(self, input_file, output_dir=None, program='blastp', num_threads=32):
        self.input_file = input_file
        self.output_dir = output_dir
        self.db = path

        self.program = program
        self.num_threads = num_threads
        self.outfmt = 6

    def __repr__(self):
        return "Blast({}".format(self.__dict__)

    def run(self):
        # Runs BLAST algorithm.
        logger.info("Runs BLAST algorithm")

        if not os.path.exists(self.output_dir):
            logger.info("Making temp.alignment folder {}".format(self.output_dir))
            os.makedirs(self.output_dir)

        cmd = "{program} -query {input} -db {db} -num_threads {num_threads} -outfmt {outfmt} -out {output_file}" \
            .format(
                program=self.program,
                input=self.input_file,
                db=os.path.join(self.db, "protein.db"),
                num_threads=self.num_threads,
                outfmt=self.outfmt,
                output_file=os.path.join(self.output_dir, "temp.alignment.txt")
            )
        logger.debug(cmd)
        os.system(cmd)

        # move final file
        original = os.path.join(self.output_dir, "temp.alignment.txt")
        target = os.path.join(os.path.abspath(os.path.join(self.output_dir, os.path.pardir)), "final.alignment.txt")
        os.system("ln -fs {final} {soft_link}".format(final = original, soft_link = target))