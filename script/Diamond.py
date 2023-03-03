from script.settings import *

class Diamond(object):

    def __init__(self,input_file, output_dir=None, program = 'blastp', num_threads=32):
        self.input_file = input_file
        self.output_dir = output_dir
        self.db = path
        self.data_path = data_path

        self.program = program
        self.num_threads = num_threads
        self.index_chunks = 1
        self.block_size = 1
        self.outfmt = 6

    def __repr__(self):
        return "Diamond({}".format(self.__dict__)

    def run(self):
        # Runs DIAMOND algorithm.
        logger.info("Runs DIAMOND algorithm")

        if not os.path.exists(self.output_dir):
            logger.info("Making temp.alignment folder {}".format(self.output_dir))
            os.makedirs(self.output_dir)

        cmd = "diamond {program} --in {in_ref} --db {db} --query {input} --outfmt {outfmt} --out {output_file} --threads {num_threads} --index-chunks {index_chunks} --block-size {block_size} --salltitles --quiet --more-sensitive" \
            .format(
                program=self.program,
                in_ref=os.path.join(self.data_path,"proteindb.fsa"),
                db=os.path.join(self.db,"protein.db"),
                input=self.input_file,
                output_file=os.path.join(self.output_dir,"temp.alignment.txt"),
                num_threads=self.num_threads,
                index_chunks=self.index_chunks,
                block_size=self.block_size,
                outfmt=self.outfmt
            )

        logger.debug(cmd)
        os.system(cmd)

        # move final file
        original = os.path.join(self.output_dir, "temp.alignment.txt")
        target = os.path.join(os.path.abspath(os.path.join(self.output_dir, os.path.pardir)), "final.alignment.txt")
        os.system("ln -fs {final} {soft_link}".format(final = original, soft_link = target))