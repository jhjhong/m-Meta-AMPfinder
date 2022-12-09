
from settings import *
from AMP import*
from progress.bar import Bar
import argparse

class Main(object):
    def __init__(self, api=False):
        # """
        self.cpu_count = os.cpu_count()

    def main(self):
        parser = self.main_args()
        args = parser.parse_args(sys.argv[2:])
        self.main_run(args)


    def main_args(self):
        parser = argparse.ArgumentParser(prog="python ampfinder",
                                         description="{} - {} - main".format(APP_NAME, SOFTWARE_VERSION))
        parser.add_argument('-i', '--input_sequence', dest="input_sequence", required=True, \
                            help='input file must be in FASTA (contig and protein) format! e.g myFile.fasta')
        parser.add_argument('-o', '--output_file', dest="output_file", required=True,
                            help="output folder and base filename")
        parser.add_argument('-t', '--input_type', dest="input_type",
                            type=str.lower,
                            default="contig", choices=['contig', 'protein'],
                            required=False,
                            help='specify data input type (default = contig)')
        parser.add_argument('-a', '--alignment_tool', dest="aligner",
                            type=str.upper,
                            choices=['DIAMOND', 'BLAST'],
                            default="BLAST",
                            help="specify alignment tool (default = BLAST)")
        parser.add_argument('-n', '--num_threads', dest="threads", type=int,
                            default=self.cpu_count,
                            help="number of threads (CPUs) to use in the BLAST search (default={})".format(
                                self.cpu_count))

        parser.add_argument('-v', '--version', action='version', version="{}".format(SOFTWARE_VERSION),
                            help="prints software version number")
        return parser

    def main_run(self, args):
        amp_obj = AMP(**vars(args))
        amp_obj.run()

if __name__ == '__main__':
    m = MainBase()
    m.main()

