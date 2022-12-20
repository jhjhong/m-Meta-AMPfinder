
from script.settings import *
import script.load
# from script.AMP import*
from progress.bar import Bar
import argparse

class Main(object):
    def __init__(self):
        self.cpu_count = os.cpu_count()
        USAGE='''%(prog)s <command> [<args>]

            commands are:

                Database:
                ---------------------------------------------------------------------------------------
                auto_load   Automatically loads blast database, and installed mACPfinder databases
                load        Loads mACPfinder database
                install     Install related mACPfinder alignment database
                clean       Removes current databases and temporary files

                
                m(Meta)ACPfinder:
                ---------------------------------------------------------------------------------------
                main        Runs mACPfinder application

            Examples:
            
                parse blast database:
                macpfinder load --input_json ampfinder_210729.json

                install blast/diamond database:
                macpfinder install

                run Macrel on peptides:  
                macrel peptides --fasta example_seqs/expep.faa.gz --output out_peptides
                
                run Macrel on contigs:
                macrel contigs --fasta example_seqs/excontigs.fna.gz --output out_contigs
                
                run Macrel on paired-end reads:
                macrel reads -1 example_seqs/R1.fq.gz -2 example_seqs/R2.fq.gz --output out_metag --outtag example_metag
                
                run Macrel to get abundance profiles: 
                macrel abundance -1 example_seqs/R1.fq.gz --fasta example_seqs/ref.faa.gz --output out_abundance --outtag example_abundance

                removes databases and temporary files:
                macpfinder clean
                
                For more information,please read the docs: https://macrel.readthedocs.io/en/latest/
               '''

        parser = argparse.ArgumentParser(prog="macpfinder", description='{} - {}'.format(APP_NAME, SOFTWARE_VERSION), epilog=SOFTWARE_SUMMARY, usage=USAGE)
        parser.add_argument('command', choices=['main', 'load', 'auto_load', 'install', 'clean'], help='Subcommand to run')

        args=parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            logger.info("Unrecognized command: {}".format(args.command))
            exit("Error: Unrecognized command: {}".format(args.command))
        getattr(self, args.command)()


    def load(self):
        parser = self.load_args()
        args = parser.parse_args(sys.argv[2:])
        self.load_run(args)

    def load_args(self):
        parser = script.load.create_parser()
        return parser

    def load_run(self, args):
        script.load.main(args)
    
    # def install(self):
    #     parser = self.install_args()
    #     args = parser.parse_args(sys.argv[2:])
    #     print(self.install_run(args))

    # def install_args(self):
    #     parser = argparse.ArgumentParser(prog="macpfinder install", description="{} - {} - Install".format(APP_NAME, SOFTWARE_VERSION))
    #     parser.add_argument('-v','--version',action='store_true', required=True, help = "prints data version number")
    #     parser.add_argument('--local', dest="local_database", action='store_true', help="use local database (default: uses database in executable directory)")
    #     parser.add_argument('--all', action='store_true', help="data version number used for `rgi bwt` and `rgi main` (default: rgi main)")
    #     return parser

    # def install_run(self, args):
    #     obj = Install(args.version)
    #     obj.run()

"""
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
"""

if __name__ == '__main__':
    m = MainBase()
    m.main()

