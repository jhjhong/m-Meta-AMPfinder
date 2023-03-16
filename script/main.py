from script.settings import *
from script.AMP import *
import script.load
import script.clean


from progress.bar import Bar
import argparse

class Main(object):
    def __init__(self):
        self.cpu_count = os.cpu_count()
        USAGE='''%(prog)s <command> [<args>]

            commands are:

                Database:
                ---------------------------------------------------------------------------------------
                load        Loads mACPfinder database
                clean       Removes current databases and temporary files

                
                m(Meta)ACPfinder:
                ---------------------------------------------------------------------------------------
                main        Runs mACPfinder application

            Examples:
            
                parse blast database:
                macpfinder load --input_json ampfinder_210729.json

                run mACPfinder on paired-end reads:
                macpfinder main --short1 examples/example_read1.fastq --short2 examples/example_read2.fastq --input_type read --output_dir test

                run mACPfinder on genome contigs:
                macpfinder main --input_fasta examples/example_contig.fsa --input_type contig --output_dir test

                run mACPfinder on protein peptides:
                macpfinder main --input_fasta examples/example_prot.fsa --input_type peptide --output_dir test


                removes databases and temporary files:
                macpfinder clean
                
                For more information,please read the docs: https://github.com/jhjhong/mACPfinder
               '''

        parser = argparse.ArgumentParser(prog="macpfinder", description='{} - {}'.format(APP_NAME, SOFTWARE_VERSION), epilog=SOFTWARE_SUMMARY, usage=USAGE)
        parser.add_argument('command', choices=['main', 'load', 'clean'], help='Subcommand to run')

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


    def clean(self):
        parser = self.clean_args()
        args = parser.parse_args(sys.argv[2:])
        self.clean_run(args)

    def clean_args(self):
        parser = script.clean.create_parser()
        return parser

    def clean_run(self, args):
        script.clean.main(args)
    
# ====================================================================================
# main 
# ====================================================================================

    def main(self):
        parser = self.main_args()
        if len(sys.argv) == 2:
            parser.print_help(file=sys.stderr)
            sys.exit(1)
        args = parser.parse_args(sys.argv[2:])
        self.main_run(args)

    def main_args(self):
        parser = argparse.ArgumentParser(prog="macpfinder main",
                                         description="{} - {} - main".format(APP_NAME, SOFTWARE_VERSION))
        parser.add_argument('-1', '--short1', required=False, \
                            help='FASTQ file (metagenomics reads required) of first short reads in each pair')
        parser.add_argument('-2', '--short2', required=False, \
                            help='FASTQ file (metagenomics reads required) of second short reads in each pair')
        parser.add_argument('-i', '--input_fasta', dest="input_sequence", required=False, \
                            help='input file must be in FASTA (contig and peptide required) format!')
        parser.add_argument('-o', '--output_dir', required=True,
                            help="path to the output directory (required)")
        parser.add_argument('-t', '--input_type', dest="input_type",
                            type=str.lower,
                            choices=['read', 'contig', 'peptide'],
                            required=True,
                            help='specify data input type')
        parser.add_argument('-n', '--num_threads', type=int, dest="threads",
                            default=self.cpu_count,
                            help="number of threads (CPUs) to use in the BLAST search (default={})".format(self.cpu_count))
        parser.add_argument('--assembler', dest="assembler",
                            type=str.lower,
                            choices=['megahit', 'metaspades'],
                            default="megahit",
                            help="specify assembler tool (default = megahit)")
        parser.add_argument('--alignment_tool', dest="aligner",
                            type=str.upper,
                            choices=['DIAMOND', 'BLAST'],
                            default="DIAMOND",
                            help="specify alignment tool (default = BLAST)")
        parser.add_argument('-v', '--version', action='version', version="{}".format(SOFTWARE_VERSION),
                            help="show mACPfinder software version number")
        parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
        return parser

    def main_run(self, args):
        logger.info(vars(args))
        amp_obj = AMP(**vars(args))
        amp_obj.run()

if __name__ == '__main__':
    m = MainBase()
    m.main()

