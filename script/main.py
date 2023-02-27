from script.settings import *
from script.load import *
from script.AMP import*

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

                run mACPfinder on paired-end reads:
                macpfinder main -1 example_seqs/R1.fq.gz -2 example_seqs/R2.fq.gz -t read --output out_metag --outtag example_metag

                

                removes databases and temporary files:
                macpfinder clean
                
                For more information,please read the docs: https://github.com/jhjhong/mACPfinder
               '''
                # run mACPfinder on peptides:  
                # macpfinder peptides --fasta example_seqs/expep.faa.gz --output out_peptides
                
                # run mACPfinder on contigs:
                # macpfinder contigs --fasta example_seqs/excontigs.fna.gz --output out_contigs
                
                # run mACPfinder on paired-end reads:
                # macpfinder reads -1 example_seqs/R1.fq.gz -2 example_seqs/R2.fq.gz --output out_metag --outtag example_metag

        parser = argparse.ArgumentParser(prog="macpfinder", description='{} - {}'.format(APP_NAME, SOFTWARE_VERSION), epilog=SOFTWARE_SUMMARY, usage=USAGE)
        # parser.add_argument('command', choices=['peptides', 'contigs', 'reads', 'load', 'auto_load', 'clean'], help='Subcommand to run')
        parser.add_argument('command', choices=['main', 'load', 'auto_load', 'clean'], help='Subcommand to run')

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
                            default="BLAST",
                            help="specify alignment tool (default = BLAST)")
        parser.add_argument('-v', '--version', action='version', version="{}".format(SOFTWARE_VERSION),
                            help="show mACPfinder software version number")
        parser.add_argument('--debug', dest="debug", action="store_true", help="debug mode")
        return parser

    def main_run(self, args):
        print(vars(args))
        # args = ['gtdbtk', 'classify_wf', '--genome_dir', genome_test_dir,
        #             '--out_dir', output_dir, '--cpus', str(options.cpus), '-f']
        #     self.logger.info('Command: {}'.format(' '.join(args)))
        amp_obj = AMP(**vars(args))
        amp_obj.run()


    # def contigs(self):
    #     parser = self.contigs_args()
    #     if len(sys.argv) == 2:
    #         parser.print_help(file=sys.stderr)
    #         sys.exit(1)
    #     args = parser.parse_args(sys.argv[2:])
    #     self.contigs_run(args)

    # def contigs_args(self):
    #     parser = argparse.ArgumentParser(prog="macpfinder contigs",
    #                                      description="{} - {} - main".format(APP_NAME, SOFTWARE_VERSION))
    #     parser.add_argument('-i', '--input_fasta', dest="input_sequence", required=True, \
    #                         help='input file must be in FASTA (contig and protein) format! e.g myFile.fasta')
    #     parser.add_argument('-o', '--output_dir', required=True,
    #                         help="path to the output directory")
    #     parser.add_argument('-a', '--alignment_tool',
    #                         type=str.upper,
    #                         choices=['DIAMOND', 'BLAST'],
    #                         default="BLAST",
    #                         help="specify alignment tool (default = BLAST)")
    #     parser.add_argument('-n', '--num_threads', type=int,
    #                         default=self.cpu_count,
    #                         help="number of threads (CPUs) to use in the BLAST search (default={})".format(
    #                             self.cpu_count))
    #     parser.add_argument('-v', '--version', action='version', version="{}".format(SOFTWARE_VERSION),
    #                         help="show mACPfinder software version number")
    #     return parser

    # def contigs_run(self, args):
    #     amp_obj = AMP(**vars(args))
    #     amp_obj.run()
    

    # def peptides(self):
    #     parser = self.peptides_args()
    #     if len(sys.argv) == 2:
    #         parser.print_help(file=sys.stderr)
    #         sys.exit(1)
    #     args = parser.parse_args(sys.argv[2:])
    #     self.peptides_run(args)

    # def peptides_args(self):
    #     parser = argparse.ArgumentParser(prog="macpfinder peptides",
    #                                      description="{} - {} - main".format(APP_NAME, SOFTWARE_VERSION))
    #     parser.add_argument('-i', '--input_fasta', dest="input_sequence", required=True, \
    #                         help='input file must be in FASTA (contig and protein) format! e.g myFile.fasta')
    #     parser.add_argument('-o', '--output_dir', dest="output_dir", required=True,
    #                         help="path to the output directory")
    #     parser.add_argument('-a', '--alignment_tool', dest="aligner",
    #                         type=str.upper,
    #                         choices=['DIAMOND', 'BLAST'],
    #                         default="BLAST",
    #                         help="specify alignment tool (default = BLAST)")
    #     parser.add_argument('-n', '--num_threads', dest="threads", type=int,
    #                         default=self.cpu_count,
    #                         help="number of threads (CPUs) to use in the BLAST search (default={})".format(
    #                             self.cpu_count))
    #     parser.add_argument('-v', '--version', action='version', version="{}".format(SOFTWARE_VERSION),
    #                         help="show mACPfinder software version number")
    #     return parser

    # def peptides_run(self, args):
    #     amp_obj = AMP(**vars(args))
    #     amp_obj.run()


if __name__ == '__main__':
    m = MainBase()
    m.main()

