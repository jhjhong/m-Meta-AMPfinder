
from script.settings import *
from script.AMP import*
from progress.bar import Bar
import argparse

class Main(object):
    def __init__(self, api=False):
        # """
        self.cpu_count = os.cpu_count()
        USAGE='''%(prog)s <command> [<args>]
            commands are:
               ---------------------------------------------------------------------------------------
               Database
               ---------------------------------------------------------------------------------------
               auto_load Automatically loads CARD database, annotations and k-mer database
               load      Loads CARD database, annotations and k-mer database
               clean     Removes BLAST databases and temporary files
               database  Information on installed card database
               galaxy    Galaxy project wrapper
               ---------------------------------------------------------------------------------------
               Genomic
               ---------------------------------------------------------------------------------------
               main     Runs rgi application
               tab      Creates a Tab-delimited from rgi results
               parser   Creates categorical JSON files RGI wheel visualization
               heatmap  Heatmap for multiple analysis
               ---------------------------------------------------------------------------------------
               Metagenomic
               ---------------------------------------------------------------------------------------
               bwt                   Align reads to CARD and in silico predicted allelic variants (beta)
               ---------------------------------------------------------------------------------------
               Baits validation
               ---------------------------------------------------------------------------------------
               tm                    Baits Melting Temperature
               ---------------------------------------------------------------------------------------
               Annotations
               ---------------------------------------------------------------------------------------
               card_annotation       Create fasta files with annotations from card.json
               wildcard_annotation   Create fasta files with annotations from variants
               baits_annotation      Create fasta files with annotations from baits (experimental)
               remove_duplicates     Removes duplicate sequences (experimental)
               ---------------------------------------------------------------------------------------
               Pathogen of origin
               ---------------------------------------------------------------------------------------
               kmer_build            Build AMR specific k-mers database used for pathogen of origin (beta)
               kmer_query            Query sequences against AMR k-mers database to predict pathogen of origin (beta)
               '''

        parser = argparse.ArgumentParser(prog="acpfinder", description='{} - {}'.format(APP_NAME, SOFTWARE_VERSION), epilog=SOFTWARE_SUMMARY, usage=USAGE)
        parser.add_argument('command', choices=['main', 'tab', 'parser', 'load', 'auto_load'], help='Subcommand to run')

        if api == False:
            args=parser.parse_args(sys.argv[1:2])
            if not hasattr(self, args.command):
                logger.info("Unrecognized command: {}".format(args.command))
                exit("Error: Unrecognized command: {}".format(args.command))
            getattr(self, args.command)()

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

