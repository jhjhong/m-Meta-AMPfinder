from script.settings import *
from script.AMP import *
from progress.bar import Bar
from script.Database import *
import argparse



def database(self):
        parser = self.database_args()
        args = parser.parse_args(sys.argv[2:])
        self.databse_run(args)

def create_parser():
    parser = argparse.ArgumentParser(prog="macpfinder load", description="{} - {} - Load".format(APP_NAME, SOFTWARE_VERSION))
    parser.add_argument('-i', '--input_json', required=True, help='input file must be in json format! e.g ampfinder_210729.json')
    return parser

def run():
	parser = create_parser()
	args = parser.parse_args()
	main(args)

if __name__ == "__main__":
	run()

"""
class makedatabase(object):
    def __init__(self):
        self.cpu_count = os.cpu_count()

    def database(self):
        parser = self.database_args()
        args = parser.parse_args(sys.argv[2:])
        self.databse_run(args)

    def database_args(self):
        parser = argparse.ArgumentParser(prog="macpfinder load",
                                         description="{} - {} - Database".format(APP_NAME, SOFTWARE_VERSION))
        parser.add_argument('-i', '--input_database', dest="input_database", required=True, \
                            help='input file must be in json format! e.g database.json')
        return parser

    def databse_run(self, args):
        database_obj = Database(args)
        database_obj.run()


if __name__ == '__main__':
    m = makedatabase()
    m.database()
"""