from settings import *
from AMP import *
from progress.bar import Bar
from Database import *
import argparse


class makedatabase(object):
    def __init__(self):
<<<<<<< HEAD
        self.cpu_count = os.cpu_count()

=======
        # """
        self.cpu_count = os.cpu_count()


>>>>>>> eacbca3da2cb9e7278c934f260583853a33684a4
    def database(self):
        parser = self.database_args()
        args = parser.parse_args(sys.argv[2:])
        self.databse_run(args)

    def database_args(self):
        parser = argparse.ArgumentParser(prog="amp database",
                                         description="{} - {} - Database".format(APP_NAME, SOFTWARE_VERSION))
        parser.add_argument('-i', '--input_database', dest="input_database", required=True, \
                            help='input file must be in json format! e.g database.json')
        return parser

    def databse_run(self, args):
        database_obj = Database(args)
        database_obj.run()


if __name__ == '__main__':
<<<<<<< HEAD
    m = makedatabase()
    m.database()
=======
    time_start = time.time()
    m = makedatabase()
    m.database()
    time_end = time.time()
    print('totally cost', time_end - time_start)
>>>>>>> eacbca3da2cb9e7278c934f260583853a33684a4
