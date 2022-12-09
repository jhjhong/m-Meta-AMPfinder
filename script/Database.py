import argparse
import json
import os
import re
from Base import *
from settings import *


class Database():
    def __init__(self,args):
        self.db = PATH
        self.args = args
        self.stdout = "2>&1 >> /dev/null"  # "2> /dev/null"

    def __repr__(self):
        return "Database({}".format(self.__dict__)

    def run(self):
        self.build_databases()

    def build_databases(self):
        self.write_fasta_from_json()
        self.make_blast_database()
        self.make_diamond_database()

    def write_fasta_from_json(self):
        working_directory = os.getcwd()
        with open(os.path.join(working_directory, self.args.input_database), 'r') as jfile:
            data = json.load(jfile)
            with open(os.path.join(working_directory, "database.fasta"), 'a') as fout:
                for i in range(0, len(data)):
                    header = (data[i]['dbAMP_ID'])
                    sequence = data[i]['Seq']
                    fout.write(">{}\n".format(header))
                    fout.write("{}\n".format(sequence))


    def make_blast_database(self):
        if os.path.isfile(os.path.join(self.db, "database.fasta")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_all.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_all.db.phr")) == True and os.path.exists(
            os.path.join(self.db, "protein_all.db.pin")) == True \
                and os.path.exists(os.path.join(self.db, "protein_all.db.psq")) == True:
            pass
        else:
            os.system('makeblastdb -in {} -dbtype prot -out {} '.format(os.path.join(self.db, "database.fasta"),
                                                                                os.path.join(self.db, "protein_all.db"),
                                                                                ))


    def make_diamond_database(self):
        if os.path.isfile(os.path.join(self.db, "proteindb_all.fsa")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_all.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_all.db.dmnd")) == True:
            logger.info("diamond DB exists")
            pass
        else:
            os.system('diamond makedb --quiet --in {} --db {}'.format(os.path.join(self.db, "database.fasta"),
                                                                               os.path.join(self.db, "protein_all.db"),
                                                                               ))
