import argparse
import json
import os
import re
from script.Base import *
from script.settings import *


class Database(object):
    def __init__(self):
        self.stdout = "2>&1 >> /dev/null"  # "2> /dev/null"
        self.build_databases()

    def build_databases(self):
        # self.write_fasta_from_json()
        self.make_blast_database()
        self.make_diamond_database()


    def make_blast_database(self):
        # Build BLAST database from a FASTA file.
        if os.path.isfile(os.path.join(data_path, "proteindb.fsa")) == True \
            and os.path.exists(os.path.join(data_path, "proteindb.fsa")) == True \
            and os.path.exists(os.path.join(path, "protein.db.phr")) == True \
            and os.path.exists(os.path.join(path, "protein.db.pin")) == True \
            and os.path.exists(os.path.join(path, "protein.db.psq")) == True:
            logger.info("blast DB exists")
            pass
        else:
            logger.info("create blast DB.")
            # os.system('makeblastdb -in {} -dbtype prot -out {}'.format(os.path.join(data_path, "proteindb.fsa"), os.path.join(path, "protein_all.db")))
            os.system('makeblastdb -in {in_file} -dbtype prot -out {out_file} {stdout}'.format(in_file = os.path.join(data_path, "proteindb.fsa"), out_file = os.path.join(path, "protein.db"), stdout = self.stdout))


    def make_diamond_database(self):
        # Build DIAMOND database from a FASTA file.
        if os.path.isfile(os.path.join(data_path, "proteindb.fsa")) == True \
            and os.path.exists(os.path.join(data_path, "proteindb.fsa")) == True \
            and os.path.exists(os.path.join(path, "protein.db.dmnd")) == True:
            logger.info("diamond DB exists")
            pass
        else:
            logger.info("create diamond DB.")
            # os.system('diamond makedb --quiet --in {} --db {}'.format(os.path.join(data_path, "proteindb.fsa"),os.path.join(path, "protein_all.db")))
            os.system('diamond makedb --quiet --in {in_file} --db {out_file} {stdout}'.format(in_file = os.path.join(data_path, "proteindb.fsa"), out_file = os.path.join(path, "protein.db"), stdout = self.stdout))
