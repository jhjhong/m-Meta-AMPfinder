import argparse
import json
import os
import re
from script.Base import *
from script.settings import *


class Database():
    def __init__(self):
        self.build_databases()
        self.stdout = "2>&1 >> /dev/null"  # "2> /dev/null"

    def build_databases(self):
        # self.write_fasta_from_json()
        self.make_blast_database()
        self.make_diamond_database()


    def make_blast_database(self):
        # Build BLAST database from a FASTA file.
        if os.path.isfile(os.path.join(data_path, "proteindb.fsa")) == True \
            and os.path.exists(os.path.join(path, "proteindb_all.fsa")) == True \
            and os.path.exists(os.path.join(path, "protein_all.db.phr")) == True \
            and os.path.exists(os.path.join(path, "protein_all.db.pin")) == True \
            and os.path.exists(os.path.join(path, "protein_all.db.psq")) == True:
            logger.info("blast DB exists")
            pass
        else:
            logger.info("create blast DB.")
            os.system('makeblastdb -in {} -dbtype prot -out {} {stdout}'.format(os.path.join(data_path, "proteindb.fsa"), os.path.join(path, "protein_all.db"),stdout=self.stdout))


    def make_diamond_database(self):
        # Build DIAMOND database from a FASTA file.
        if os.path.isfile(os.path.join(path, "proteindb_all.fsa")) == True \
            and os.path.exists(os.path.join(path, "proteindb_all.fsa")) == True \
            and os.path.exists(os.path.join(path, "protein_all.db.dmnd")) == True:
            logger.info("diamond DB exists")
            pass
        else:
            logger.info("create diamond DB.")
            os.system('diamond makedb --quiet --in {} --db {} {stdout}'.format(os.path.join(self.db, "proteindb.fsa"),os.path.join(self.db, "protein_all.db"),stdout=self.stdout))
