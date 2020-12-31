import argparse
import json
import os
import re
from settings import *


class Database(object):
    """Class to create BLAST databases from a card.json file."""

    def __init__(self):
        """Creates Database object."""
        self.db = PATH
        self.stdout = "2>&1 >> /dev/null"  # "2> /dev/null"

    def __repr__(self):
        """Returns Database class full object."""
        return "Database({}".format(self.__dict__)


    def build_databases(self):
        """Build BLAST and DIAMOND databases."""
        self.write_fasta_from_json()
        self.make_blast_database()
        self.make_diamond_database()

    def write_fasta_from_json(self):
        working_directory = os.getcwd()
        """
        reading ampfinder.json
        """
        with open(os.path.join(working_directory, "ampfinder.json"), 'r') as jfile:
            print(os.path.join(working_directory, "ampfinder.json"))
            data = json.load(jfile)

        """
        write database fasta (FASTA format)
        """
        with open(os.path.join(working_directory, "database.fasta"), 'w') as fout:
            for i in range(0, len(data)):
                Taxonomy = data[i]['dbAMP_Taxonomy']
                numList = [m.start() for m in re.finditer(",", Taxonomy)]
                if (len(numList) != 0):
                    Taxonomy = Taxonomy[numList[-1] + 1:]
                    numList = [m.start() for m in re.finditer("&&", Taxonomy)]
                    if (len(numList) != 0):
                        Taxonomy = Taxonomy[0:numList[0]]

                numList = [m.start() for m in re.finditer("&&", Taxonomy)]
                if (len(numList) != 0):
                    Taxonomy = Taxonomy[0:numList[0]]
                Taxonomy = Taxonomy.replace('.', '')
                try:
                    header = ("dbAMP_ID:{}|dbAMP_Taxonomy:{}".format(
                        data[i]['dbAMP_ID'],
                        Taxonomy
                    ))
                    sequence = data[i]['Sequence']
                    fout.write(">{}\n".format(header))
                    fout.write("{}\n".format(sequence))

                except Exception as e:
                    print("No dbAMP {}. Omitting this and keep running.")



    def make_blast_database(self):
        """Build BLAST database from a FASTA file."""
        if os.path.isfile(os.path.join(self.db, "database.fasta")) == True and os.path.exists(
                os.path.join(self.db, "proteindb.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein.db.phr")) == True and os.path.exists(
            os.path.join(self.db, "protein.db.pin")) == True \
                and os.path.exists(os.path.join(self.db, "protein.db.psq")) == True:
            pass
        else:
            os.system('makeblastdb -in {} -dbtype prot -out {} '.format(os.path.join(self.db, "database.fasta"),
                                                                                os.path.join(self.db, "protein.db"),
                                                                                ))

    def make_diamond_database(self):
        """Build DIAMOND database from a FASTA file."""
        if os.path.isfile(os.path.join(self.db, "proteindb.fsa")) == True and os.path.exists(
                os.path.join(self.db, "proteindb.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein.db.dmnd")) == True:
            logger.info("diamond DB exists")
            pass
        else:
            os.system('diamond makedb --quiet --in {} --db {}'.format(os.path.join(self.db, "database.fasta"),
                                                                               os.path.join(self.db, "protein.db"),
                                                                               ))



db_obj = Database()
db_obj.build_databases()
