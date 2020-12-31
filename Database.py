import argparse
import json
import os
import re
from Base import *
from settings import *


class Database():
    """Class to create BLAST databases from a card.json file."""

    def __init__(self,args):
        """Creates Database object."""
        self.db = PATH
        self.args = args
        self.stdout = "2>&1 >> /dev/null"  # "2> /dev/null"

    def __repr__(self):
        """Returns Database class full object."""
        return "Database({}".format(self.__dict__)

    def run(self):
        self.build_databases()

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
        with open(os.path.join(working_directory, self.args.input_database), 'r') as jfile:
            print(os.path.join(working_directory, "ampfinder.json"))
            data = json.load(jfile)

        """
        write database fasta (FASTA format)
        """
        if os.path.exists("database_Euk.fasta"):
            os.remove("database_Euk.fasta")

        if os.path.exists("database_Fgl.fasta"):
            os.remove("database_Fgl.fasta")

        if os.path.exists("database_Bac.fasta"):
            os.remove("database_Bac.fasta")

        if os.path.exists("database_Arc.fasta"):
            os.remove("database_Arc.fasta")

        if os.path.exists("database_Vir.fasta"):
            os.remove("database_Vir.fasta")

        if os.path.exists("database_all.fasta"):
            os.remove("database_all.fasta")

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

            Taxclass = data[i]['dbAMP_Taxclass']
            if("Euk" in Taxclass):
                with open(os.path.join(working_directory, "database_Euk.fasta"), 'a') as fout:
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

            if ("Fgl" in Taxclass):
                with open(os.path.join(working_directory, "database_Fgl.fasta"), 'a') as fout:
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

            if ("Bac" in Taxclass):
                with open(os.path.join(working_directory, "database_Bac.fasta"), 'a') as fout:
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

            if ("Arc" in Taxclass):
                with open(os.path.join(working_directory, "database_Arc.fasta"), 'a') as fout:
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

            if ("Vir" in Taxclass):
                with open(os.path.join(working_directory, "database_Vir.fasta"), 'a') as fout:
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

            with open(os.path.join(working_directory, "database_all.fasta"), 'a') as fout:
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
        if os.path.isfile(os.path.join(self.db, "database_Euk.fasta")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_Euk.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Euk.db.phr")) == True and os.path.exists(
            os.path.join(self.db, "protein_Euk.db.pin")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Euk.db.psq")) == True:
            pass
        else:
            os.system('makeblastdb -in {} -dbtype prot -out {} '.format(os.path.join(self.db, "database_Euk.fasta"),
                                                                                os.path.join(self.db, "protein_Euk.db"),
                                                                                ))

        if os.path.isfile(os.path.join(self.db, "database_Fgl.fasta")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_Fgl.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Fgl.db.phr")) == True and os.path.exists(
            os.path.join(self.db, "protein_Fgl.db.pin")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Fgl.db.psq")) == True:
            pass
        else:
            os.system('makeblastdb -in {} -dbtype prot -out {} '.format(os.path.join(self.db, "database_Fgl.fasta"),
                                                                                os.path.join(self.db, "protein_Fgl.db"),
                                                                                ))

        if os.path.isfile(os.path.join(self.db, "database_Bac.fasta")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_Bac.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Bac.db.phr")) == True and os.path.exists(
            os.path.join(self.db, "protein_Bac.db.pin")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Bac.db.psq")) == True:
            pass
        else:
            os.system('makeblastdb -in {} -dbtype prot -out {} '.format(os.path.join(self.db, "database_Bac.fasta"),
                                                                                os.path.join(self.db, "protein_Bac.db"),
                                                                                ))

        if os.path.isfile(os.path.join(self.db, "database_Vir.fasta")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_Vir.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Vir.db.phr")) == True and os.path.exists(
            os.path.join(self.db, "protein_Vir.db.pin")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Vir.db.psq")) == True:
            pass
        else:
            os.system('makeblastdb -in {} -dbtype prot -out {} '.format(os.path.join(self.db, "database_Vir.fasta"),
                                                                                os.path.join(self.db, "protein_Vir.db"),
                                                                                ))

        if os.path.isfile(os.path.join(self.db, "database_all.fasta")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_all.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_all.db.phr")) == True and os.path.exists(
            os.path.join(self.db, "protein_all.db.pin")) == True \
                and os.path.exists(os.path.join(self.db, "protein_all.db.psq")) == True:
            pass
        else:
            os.system('makeblastdb -in {} -dbtype prot -out {} '.format(os.path.join(self.db, "database_all.fasta"),
                                                                                os.path.join(self.db, "protein_all.db"),
                                                                                ))


    def make_diamond_database(self):
        """Build DIAMOND database from a FASTA file."""
        if os.path.isfile(os.path.join(self.db, "proteindb_Euk.fsa")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_Euk.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Euk.db.dmnd")) == True:
            logger.info("diamond DB exists")
            pass
        else:
            os.system('diamond makedb --quiet --in {} --db {}'.format(os.path.join(self.db, "database_Euk.fasta"),
                                                                               os.path.join(self.db, "protein_Euk.db"),
                                                                               ))
        if os.path.isfile(os.path.join(self.db, "proteindb_Fgl.fsa")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_Fgl.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Fgl.db.dmnd")) == True:
            logger.info("diamond DB exists")
            pass
        else:
            os.system('diamond makedb --quiet --in {} --db {}'.format(os.path.join(self.db, "database_Fgl.fasta"),
                                                                               os.path.join(self.db, "protein_Fgl.db"),
                                                                               ))

        if os.path.isfile(os.path.join(self.db, "proteindb_Bac.fsa")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_Bac.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Bac.db.dmnd")) == True:
            logger.info("diamond DB exists")
            pass
        else:
            os.system('diamond makedb --quiet --in {} --db {}'.format(os.path.join(self.db, "database_Bac.fasta"),
                                                                               os.path.join(self.db, "protein_Bac.db"),
                                                                               ))

        if os.path.isfile(os.path.join(self.db, "proteindb_Vir.fsa")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_Vir.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_Vir.db.dmnd")) == True:
            logger.info("diamond DB exists")
            pass
        else:
            os.system('diamond makedb --quiet --in {} --db {}'.format(os.path.join(self.db, "database_Vir.fasta"),
                                                                               os.path.join(self.db, "protein_Vir.db"),
                                                                               ))

        if os.path.isfile(os.path.join(self.db, "proteindb_all.fsa")) == True and os.path.exists(
                os.path.join(self.db, "proteindb_all.fsa")) == True \
                and os.path.exists(os.path.join(self.db, "protein_all.db.dmnd")) == True:
            logger.info("diamond DB exists")
            pass
        else:
            os.system('diamond makedb --quiet --in {} --db {}'.format(os.path.join(self.db, "database_all.fasta"),
                                                                               os.path.join(self.db, "protein_all.db"),
                                                                               ))