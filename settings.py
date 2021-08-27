import os
import sys
import logging
import json

from Bio.Blast import NCBIXML
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio import SeqIO

def determine_path():
    try:
        root = __file__
        if os.path.islink(root):
            root = os.path.realpath(root)
        return os.path.dirname(os.path.abspath(root))
    except:
<<<<<<< HEAD
        sys.exit()

=======
        print("I'm sorry, but something is wrong.")
        print("There is no __file__ variable. Please contact the author.")
        sys.exit()

# ====================================================================================
# FILEPATHS
# ====================================================================================

>>>>>>> eacbca3da2cb9e7278c934f260583853a33684a4
script_path = determine_path()

path = os.path.join(script_path, "_db/")
data_path = os.path.join(script_path, "_data/")

<<<<<<< HEAD
APP_NAME="AMPfinder"
SOFTWARE_VERSION = "0.0.3"
SOFTWARE_SUMMARY = 'ampfinder'
=======
APP_NAME="AMP Identifier"
SOFTWARE_VERSION = "0.0.1"
SOFTWARE_SUMMARY = 'amp'
>>>>>>> eacbca3da2cb9e7278c934f260583853a33684a4

PATH = os.getcwd()
