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
        print("I'm sorry, but something is wrong.")
        print("There is no __file__ variable. Please contact the author.")
        sys.exit()

# ====================================================================================
# FILEPATHS
# ====================================================================================

script_path = determine_path()

path = os.path.join(script_path, "_db/")
data_path = os.path.join(script_path, "_data/")

APP_NAME="AMP Identifier"
SOFTWARE_VERSION = "0.0.1"
SOFTWARE_SUMMARY = 'amp'

PATH = os.getcwd()
