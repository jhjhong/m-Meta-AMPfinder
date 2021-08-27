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
        sys.exit()

script_path = determine_path()

path = os.path.join(script_path, "_db/")
data_path = os.path.join(script_path, "_data/")

APP_NAME="AMPfinder"
SOFTWARE_VERSION = "0.0.3"
SOFTWARE_SUMMARY = 'ampfinder'

PATH = os.getcwd()
