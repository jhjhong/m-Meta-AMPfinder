import os
import sys
import logging
import json

from Bio.Blast import NCBIXML
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

# ====================================================================================
# LOGGING CONFIG
# ====================================================================================
level = logging.WARNING
logger = logging.getLogger(__name__)
logger.setLevel(level)

# detailed log
# formatter = logging.Formatter('%(levelname)s %(asctime)s : (%(filename)s::%(funcName)s::%(lineno)d) : %(message)s')
# basic log
formatter = logging.Formatter('%(levelname)s %(asctime)s : %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

logger.addHandler(stream_handler)

APP_NAME="Meta-ACPfinder"
SOFTWARE_VERSION = "0.0.1"
SOFTWARE_SUMMARY = 'Use the Meta-ACPfinder to predict candidate ACP from protein or nucleotide'

PATH = os.getcwd()
