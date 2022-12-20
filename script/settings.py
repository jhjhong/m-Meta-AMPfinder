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
        return os.path.dirname(os.path.abspath(os.path.join(root, os.path.pardir)))
    except:
        sys.exit()

# ====================================================================================
# FILEPATHS
# ====================================================================================
script_path = determine_path()
db_path = os.path.join(script_path, "db/")

path = os.path.join(db_path, "_db/")
data_path = os.path.join(db_path, "_data/")

# ====================================================================================
# LOGGING CONFIG
# ====================================================================================
level = logging.WARNING
logger = logging.getLogger(__name__)
logger.setLevel(level)

# detailed log
# formatter = logging.Formatter('%(levelname)s %(asctime)s : (%(filename)s::%(funcName)s::%(lineno)d) : %(message)s')
# basic log
formatter = logging.Formatter('%(levelname)s %(asctime)s : (%(filename)s::%(funcName)s::%(lineno)d) : %(message)s')

stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)

logger.addHandler(stream_handler)

# ====================================================================================
# APP INFO
# ====================================================================================
APP_NAME="m(Meta)ACPfinder"
SOFTWARE_VERSION = "v0.0.1"
SOFTWARE_SUMMARY = 'Use the m(Meta)ACPfinder to predict candidate ACP from protein or nucleotide'

PATH = os.getcwd()
