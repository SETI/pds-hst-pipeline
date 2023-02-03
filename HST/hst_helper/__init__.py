##########################################################################################
# hst_helper/__init__.py
##########################################################################################
import os
from collections import defaultdict

# default start and end date of observation in query mast constraints
START_DATE = (1900, 1, 1)
END_DATE = (2025, 1, 1)

# default number of connection retry when connecting to mast failed.
RETRY = 1

# file directories
HST_DIR = {'staging': os.environ['HST_STAGING'],
           'pipeline': os.environ['HST_PIPELINE'],
           'bundles': os.environ['HST_BUNDLES']}

# suffixes for proposal files
DOCUMENT_EXT = ('apt', 'pdf', 'pro', 'prop')
DOCUMENT_EXT_FOR_CITATION_INFO = ('apt', 'pro')

# File names
PROGRAM_INFO_FILE = 'program-info.txt'
PRODUCTS_FILE = 'products.txt'
TRL_CHECKSUMS_FILE = 'trl_checksums.txt'

# Instrument ids dictionary, keyed by propoposal id and store the list of instrument ids
INST_ID_DICT = defaultdict(set)
