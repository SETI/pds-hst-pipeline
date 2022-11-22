##########################################################################################
# hst_general/__init__.py
##########################################################################################
import os
# default start and end date of observation in query mast constraints
START_DATE = (1900, 1, 1)
END_DATE = (2025, 1, 1)

# default number of connection retry when connecting to mast failed.
RETRY = 1

# file directories
HST_DIR = {"staging": os.environ["HST_STAGING"],
           "pipeline": os.environ["HST_PIPELINE"],
           "bundle": os.environ["HST_BUNDLES"]}
HST_STAGING_DIR = os.environ["HST_STAGING"]
HST_PIPELINE_DIR = os.environ["HST_PIPELINE"]
HST_BUNDLE_DIR = os.environ["HST_BUNDLES"]
