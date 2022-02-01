import os
import sys
from pdart.astroquery.astroquery import MastSlice, CustomizedQueryMastSlice, ProductSet
from pdart.pipeline.directories import DevDirectories
from pdart.pipeline.suffix_info import (  # type: ignore
    ACCEPTED_SUFFIXES,
    PART_OF_ACCEPTED_SUFFIXES,
)

# Download from Mast
if __name__ == "__main__":
    TWD = os.environ["TMP_WORKING_DIR"]
    assert len(sys.argv) == 2, sys.argv
    proposal_id = int(sys.argv[1])
    dirs = DevDirectories(TWD + "/shm_spt_from_mast")
    working_dir = dirs.working_dir(proposal_id)
    slice = CustomizedQueryMastSlice((1900, 1, 1), (2025, 1, 1), proposal_id)
    product_set = slice.to_product_set(proposal_id, True)
    product_set.download(working_dir)
