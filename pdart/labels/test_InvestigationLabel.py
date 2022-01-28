import unittest

from typing import List, Tuple

from pdart.citations import Citation_Information
from pdart.db.BundleDB import create_bundle_db_in_memory
from pdart.labels.InvestigationLabel import make_investigation_label
from pdart.labels.utils import assert_golden_file_equal

_BUNDLE_LIDVID = "urn:nasa:pds:hst_09059::1.0"
_COLLECTION_LIDVID = "urn:nasa:pds:hst_09059:data_acs_raw::1.2"
_FITS_PRODUCT_LIDVID = "urn:nasa:pds:hst_09059:data_acs_raw:j6gp01lzq_raw::1.2"


class Test_InvestigationLabel(unittest.TestCase):
    def setUp(self) -> None:
        self.db = create_bundle_db_in_memory()
        self.db.create_tables()
        self.db.create_bundle(_BUNDLE_LIDVID)
        self.db.create_other_collection(_COLLECTION_LIDVID, _BUNDLE_LIDVID)
        self.db.create_fits_product(_FITS_PRODUCT_LIDVID, _COLLECTION_LIDVID)
        self.info = Citation_Information.create_test_citation_information()

        # Create start/stop time in db for testing purpose
        self.db.update_fits_product_time(
            _FITS_PRODUCT_LIDVID,
            "2005-01-19T14:58:56Z",
            "2005-01-19T15:41:05Z",
        )

        # Create target identifications db for testing purpose
        target_id = "09059_1"
        target_identifications: List[Tuple] = [
            (
                "762 Pulcova",
                [
                    "(762) Pulcova",
                    "(762) 1913SQ",
                    "1913SQ",
                    "(762) 1952 QM1",
                    "1952 QM1",
                    "Pulcova",
                    "Minor Planet 762",
                    "NAIF ID 2000762",
                ],
                "Asteroid",
                [],
                "urn:nasa:pds:context:target:asteroid.762_pulcova",
            )
        ]
        self.db.add_record_to_target_identification_db(
            target_id, target_identifications
        )

        # Create wavelength range in db for testing purpose
        self.db.update_wavelength_range(
            _FITS_PRODUCT_LIDVID,
            ["Ultraviolet", "Visible", "Near Infrared", "Infrared"],
        )

        context_products_for_test = [
            (
                "urn:nasa:pds:context:investigation:individual.hst_09059::1.0",
                "investigation",
            ),
            (
                "urn:nasa:pds:context:instrument_host:spacecraft.hst::1.0",
                "instrument_host",
            ),
            ("urn:nasa:pds:context:instrument:hst.acs::1.0", "instrument"),
            ("urn:nasa:pds:context:target:asteroid.762_pulcova::1.0", "target"),
        ]
        for id, ref_type in context_products_for_test:
            self.db.create_context_product(id, ref_type)

    def test_make_bundle_label(self) -> None:
        label = make_investigation_label(self.db, _BUNDLE_LIDVID, self.info, True, True)
        assert_golden_file_equal(self, "test_InvestigationLabel.golden.xml", label)
