import unittest

from fs.path import join

from pdart.pds4.Archives import get_any_archive
from pdart.pds4.Bundle import *
from pdart.pds4.Collection import *
from pdart.pds4.HstFilename import *
from pdart.pds4.Product import *
from pdart.pds4.LID import LID


class TestProduct(unittest.TestCase):
    def test_init(self):
        # type: () -> None
        arch = get_any_archive()
        lid = LID('urn:nasa:pds:bundle:collection:uproduct')
        p = Product(arch, lid)
        self.assertEquals(lid, p.lid)

        # check that creation of bundle fails with collection LID
        lid = LID('urn:nasa:pds:bundle:collection')
        try:
            Collection(arch, lid)
            self.assertTrue(False)
        except Exception:
            pass

    def test_absolute_filepath(self):
        # type: () -> None
        arch = get_any_archive()
        for p in arch.products():
            if not p.is_document_product():
                lid = p.lid
                visit = p.visit()

                actual_fp = p.absolute_filepath()
                expected_fps = [join(arch.root,
                                     lid.bundle_id,
                                     lid.collection_id,
                                     ('visit_%s' % visit),
                                     lid.product_id + ext)
                                for ext in Product.FILE_EXTS]

                assert actual_fp in expected_fps
                # We only check the first product
                return

    def test_bundle(self):
        # type: () -> None
        arch = get_any_archive()
        lid = LID('urn:nasa:pds:bundle:collection:uproduct')
        p = Product(arch, lid)
        self.assertEquals(Bundle(arch, LID('urn:nasa:pds:bundle')),
                          p.bundle())

    def test_collection(self):
        # type: () -> None
        arch = get_any_archive()
        lid = LID('urn:nasa:pds:bundle:collection:uproduct')
        p = Product(arch, lid)
        self.assertEquals(Collection(arch,
                                     LID('urn:nasa:pds:bundle:collection')),
                          p.collection())

    def test_document_product(self):
        # type: () -> None
        arch = get_any_archive()
        lid = LID('urn:nasa:pds:bundle:document:product')
        p = Product(arch, lid)
        self.assertEquals(p.absolute_filepath(),
                          join(arch.root, 'bundle',
                               'document', 'product'))
        with self.assertRaises(Exception):
            p.visit_filepath()
        with self.assertRaises(Exception):
            p.visit()
        self.assertTrue(p.absolute_filepath_is_directory())

    def test_non_document_product(self):
        # type: () -> None
        arch = get_any_archive()
        lid = LID('urn:nasa:pds:bundle:collection:product')
        p = Product(arch, lid)
        self.assertFalse(p.absolute_filepath_is_directory())
