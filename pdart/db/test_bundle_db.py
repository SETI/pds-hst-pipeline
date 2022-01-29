import os
import tempfile
import unittest
from typing import Set, cast

from pdart.db.bundle_db import create_bundle_db_in_memory
from pdart.db.sql_alch_tables import (
    Base,
    Bundle,
    OtherCollection,
    switch_on_collection_subtype,
    CitationInfo,
)
from pdart.db.utils import file_md5
from pdart.pds4.LIDVID import LIDVID

_TABLES: Set[str] = {
    "associations",
    "bundles",
    "bundle_collection_links",
    "collections",
    "collection_product_links",
    "context_collections",
    "document_collections",
    "schema_collections",
    "other_collections",
    "products",
    "browse_products",
    "context_products",
    "schema_products",
    "document_products",
    "fits_products",
    "files",
    "fits_files",
    "bad_fits_files",
    "browse_files",
    "document_files",
    "hdus",
    "cards",
    "bundle_labels",
    "collection_labels",
    "collection_inventories",
    "product_labels",
    "proposal_info",
    "target_identification",
    "citation_info",
    "context_product_labels",
}


class TestBundleDB(unittest.TestCase):
    def setUp(self) -> None:
        self.db = create_bundle_db_in_memory()
        self.db.create_tables()
        (handle, filepath) = tempfile.mkstemp(suffix="_raw.fits", prefix="i")
        os.write(handle, os.urandom(32))
        os.close(handle)
        self.dummy_os_filepath = filepath

    def tearDown(self) -> None:
        os.remove(self.dummy_os_filepath)

    def test_bundle_db(self) -> None:
        db = self.db
        self.assertTrue(db.is_open())
        metadata = Base.metadata
        self.assertEqual(set(metadata.tables.keys()), _TABLES)
        db.close()
        self.assertFalse(db.is_open())

    ############################################################

    def _context_collection_exists(self, collection_lidvid: str) -> bool:
        if self.db.collection_exists(collection_lidvid):
            coll = self.db.get_collection(collection_lidvid)
            return switch_on_collection_subtype(coll, True, False, False, False)
        else:
            return False

    def _document_collection_exists(self, collection_lidvid: str) -> bool:
        if self.db.collection_exists(collection_lidvid):
            coll = self.db.get_collection(collection_lidvid)
            return switch_on_collection_subtype(coll, False, True, False, False)
        else:
            return False

    def _schema_collection_exists(self, collection_lidvid: str) -> bool:
        if self.db.collection_exists(collection_lidvid):
            coll = self.db.get_collection(collection_lidvid)
            return switch_on_collection_subtype(coll, False, False, True, False)
        else:
            return False

    def _other_collection_exists(self, collection_lidvid: str) -> bool:
        if self.db.collection_exists(collection_lidvid):
            coll = self.db.get_collection(collection_lidvid)
            return switch_on_collection_subtype(coll, False, False, False, True)
        else:
            return False

    ############################################################

    def test_create_bundle(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.assertTrue(
            self.db.session.query(Bundle).filter(Bundle.lidvid == bundle_lidvid).count()
            == 0
        )
        self.assertFalse(self.db.bundle_exists(bundle_lidvid))
        self.db.create_bundle(bundle_lidvid)
        self.assertTrue(
            self.db.session.query(Bundle).filter(Bundle.lidvid == bundle_lidvid).count()
            == 1
        )
        self.assertTrue(self.db.bundle_exists(bundle_lidvid))
        self.db.create_bundle(bundle_lidvid)
        self.assertTrue(self.db.bundle_exists(bundle_lidvid))

    def test_bundle_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.assertFalse(self.db.bundle_exists(bundle_lidvid))
        self.db.create_bundle(bundle_lidvid)
        self.assertTrue(self.db.bundle_exists(bundle_lidvid))

    def test_get_bundle(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        bundle = self.db.get_bundle(bundle_lidvid)
        self.assertEqual(bundle_lidvid, bundle.lidvid)
        self.assertEqual(99999, bundle.proposal_id)

    def test_get_bundle_collections(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        other_collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        doc_collection_lidvid = "urn:nasa:pds:hst_99999:document::1.1"
        ctxt_collection_lidvid = "urn:nasa:pds:hst_99999:context::1.1"
        sche_collection_lidvid = "urn:nasa:pds:hst_99999:schema::1.1"

        # test empty
        self.db.create_bundle(bundle_lidvid)
        self.assertFalse(self.db.get_bundle_collections(bundle_lidvid))

        # test inserting one
        self.db.create_other_collection(other_collection_lidvid, bundle_lidvid)
        self.assertEqual(
            {other_collection_lidvid},
            {c.lidvid for c in self.db.get_bundle_collections(bundle_lidvid)},
        )

        # test inserting a different kind
        self.db.create_document_collection(doc_collection_lidvid, bundle_lidvid)
        self.assertEqual(
            {other_collection_lidvid, doc_collection_lidvid},
            {c.lidvid for c in self.db.get_bundle_collections(bundle_lidvid)},
        )

        # test inserting a third different kind
        self.db.create_context_collection(ctxt_collection_lidvid, bundle_lidvid)
        self.assertEqual(
            {other_collection_lidvid, doc_collection_lidvid, ctxt_collection_lidvid},
            {c.lidvid for c in self.db.get_bundle_collections(bundle_lidvid)},
        )

        # test inserting a fourth kind
        self.db.create_schema_collection(sche_collection_lidvid, bundle_lidvid)
        self.assertEqual(
            {
                other_collection_lidvid,
                doc_collection_lidvid,
                ctxt_collection_lidvid,
                sche_collection_lidvid,
            },
            {c.lidvid for c in self.db.get_bundle_collections(bundle_lidvid)},
        )

    ############################################################

    def test_create_context_collection(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:context::1.8"
        self.assertFalse(self._context_collection_exists(collection_lidvid))

        self.db.create_context_collection(collection_lidvid, bundle_lidvid)
        self.assertTrue(self._context_collection_exists(collection_lidvid))

        self.db.create_context_collection(collection_lidvid, bundle_lidvid)
        self.assertTrue(self._context_collection_exists(collection_lidvid))

    def test_create_document_collection(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:document::1.8"
        self.assertFalse(self._document_collection_exists(collection_lidvid))

        self.db.create_document_collection(collection_lidvid, bundle_lidvid)
        self.assertTrue(self._document_collection_exists(collection_lidvid))

        self.db.create_document_collection(collection_lidvid, bundle_lidvid)
        self.assertTrue(self._document_collection_exists(collection_lidvid))

    def test_create_other_collection(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.assertFalse(self._other_collection_exists(collection_lidvid))

        self.db.create_other_collection(collection_lidvid, bundle_lidvid)
        self.assertTrue(self._other_collection_exists(collection_lidvid))

        self.db.create_other_collection(collection_lidvid, bundle_lidvid)
        self.assertTrue(self._other_collection_exists(collection_lidvid))

    def test_collection_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        self.db.create_bundle(bundle_lidvid)
        self.assertFalse(self.db.collection_exists(collection_lidvid))
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        self.assertTrue(self.db.collection_exists(collection_lidvid))

    def test_document_collection_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        other_collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        doc_collection_lidvid = "urn:nasa:pds:hst_99999:document::1.1"

        self.db.create_other_collection(other_collection_lidvid, bundle_lidvid)
        self.assertFalse(self._document_collection_exists(doc_collection_lidvid))

        self.db.create_document_collection(doc_collection_lidvid, bundle_lidvid)
        self.assertTrue(self._document_collection_exists(doc_collection_lidvid))

    def test_other_collection_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        other_collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        doc_collection_lidvid = "urn:nasa:pds:hst_99999:document::1.1"

        self.db.create_document_collection(doc_collection_lidvid, bundle_lidvid)
        self.assertFalse(self._other_collection_exists(other_collection_lidvid))

        self.db.create_other_collection(other_collection_lidvid, bundle_lidvid)
        self.assertTrue(self._other_collection_exists(other_collection_lidvid))

    def test_get_collection(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        collection = self.db.get_collection(collection_lidvid)
        self.assertTrue(isinstance(collection, OtherCollection))
        nd_collection = cast(OtherCollection, collection)
        self.assertEqual("acs", nd_collection.instrument)
        self.assertEqual("raw", nd_collection.suffix)

    def test_get_collection_products(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        p1_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:j6gp01lzq::1.8"
        p2_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:j6gp02lzq::1.8"
        p3_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:j6gp03lzq::1.8"

        self.db.create_bundle(bundle_lidvid)
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)
        self.assertFalse(self.db.get_collection_products(collection_lidvid))

        self.db.create_fits_product(p1_lidvid, collection_lidvid)

        self.assertEqual(
            {p1_lidvid},
            {p.lidvid for p in self.db.get_collection_products(collection_lidvid)},
        )

        self.db.create_fits_product(p2_lidvid, collection_lidvid)

        self.assertEqual(
            {p1_lidvid, p2_lidvid},
            {p.lidvid for p in self.db.get_collection_products(collection_lidvid)},
        )

        self.db.create_fits_product(p3_lidvid, collection_lidvid)

        self.assertEqual(
            {p1_lidvid, p2_lidvid, p3_lidvid},
            {p.lidvid for p in self.db.get_collection_products(collection_lidvid)},
        )

    ############################################################

    def test_create_browse_product(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        data_collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(data_collection_lidvid, bundle_lidvid)

        browse_collection_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw::1.8"
        self.db.create_other_collection(browse_collection_lidvid, bundle_lidvid)

        fits_product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::8.1"
        browse_product_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw:p::8.1"

        with self.assertRaises(Exception):
            self.db.create_browse_product(
                browse_product_lidvid, fits_product_lidvid, browse_collection_lidvid
            )

        self.db.create_fits_product(fits_product_lidvid, data_collection_lidvid)
        self.assertFalse(self.db.browse_product_exists(browse_product_lidvid))

        self.db.create_browse_product(
            browse_product_lidvid, fits_product_lidvid, browse_collection_lidvid
        )
        self.assertTrue(self.db.browse_product_exists(browse_product_lidvid))

        self.db.create_browse_product(
            browse_product_lidvid, fits_product_lidvid, browse_collection_lidvid
        )
        self.assertTrue(self.db.browse_product_exists(browse_product_lidvid))

    def test_create_document_product(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:document::1.8"
        self.db.create_document_collection(collection_lidvid, bundle_lidvid)
        product_lidvid = "urn:nasa:pds:hst_99999:document:important_stuff::1.8"
        self.assertFalse(self.db.document_product_exists(product_lidvid))

        self.db.create_document_product(product_lidvid, collection_lidvid)
        self.assertTrue(self.db.document_product_exists(product_lidvid))

    def test_create_fits_product(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::8.1"
        self.assertFalse(self.db.fits_product_exists(product_lidvid))

        self.db.create_fits_product(product_lidvid, collection_lidvid)
        self.assertTrue(self.db.fits_product_exists(product_lidvid))

        self.db.create_fits_product(product_lidvid, collection_lidvid)
        self.assertTrue(self.db.fits_product_exists(product_lidvid))

    def test_product_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        data_coll_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        browse_coll_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw::1.8"
        doc_coll_lidvid = "urn:nasa:pds:hst_99999:document::1.8"
        self.db.create_other_collection(data_coll_lidvid, bundle_lidvid)
        self.db.create_other_collection(browse_coll_lidvid, bundle_lidvid)
        self.db.create_document_collection(doc_coll_lidvid, bundle_lidvid)

        fits_prod_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::1.8"
        browse_prod_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw:p::1.8"
        doc_prod_lidvid = "urn:nasa:pds:hst_99999:document:specs::1.8"

        self.assertFalse(self.db.product_exists(fits_prod_lidvid))
        self.assertFalse(self.db.product_exists(browse_prod_lidvid))
        self.assertFalse(self.db.product_exists(doc_prod_lidvid))

        self.db.create_fits_product(fits_prod_lidvid, data_coll_lidvid)

        self.assertTrue(self.db.product_exists(fits_prod_lidvid))
        self.assertFalse(self.db.product_exists(browse_prod_lidvid))
        self.assertFalse(self.db.product_exists(doc_prod_lidvid))

        self.db.create_browse_product(
            browse_prod_lidvid, fits_prod_lidvid, browse_coll_lidvid
        )

        self.assertTrue(self.db.product_exists(fits_prod_lidvid))
        self.assertTrue(self.db.product_exists(browse_prod_lidvid))
        self.assertFalse(self.db.product_exists(doc_prod_lidvid))

        self.db.create_document_product(doc_prod_lidvid, doc_coll_lidvid)

        self.assertTrue(self.db.product_exists(fits_prod_lidvid))
        self.assertTrue(self.db.product_exists(browse_prod_lidvid))
        self.assertTrue(self.db.product_exists(doc_prod_lidvid))

    def test_browse_product_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw::1.8"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)
        fits_prod_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::1.8"
        product_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw:p::1.8"

        self.assertFalse(self.db.browse_product_exists(product_lidvid))

        self.db.create_fits_product(fits_prod_lidvid, collection_lidvid)
        self.db.create_browse_product(
            product_lidvid, fits_prod_lidvid, collection_lidvid
        )
        self.assertTrue(self.db.browse_product_exists(product_lidvid))

    def test_context_product_exists(self) -> None:
        product_lidvid = "urn:nasa:pds:hst_99999:context:important_stuff::1.8"

        self.assertFalse(self.db.context_product_exists(product_lidvid))
        self.db.create_context_product(product_lidvid, "important")
        self.assertTrue(self.db.context_product_exists(product_lidvid))

    def test_document_product_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:document::1.8"
        self.db.create_document_collection(collection_lidvid, bundle_lidvid)
        product_lidvid = "urn:nasa:pds:hst_99999:document:important_stuff::1.8"

        self.assertFalse(self.db.document_product_exists(product_lidvid))
        self.db.create_document_product(product_lidvid, collection_lidvid)
        self.assertTrue(self.db.document_product_exists(product_lidvid))

    def test_fits_product_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)
        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::1.8"

        self.assertFalse(self.db.fits_product_exists(product_lidvid))
        self.db.create_fits_product(product_lidvid, collection_lidvid)
        self.assertTrue(self.db.fits_product_exists(product_lidvid))

    def test_get_product(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::8.1"
        self.db.create_fits_product(product_lidvid, collection_lidvid)

        product = self.db.get_product(product_lidvid)
        self.assertEqual(product_lidvid, product.lidvid)

    def test_get_product_files(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::8.1"
        self.db.create_fits_product(product_lidvid, collection_lidvid)

        self.assertFalse(self.db.get_product_files(product_lidvid))

        self.db.create_fits_file(
            self.dummy_os_filepath, "p_raw.fits", product_lidvid, 1
        )
        self.assertEqual(1, len(self.db.get_product_files(product_lidvid)))

    ############################################################

    def test_create_browse_file(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        data_collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(data_collection_lidvid, bundle_lidvid)

        browse_collection_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw::1.8"
        self.db.create_other_collection(browse_collection_lidvid, bundle_lidvid)

        fits_product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::8.1"
        self.db.create_fits_product(fits_product_lidvid, data_collection_lidvid)
        browse_product_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw:p::8.1"
        self.db.create_browse_product(
            browse_product_lidvid, fits_product_lidvid, browse_collection_lidvid
        )

        basename = "file.jpg"
        self.assertFalse(self.db.browse_file_exists(basename, browse_product_lidvid))

        self.db.create_browse_file(
            self.dummy_os_filepath, basename, browse_product_lidvid, 1
        )
        self.assertTrue(self.db.browse_file_exists(basename, browse_product_lidvid))

        self.db.create_browse_file(
            self.dummy_os_filepath, basename, browse_product_lidvid, 1
        )
        self.assertTrue(self.db.browse_file_exists(basename, browse_product_lidvid))

    def test_create_document_file(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:document::1.8"
        self.db.create_document_collection(collection_lidvid, bundle_lidvid)

        product_lidvid = "urn:nasa:pds:hst_99999:document:p::8.1"
        self.db.create_document_product(product_lidvid, collection_lidvid)

        filename = "stuff.pdf"
        self.assertFalse(self.db.document_file_exists(filename, product_lidvid))

        self.db.create_document_file(self.dummy_os_filepath, filename, product_lidvid)
        self.assertTrue(self.db.document_file_exists(filename, product_lidvid))

    def test_create_fits_file(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::8.1"
        self.db.create_fits_product(product_lidvid, collection_lidvid)

        basename = "file.fits"
        self.assertFalse(self.db.fits_file_exists(basename, product_lidvid))

        self.db.create_fits_file(self.dummy_os_filepath, basename, product_lidvid, 1)
        self.assertTrue(self.db.fits_file_exists(basename, product_lidvid))

        self.db.create_fits_file(self.dummy_os_filepath, basename, product_lidvid, 1)
        self.assertTrue(self.db.fits_file_exists(basename, product_lidvid))

    def test_file_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        raw_collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        browse_collection_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw::1.8"
        doc_collection_lidvid = "urn:nasa:pds:hst_99999:document::1.8"
        self.db.create_other_collection(raw_collection_lidvid, bundle_lidvid)
        self.db.create_other_collection(browse_collection_lidvid, bundle_lidvid)
        self.db.create_document_collection(doc_collection_lidvid, bundle_lidvid)
        bad_fits_prod_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:x::1.8"
        fits_prod_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::1.8"
        browse_prod_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw:p::1.8"
        doc_prod_lidvid = "urn:nasa:pds:hst_99999:document:specs::1.8"
        self.db.create_fits_product(bad_fits_prod_lidvid, raw_collection_lidvid)
        self.db.create_fits_product(fits_prod_lidvid, raw_collection_lidvid)
        self.db.create_browse_product(
            browse_prod_lidvid, fits_prod_lidvid, browse_collection_lidvid
        )
        self.db.create_document_product(doc_prod_lidvid, doc_collection_lidvid)

        bad_fits_filename = "x_raw.fits"
        fits_filename = "p_raw.fits"
        browse_filename = "p_raw.jpg"
        doc_filename = "specs.pdf"
        doc_filename_2 = "specs.txt"

        self.assertFalse(self.db.file_exists(bad_fits_filename, bad_fits_prod_lidvid))
        self.assertFalse(self.db.file_exists(fits_filename, fits_prod_lidvid))
        self.assertFalse(self.db.file_exists(browse_filename, browse_prod_lidvid))
        self.assertFalse(self.db.file_exists(doc_filename, doc_prod_lidvid))
        self.assertFalse(self.db.file_exists(doc_filename_2, doc_prod_lidvid))

        self.db.create_bad_fits_file(
            self.dummy_os_filepath, bad_fits_filename, bad_fits_prod_lidvid, "kablooey!"
        )

        self.assertTrue(self.db.file_exists(bad_fits_filename, bad_fits_prod_lidvid))
        self.assertFalse(self.db.file_exists(fits_filename, fits_prod_lidvid))
        self.assertFalse(self.db.file_exists(browse_filename, browse_prod_lidvid))
        self.assertFalse(self.db.file_exists(doc_filename, doc_prod_lidvid))
        self.assertFalse(self.db.file_exists(doc_filename_2, doc_prod_lidvid))

        self.db.create_fits_file(
            self.dummy_os_filepath, fits_filename, fits_prod_lidvid, 666
        )

        self.assertTrue(self.db.file_exists(bad_fits_filename, bad_fits_prod_lidvid))
        self.assertTrue(self.db.file_exists(fits_filename, fits_prod_lidvid))
        self.assertFalse(self.db.file_exists(browse_filename, browse_prod_lidvid))
        self.assertFalse(self.db.file_exists(doc_filename, doc_prod_lidvid))
        self.assertFalse(self.db.file_exists(doc_filename_2, doc_prod_lidvid))

        self.db.create_browse_file(
            self.dummy_os_filepath, browse_filename, browse_prod_lidvid, 666
        )

        self.assertTrue(self.db.file_exists(bad_fits_filename, bad_fits_prod_lidvid))
        self.assertTrue(self.db.file_exists(fits_filename, fits_prod_lidvid))
        self.assertTrue(self.db.file_exists(browse_filename, browse_prod_lidvid))
        self.assertFalse(self.db.file_exists(doc_filename, doc_prod_lidvid))
        self.assertFalse(self.db.file_exists(doc_filename_2, doc_prod_lidvid))

        self.db.create_document_file(
            self.dummy_os_filepath, doc_filename, doc_prod_lidvid
        )

        self.assertTrue(self.db.file_exists(bad_fits_filename, bad_fits_prod_lidvid))
        self.assertTrue(self.db.file_exists(fits_filename, fits_prod_lidvid))
        self.assertTrue(self.db.file_exists(browse_filename, browse_prod_lidvid))
        self.assertTrue(self.db.file_exists(doc_filename, doc_prod_lidvid))
        self.assertFalse(self.db.file_exists(doc_filename_2, doc_prod_lidvid))

        self.db.create_document_file(
            self.dummy_os_filepath, doc_filename_2, doc_prod_lidvid
        )

        self.assertTrue(self.db.file_exists(bad_fits_filename, bad_fits_prod_lidvid))
        self.assertTrue(self.db.file_exists(fits_filename, fits_prod_lidvid))
        self.assertTrue(self.db.file_exists(browse_filename, browse_prod_lidvid))
        self.assertTrue(self.db.file_exists(doc_filename, doc_prod_lidvid))
        self.assertTrue(self.db.file_exists(doc_filename_2, doc_prod_lidvid))

    def test_bad_fits_file_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::8.1"
        self.db.create_fits_product(product_lidvid, collection_lidvid)

        filename = "p_raw.fits"
        self.assertFalse(self.db.bad_fits_file_exists(filename, product_lidvid))

        self.db.create_bad_fits_file(
            self.dummy_os_filepath, filename, product_lidvid, "kablooey!"
        )
        self.assertTrue(self.db.bad_fits_file_exists(filename, product_lidvid))

    def test_browse_file_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        raw_collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        collection_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw::1.8"
        self.db.create_other_collection(raw_collection_lidvid, bundle_lidvid)
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)
        fits_product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::1.8"
        product_lidvid = "urn:nasa:pds:hst_99999:browse_acs_raw:p::1.8"
        self.db.create_fits_product(fits_product_lidvid, raw_collection_lidvid)
        self.db.create_browse_product(
            product_lidvid, fits_product_lidvid, collection_lidvid
        )
        basename = "p.jpg"
        self.assertFalse(self.db.browse_file_exists(basename, product_lidvid))
        self.db.create_browse_file(
            self.dummy_os_filepath, basename, product_lidvid, 1024
        )
        self.assertTrue(self.db.browse_file_exists(basename, product_lidvid))

    def test_document_file_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:document::1.8"
        self.db.create_document_collection(collection_lidvid, bundle_lidvid)

        product_lidvid = "urn:nasa:pds:hst_99999:document:p::8.1"
        self.db.create_document_product(product_lidvid, collection_lidvid)

        filename = "stuff.pdf"
        self.assertFalse(self.db.document_file_exists(filename, product_lidvid))

        self.db.create_document_file(self.dummy_os_filepath, filename, product_lidvid)
        self.assertTrue(self.db.document_file_exists(filename, product_lidvid))

    def test_fits_file_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::8.1"
        self.db.create_fits_product(product_lidvid, collection_lidvid)

        filename = "p_raw.fits"
        self.assertFalse(self.db.fits_file_exists(filename, product_lidvid))

        self.db.create_fits_file(self.dummy_os_filepath, filename, product_lidvid, 1)
        self.assertTrue(self.db.fits_file_exists(filename, product_lidvid))

    def test_get_file(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.8"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:p::8.1"
        self.db.create_fits_product(product_lidvid, collection_lidvid)

        filename = "p_raw.fits"
        with self.assertRaises(Exception):
            self.db.get_file(filename, product_lidvid)

        self.db.create_fits_file(self.dummy_os_filepath, filename, product_lidvid, 1)
        self.assertTrue(self.db.get_file(filename, product_lidvid))

    ############################################################

    def test_create_bundle_label(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        self.assertFalse(self.db.bundle_label_exists(bundle_lidvid))
        basename = "bundle.xml"
        self.db.create_bundle_label(self.dummy_os_filepath, basename, bundle_lidvid)
        self.assertTrue(self.db.bundle_label_exists(bundle_lidvid))
        self.db.create_bundle_label(self.dummy_os_filepath, basename, bundle_lidvid)
        self.assertTrue(self.db.bundle_label_exists(bundle_lidvid))

    def test_bundle_label_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        self.assertFalse(self.db.bundle_label_exists(bundle_lidvid))
        basename = "bundle.xml"
        self.db.create_bundle_label(self.dummy_os_filepath, basename, bundle_lidvid)
        self.assertTrue(self.db.bundle_label_exists(bundle_lidvid))

    def test_get_bundle_label(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)

        basename = "bundle.xml"
        self.db.create_bundle_label(self.dummy_os_filepath, basename, bundle_lidvid)
        bundle_label = self.db.get_bundle_label(bundle_lidvid)
        self.assertEqual(bundle_lidvid, bundle_label.bundle_lidvid)
        self.assertEqual(basename, bundle_label.basename)
        self.assertEqual(file_md5(self.dummy_os_filepath), bundle_label.md5_hash)

    def test_create_collection_label(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        self.assertFalse(self.db.collection_label_exists(collection_lidvid))
        basename = "collection.xml"
        self.db.create_collection_label(
            self.dummy_os_filepath, basename, collection_lidvid
        )
        self.assertTrue(self.db.collection_label_exists(collection_lidvid))
        self.db.create_collection_label(
            self.dummy_os_filepath, basename, collection_lidvid
        )
        self.assertTrue(self.db.collection_label_exists(collection_lidvid))

    def test_collection_label_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        self.assertFalse(self.db.collection_label_exists(collection_lidvid))
        basename = "collection.xml"
        self.db.create_collection_label(
            self.dummy_os_filepath, basename, collection_lidvid
        )
        self.assertTrue(self.db.collection_label_exists(collection_lidvid))

        # TODO Fix this.

    #     def test_get_collection_label(self) -> None:
    #         bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
    #         self.db.create_bundle(bundle_lidvid)
    #         collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
    #         self.db.create_other_collection(collection_lidvid, bundle_lidvid)
    #
    #         basename = get_collection_label_name(self.db, collection_lidvid)
    #         self.db.create_collection_label(
    #             self.dummy_os_filepath, basename, collection_lidvid
    #         )
    #         collection_label = self.db.get_collection_label(collection_lidvid)
    #         self.assertEqual(collection_lidvid, collection_label.collection_lidvid)
    #         self.assertEqual(basename, collection_label.basename)
    #         self.assertEqual(file_md5(self.dummy_os_filepath), collection_label.md5_hash)

    def test_create_collection_inventory(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        self.assertFalse(self.db.collection_inventory_exists(collection_lidvid))
        basename = "collection.xml"
        self.db.create_collection_inventory(
            self.dummy_os_filepath, basename, collection_lidvid
        )
        self.assertTrue(self.db.collection_inventory_exists(collection_lidvid))
        self.db.create_collection_inventory(
            self.dummy_os_filepath, basename, collection_lidvid
        )
        self.assertTrue(self.db.collection_inventory_exists(collection_lidvid))

    def test_collection_inventory_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)

        self.assertFalse(self.db.collection_inventory_exists(collection_lidvid))
        basename = "collection.xml"
        self.db.create_collection_inventory(
            self.dummy_os_filepath, basename, collection_lidvid
        )
        self.assertTrue(self.db.collection_inventory_exists(collection_lidvid))

        # TODO Fix this.

    #     def test_get_collection_inventory(self) -> None:
    #         bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
    #         self.db.create_bundle(bundle_lidvid)
    #         collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
    #         self.db.create_other_collection(collection_lidvid, bundle_lidvid)
    #
    #         basename = get_collection_inventory_name(self.db, collection_lidvid)
    #         self.db.create_collection_inventory(
    #             self.dummy_os_filepath, basename, collection_lidvid
    #         )
    #         collection_inventory = self.db.get_collection_inventory(collection_lidvid)
    #         self.assertEqual(collection_lidvid, collection_inventory.collection_lidvid)
    #         self.assertEqual(basename, collection_inventory.basename)
    #         self.assertEqual(
    #             file_md5(self.dummy_os_filepath), collection_inventory.md5_hash
    #         )

    def test_create_product_label(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)
        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:j6gp01lzq::1.8"
        self.db.create_fits_product(product_lidvid, collection_lidvid)

        self.assertFalse(self.db.product_label_exists(product_lidvid))
        basename = f"{LIDVID(product_lidvid).lid().product_id}.xml"
        self.db.create_product_label(self.dummy_os_filepath, basename, product_lidvid)
        self.assertTrue(self.db.product_label_exists(product_lidvid))
        self.db.create_product_label(self.dummy_os_filepath, basename, product_lidvid)
        self.assertTrue(self.db.product_label_exists(product_lidvid))

    def test_product_label_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)
        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:j6gp01lzq::1.8"
        self.db.create_fits_product(product_lidvid, collection_lidvid)

        self.assertFalse(self.db.product_label_exists(product_lidvid))
        basename = f"{LIDVID(product_lidvid).lid().product_id}.xml"
        self.db.create_product_label(self.dummy_os_filepath, basename, product_lidvid)
        self.assertTrue(self.db.product_label_exists(product_lidvid))

    def test_get_product_label(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_bundle(bundle_lidvid)
        collection_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw::1.1"
        self.db.create_other_collection(collection_lidvid, bundle_lidvid)
        product_lidvid = "urn:nasa:pds:hst_99999:data_acs_raw:j6gp01lzq::1.8"
        self.db.create_fits_product(product_lidvid, collection_lidvid)

        basename = "j6gp01lzq_raw.fits"
        self.db.create_product_label(self.dummy_os_filepath, basename, product_lidvid)
        product_label = self.db.get_product_label(product_lidvid)
        self.assertEqual(product_lidvid, product_label.product_lidvid)
        self.assertEqual(basename, product_label.basename)
        self.assertEqual(file_md5(self.dummy_os_filepath), product_label.md5_hash)

    ############################################################

    def test_get_card_dictionaries(self) -> None:
        # We don't test the functionality of get_card_dictionaries()
        # here, since we'd first need to populate the database with
        # the contents of the FITS file and that machinery is in
        # fits_file_db.  So the real testing happens there.  I leave
        # this empty test as the equivalent of a "this page left
        # intentionally blank" message.
        pass

    def test_exploratory(self) -> None:
        # what happens if you create tables twice?
        self.db.create_tables()
        # no exception, at least

    ############################################################

    def test_get_proposal_info(self) -> None:
        bundle_lid = "urn:nasa:pds:hst_99999"
        self.assertFalse(self.db.proposal_info_exists(bundle_lid))
        with self.assertRaises(Exception):
            self.db.get_proposal_info(bundle_lid)

        self.db.create_proposal_info(
            bundle_lid,
            proposal_title="Proposal of marriage",
            pi_name="Romeo Montague",
            author_list="R. Montague, B. Montague",
            proposal_year="1476",
            publication_year="1597",
        )
        self.assertTrue(self.db.proposal_info_exists(bundle_lid))

        not_a_bundle_lid = "urn:nasa:pds:space_1999"
        self.assertFalse(self.db.proposal_info_exists(not_a_bundle_lid))

        proposal_info = self.db.get_proposal_info(bundle_lid)
        self.assertEqual("Proposal of marriage", proposal_info.proposal_title)

    ############################################################
    def test_create_citation(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.assertTrue(
            self.db.session.query(CitationInfo)
            .filter(CitationInfo.lidvid == bundle_lidvid)
            .count()
            == 0
        )
        self.assertFalse(self.db.citation_exists(bundle_lidvid))
        self.db.create_citation(bundle_lidvid)
        self.assertTrue(
            self.db.session.query(CitationInfo)
            .filter(CitationInfo.lidvid == bundle_lidvid)
            .count()
            == 1
        )
        self.assertTrue(self.db.citation_exists(bundle_lidvid))
        self.db.create_citation(bundle_lidvid)
        self.assertTrue(self.db.citation_exists(bundle_lidvid))

    def test_citation_exists(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.assertFalse(self.db.citation_exists(bundle_lidvid))
        self.db.create_citation(bundle_lidvid)
        self.assertTrue(self.db.citation_exists(bundle_lidvid))

    def test_get_citation(self) -> None:
        bundle_lidvid = "urn:nasa:pds:hst_99999::1.1"
        self.db.create_citation(bundle_lidvid)
        citation = self.db.get_citation(bundle_lidvid)
        self.assertEqual(bundle_lidvid, citation.lidvid)
        self.assertEqual("{filename}", citation.filename)
        self.assertEqual(99999, citation.propno)
        self.assertEqual("{category}", citation.category)
        self.assertEqual(1, citation.cycle)
        self.assertEqual("{author_1},{author_2}", citation.authors)
        self.assertEqual("{title}", citation.title)
        self.assertEqual(2001, citation.submission_year)
        self.assertEqual(2001, citation.timing_year)
