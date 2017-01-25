"""Functionality to build a bundle label using a SQLite database."""
from contextlib import closing
import sys

from pdart.pds4labels.BundleLabelXml import *
from pdart.pds4labels.CitationInformation import *
from pdart.pds4labels.DatabaseCaches import *
from pdart.pds4labels.DBCalls import get_bundle_collections_db
from pdart.xml.Pretty import *
from pdart.xml.Schema import *

from typing import cast, TYPE_CHECKING
if TYPE_CHECKING:
    from typing import Iterable


def make_db_bundle_label(conn, lid, verify):
    # type: (sqlite3.Connection, unicode, bool) -> unicode
    """
    Create the label text for the bundle having this
    :class:`~pdart.pds4.LID` using the database connection.  If verify
    is True, verify the label against its XML and Schematron schemas.
    Raise an exception if either fails.
    """
    d = lookup_bundle(conn, lid)
    label_fp = d['label_filepath']
    proposal_id = d['proposal_id']

    with closing(conn.cursor()) as cursor:
        reduced_collections = \
            [make_bundle_entry_member({'lid': collection_lid})
             for (collection_lid,)
             in get_bundle_collections_db(cursor, lid)]

    label = make_label({
            'lid': lid,
            'proposal_id': str(proposal_id),
            'Citation_Information': make_placeholder_citation_information(lid),
            'Bundle_Member_Entries': combine_nodes_into_fragment(
                reduced_collections)
            }).toxml()
    label = pretty_print(label)

    with open(label_fp, 'w') as f:
        f.write(label)

    print 'bundle label for', lid
    sys.stdout.flush()

    if verify:
        verify_label_or_raise(label)

    return label
