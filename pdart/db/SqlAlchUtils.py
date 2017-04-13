import os.path

from sqlalchemy import create_engine
from sqlalchemy.orm import Session, sessionmaker

from pdart.db.SqlAlchDBName import DATABASE_NAME

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    import sqlite3

    import pdart.pds4.Bundle as B


def bundle_database_filepath(bundle):
    # type: (B.Bundle) -> unicode
    """
    Given a PDS4 :class:`~pdart.pds4.Bundle.Bundle`, create the
    filepath for the bundle's database.
    """
    return os.path.join(bundle.absolute_filepath(), DATABASE_NAME)


def create_bundle_database_session(bundle):
    # type: (B.Bundle) -> Session
    db_filepath = bundle_database_filepath(bundle)
    return create_database_session(db_filepath)


def create_database_session(db_filepath):
    # type: (unicode) -> Session
    engine = create_engine('sqlite:///' + db_filepath)
    return sessionmaker(bind=engine)()
