"""Representation of a PDS4 LID."""
import functools
import re

from typing import List


@functools.total_ordering
class LID(object):
    """Representation of a PDS4 LID."""

    def __init__(self, lid_str: str) -> None:
        """
        Create a LID object from a string, raising an exception if
        the LID string is malformed.
        """

        ids = lid_str.split(":")

        # Check requirements
        if len(lid_str) > 255:
            raise ValueError(f"{lid_str} length should be <= 255.")
        if len(ids) not in [4, 5, 6]:
            raise ValueError(f"The number of {lid_str} sections are out of range.")
        if ids[0] != "urn" or ids[1] != "nasa" or ids[2] != "pds":
            raise ValueError(f"{lid_str} doesn't start with 'urn:nasa:pds'.")
        allowed_chars_re = r"\A[-._a-z0-9]+\Z"
        for id_ in ids:
            if not re.match(allowed_chars_re, id_):
                raise ValueError(f"{id_} doesn't match the expected pattern.")

        # Assign the Id fields
        self.lid = lid_str
        self.bundle_id = ids[3]

        # ...so this indexing of ids is safe
        self.collection_id = ids[4] if len(ids) > 4 else None
        self.product_id = ids[5] if len(ids) > 5 else None

    @staticmethod
    def create_from_parts(parts: List[str]) -> "LID":
        parts_len = len(parts)
        if parts_len not in [1, 2, 3]:
            raise ValueError(f"{parts} length: {parts_len} is out of range.")
        if parts_len == 1:
            return LID(f"urn:nasa:pds:{parts[0]}")
        elif parts_len == 2:
            return LID(f"urn:nasa:pds:{parts[0]}:{parts[1]}")
        else:
            return LID(f"urn:nasa:pds:{parts[0]}:{parts[1]}:{parts[2]}")

    def parts(self) -> List[str]:
        ids = [self.bundle_id, self.collection_id, self.product_id]
        return [p for p in ids if p is not None]

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, LID):
            return NotImplemented
        return str(self) < str(other)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, LID):
            return NotImplemented
        return str(self) == str(other)

    def __hash__(self) -> int:
        return hash(self.lid)

    def __str__(self) -> str:
        return self.lid

    def __repr__(self) -> str:
        return f"LID({self.lid!r})"

    def is_product_lid(self) -> bool:
        """Return True iff the LID is a product LID."""
        return self.product_id is not None

    def is_collection_lid(self) -> bool:
        """Return True iff the LID is a collection LID."""
        return self.collection_id is not None and self.product_id is None

    def is_bundle_lid(self) -> bool:
        """Return True iff the LID is a bundle LID."""
        return self.collection_id is None

    def parent_lid(self) -> "LID":
        """
        Return a LID object for the object's parent.  Throw an error
        iff the object is a bundle LID.
        """
        if self.is_bundle_lid():
            raise ValueError(f"bundle LID {self.lid!r} has no parent LID")
        else:
            parts = self.lid.split(":")
            return LID(":".join(parts[:-1]))

    def to_browse_lid(self) -> "LID":
        """
        Convert a LID within a data collection into the corresponding
        LID in the browse collection.
        """
        if not self.collection_id:
            raise ValueError("to_browse_lid() -> None: Can't call on bundle LID.")
        collection_id_parts = self.collection_id.split("_")
        if collection_id_parts[0] != "data":
            raise ValueError(
                f"to_browse_lid: Only legal within data_ collections; had {self}."
            )
        collection_id_parts[0] = "browse"
        browse_collection_id = "_".join(collection_id_parts)

        lid_parts = self.lid.split(":")
        lid_parts[4] = browse_collection_id
        browse_collection_lid = ":".join(lid_parts)
        return LID(browse_collection_lid)

    def to_data_lid(self) -> "LID":
        """
        Convert a LID within a browse collection into the
        corresponding LID in the data collection.
        """
        if not self.collection_id:
            raise ValueError("to_data_lid() -> None: Can't call on bundle LID.")
        collection_id_parts = self.collection_id.split("_")
        if collection_id_parts[0] != "browse":
            raise ValueError(
                f"to_data_lid: Only legal within browse_ collections; had {self}."
            )
        collection_id_parts[0] = "data"
        data_collection_id = "_".join(collection_id_parts)

        lid_parts = self.lid.split(":")
        lid_parts[4] = data_collection_id
        data_collection_lid = ":".join(lid_parts)
        return LID(data_collection_lid)

    def to_other_suffixed_lid(self, suffix: str) -> "LID":
        """
        Convert a product LID into the corresponding LID for a file
        with a different suffix.
        """
        if not self.collection_id:
            raise ValueError("to_other_suffixed_lid(): Can't call on bundle LID.")
        collection_id_parts = self.collection_id.split("_")
        if collection_id_parts[0] != "data":
            raise ValueError(
                f"to_other_suffixed_lid: Only legal within "
                + f"data_ collections; had {self}."
            )
        # replace the suffix
        collection_id_parts[2] = suffix
        other_collection_id = "_".join(collection_id_parts)

        lid_parts = self.lid.split(":")
        lid_parts[4] = other_collection_id
        other_product_lid = ":".join(lid_parts)
        return LID(other_product_lid)

    def to_raw_lid(self) -> "LID":
        """
        Convert a product LID into the corresponding LID for a RAW file.
        """
        return self.to_other_suffixed_lid("raw")

    def to_shm_lid(self) -> "LID":
        """
        Convert a product LID into the corresponding LID for a SHM file.
        """
        return self.to_other_suffixed_lid("shm")

    def extend_lid(self, segment: str) -> "LID":
        """
        Create a new LID by extending this one with  another segment.
        """
        ps = self.parts()
        ps.append(segment)
        return LID.create_from_parts(ps)
