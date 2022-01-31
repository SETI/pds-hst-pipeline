from logging import Logger
from typing import Dict, ItemsView, Optional, Tuple
from pdart.pds4.lid import LID
from pdart.pds4.lidvid import LIDVID
from pdart.pds4.vid import VID


class ChangesDict(object):
    """
    This class represents a list of PDS4 components in a version and
    whether they have changed since the previous version or not.
    Internally, it's a list of LID, VID, and a bool: True means it's
    changed.

    This class is written to (and can be read from) the filesystem so
    that the information is available to multiple stages in the
    pipeline.
    """

    def __init__(
        self, changes_dict: Optional[Dict[LID, Tuple[VID, bool]]] = None
    ) -> None:
        self.changes_dict = changes_dict if changes_dict else dict()

    def contains(self, lid: LID) -> bool:
        return lid in self.changes_dict

    def set(self, lid: LID, vid: VID, changed: bool) -> None:
        self.changes_dict[lid] = (vid, changed)

    def vid(self, lid: LID) -> VID:
        return self.changes_dict[lid][0]

    def changed(self, lid: LID) -> bool:
        return self.changes_dict[lid][1]

    def items(self) -> ItemsView[LID, Tuple[VID, bool]]:
        return self.changes_dict.items()

    def parent_lidvid(self, lidvid: LIDVID) -> LIDVID:
        lid = lidvid.lid()
        if lid not in self.changes_dict:
            raise KeyError(f"lid={lid} not in changes_dict.")
        parent_lid = lid.parent_lid()
        if parent_lid not in self.changes_dict:
            raise KeyError(f"parent_lid={parent_lid} not in changes_dict.")
        parent_vid = self.vid(parent_lid)
        return LIDVID.create_from_lid_and_vid(parent_lid, parent_vid)

    def dump(self, note: str = "") -> None:
        print(f"**** ChangeDict {note} ****")
        for lid, (vid, changed) in sorted(self.changes_dict.items()):
            print(lid, vid, changed)
        print("****************")

    def log(self, logger: Logger, level: int, note: str = "") -> None:
        logger.log(level, f"**** ChangeDict {note} ****")
        for lid, (vid, changed) in sorted(self.changes_dict.items()):
            logger.log(level, f"{lid} {vid} {changed}")
        logger.log(level, "****************")

    def has_changes(self) -> bool:
        res = False
        for lid, (vid, changed) in self.changes_dict.items():
            res = res or changed
        return res


CHANGES_DICT = ChangesDict
CHANGES_DICT_NAME: str = "changes$dict.txt"


def read_changes_dict(changes_path: str) -> CHANGES_DICT:
    changes_dict = dict()
    with open(changes_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if parts:
                if len(parts) != 3:
                    raise ValueError(
                        f"Length of parts {parts} for changes_dict is not 3."
                    )
                lid, vid, changed = parts
                if changed not in ["False", "True"]:
                    raise ValueError(f"{changed} from parts has unexpected value.")
                changes_dict[LID(lid)] = (VID(vid), changed == "True")
    return ChangesDict(changes_dict)


def write_changes_dict(d: ChangesDict, changes_path: str) -> None:
    with open(changes_path, "w") as changes_file:
        for lid, (vid, changed) in sorted(d.changes_dict.items()):
            print(lid, vid, changed, file=changes_file)
