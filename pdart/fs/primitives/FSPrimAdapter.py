import os
import stat
from typing import Any, BinaryIO, Dict, List, Mapping, Optional, Tuple, cast

import fs.mode
import fs.path
import fs.subfs
from fs.base import FS
from fs.enums import ResourceType
from fs.error_tools import convert_os_errors
from fs.info import Info
from fs.permissions import Permissions

from pdart.fs.primitives.FSPrimitives import Dir, FSPrimitives, Node

_WINDOWS_PLATFORM: bool = False


def _remove_exclusive_flag(m: fs.mode.Mode) -> fs.mode.Mode:
    return fs.mode.Mode(m._mode.replace("x", "w"))


class FSPrimAdapter(FS):
    def __init__(self, fs_prims: FSPrimitives) -> None:
        FS.__init__(self)
        self.prims = fs_prims

        _meta = self._meta = {
            "case_insensitive": os.path.normcase("Aa") != "aa",
            "network": False,
            "read_only": False,
            "supports_rename": False,
            "thread_safe": True,
            "unicode_paths": False,
            "virtual": False,
            "invalid_path_chars": "\0",
        }

    def getinfo(self, path: str, namespaces: Any = None) -> Info:
        # The pyfilesystem2 documentation says namespaces should be a
        # list of strings, but the test-suite has a case expecting it
        # to succeed when it's a single string.  Geez.

        # I REALLY REALLY hate untyped languages.
        self.check()
        if not namespaces:
            namespaces = ["basic"]
        if type(namespaces) is not list:
            namespaces = [namespaces]
        node = self._resolve_path_to_node(path)
        if not node:
            raise ValueError(f"Invalid path: {path}")
        info = {}  # Dict[str, Dict[str, object]]
        info["basic"] = {
            "is_dir": self.prims.is_dir(node),
            "name": fs.path.basename(node.path),
        }
        if "details" in namespaces:
            sys_path = self.getsyspath(path)
            if sys_path:
                with convert_os_errors("getinfo", path):
                    _stat = os.stat(sys_path)
                info["details"] = self._make_details_from_stat(_stat)
            else:
                info["details"] = self._make_default_details(node)

        return Info(info)

    def listdir(self, path: str) -> List[str]:
        self.check()
        prims = self.prims
        node = self._resolve_path_to_node(path)
        if prims.is_file(node):
            raise fs.errors.DirectoryExpected(path)
        else:
            return list(prims.get_dir_children(cast(Dir, node)))

    def makedir(
        self,
        path: str,
        permissions: Optional[Permissions] = None,
        recreate: bool = False,
    ) -> fs.subfs.SubFS:
        self.check()
        parts = fs.path.iteratepath(fs.path.abspath(path))
        if not parts:  # we're looking at the root
            if recreate:
                return fs.subfs.SubFS(self, self.prims.root_node().path)
            else:
                raise fs.errors.DirectoryExists(path)
        else:
            prims = self.prims
            parent_dir_node, name = self._resolve_path_to_parent_and_name(path)
            try:
                child = prims.get_dir_child(parent_dir_node, name)
            except KeyError:
                # it doesn't exist
                child = prims.add_child_dir(parent_dir_node, name)
                return fs.subfs.SubFS(self, child.path)
            # it exists
            if prims.is_file(child):
                # TODO This is wrong, but the pyfilesystem test suite
                # asks for it...
                raise fs.errors.DirectoryExists(path)
            else:
                if recreate:
                    return fs.subfs.SubFS(self, child.path)
                else:
                    raise fs.errors.DirectoryExists(path)

    def openbin(
        self, path: str, mode: str = "r", buffering: int = -1, **options: Any
    ) -> BinaryIO:
        self.check()
        self.validatepath(path)
        if path == "/":
            # TODO  Hackish special case.  Clean this up.
            raise fs.errors.FileExpected(path)
        m = fs.mode.Mode(mode)
        m.validate_bin()
        prims = self.prims
        parent_dir_node, name = self._resolve_path_to_parent_and_name(path)

        if "t" in mode:
            raise ValueError(f"openbin() called with text mode {mode}")
        exists = prims.is_dir(parent_dir_node) and name in prims.get_children(
            parent_dir_node
        )
        if exists:
            if m.exclusive:
                raise fs.errors.FileExists(path)
            else:
                file = prims.get_dir_child(parent_dir_node, name)
        elif m.create:
            file = prims.add_child_file(parent_dir_node, name)
            # remove exclusive, since the file now exists
            m = _remove_exclusive_flag(m)
        else:
            raise fs.errors.ResourceNotFound(path)
        return cast(BinaryIO, prims.get_handle(file, m.to_platform()))

    def remove(self, path: str) -> None:
        self.check()
        prims = self.prims
        parent_dir_node, name = self._resolve_path_to_parent_and_name(path)
        try:
            dir = prims.get_dir_child(parent_dir_node, name)
        except KeyError:
            raise fs.errors.ResourceNotFound(path)
        if prims.is_file(dir):
            prims.remove_child(parent_dir_node, name)
        else:
            raise fs.errors.FileExpected(path)

    def removedir(self, path: str) -> None:
        self.check()
        prims = self.prims
        try:
            parent_dir_node, name = self._resolve_path_to_parent_and_name(path)
        except IndexError:
            raise fs.errors.RemoveRootError(path)
        try:
            dir = prims.get_dir_child(parent_dir_node, name)
        except KeyError:
            raise fs.errors.ResourceNotFound(path)
        if prims.is_dir(dir):
            if prims.get_dir_children(cast(Dir, dir)):
                raise fs.errors.DirectoryNotEmpty(path)
            else:
                prims.remove_child(parent_dir_node, name)
        else:
            raise fs.errors.DirectoryExpected(path)

    def setinfo(self, path: str, info: Mapping[str, Mapping[str, object]]) -> None:
        self.check()
        # Check for errors.
        self._resolve_path_to_node(path)
        if "details" in info:
            sys_path = self.getsyspath(path)
            if sys_path:
                details = info["details"]
                if "accessed" in details or "modified" in details:
                    accessed = cast(int, details.get("accessed"))
                    modified = cast(int, details.get("modified", accessed))
                    with convert_os_errors("setinfo", path):
                        os.utime(sys_path, (accessed, modified))

        return None

    def _resolve_path_to_node(self, path: str) -> Node:
        prims = self.prims
        node: Node = prims.root_node()
        try:
            for nm in fs.path.iteratepath(path):
                node = prims.get_dir_child(cast(Dir, node), nm)
            return node
        except KeyError:
            raise fs.errors.ResourceNotFound(path)

    def _resolve_path_to_parent_and_name(self, path: str) -> Tuple[Dir, str]:
        prims = self.prims
        node: Node = prims.root_node()
        parts = fs.path.iteratepath(path)
        try:
            for nm in parts[:-1]:
                node = prims.get_dir_child(cast(Dir, node), nm)
        except KeyError as e:
            print("e =", e)
            print("nm =", repr(nm))
            print(
                (
                    f"prims.get_dir_children({node!r}) = {prims.get_dir_children(cast(Dir, node))!r}"
                )
            )
            raise fs.errors.ResourceNotFound(path)
        return (cast(Dir, node), parts[-1])

    @classmethod
    def _make_details_from_stat(cls, stat_result: Any) -> Dict[str, Any]:
        """Make a *details* info dict from an `os.stat_result` object."""
        details: Dict[str, Any] = {
            "_write": ["accessed", "modified"],
            "accessed": stat_result.st_atime,
            "modified": stat_result.st_mtime,
            "size": stat_result.st_size,
            "type": int(cls._get_type_from_stat(stat_result)),
        }

        # On other Unix systems (such as FreeBSD), the following
        # attributes may be available (but may be only filled out if
        # root tries to use them):
        details["created"] = getattr(stat_result, "st_birthtime", None)
        ctime_key = "created" if _WINDOWS_PLATFORM else "metadata_changed"
        details[ctime_key] = stat_result.st_ctime
        return details

    def _make_default_details(self, node: Node) -> Dict[str, Any]:
        """Make a default *details* info dict"""
        prims = self.prims
        if prims.is_dir(node):
            resource_type = ResourceType.directory
        elif prims.is_file(node):
            resource_type = ResourceType.file
        else:
            resource_type = ResourceType.unknown

        details: Dict[str, Any] = {
            "accessed": None,
            "created": None,
            "modified": None,
            "size": 0,
            "type": resource_type,
        }
        return details

    @classmethod
    def _get_type_from_stat(cls, _stat: Any) -> ResourceType:
        """Get the resource type from an `os.stat_result` object."""
        st_mode = _stat.st_mode
        st_type = stat.S_IFMT(st_mode)
        return cls.STAT_TO_RESOURCE_TYPE.get(st_type, ResourceType.unknown)

    STAT_TO_RESOURCE_TYPE = {
        stat.S_IFDIR: ResourceType.directory,
        stat.S_IFCHR: ResourceType.character,
        stat.S_IFBLK: ResourceType.block_special_file,
        stat.S_IFREG: ResourceType.file,
        stat.S_IFIFO: ResourceType.fifo,
        stat.S_IFLNK: ResourceType.symlink,
        stat.S_IFSOCK: ResourceType.socket,
    }
