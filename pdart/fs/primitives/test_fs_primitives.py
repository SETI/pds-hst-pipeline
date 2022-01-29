import abc
import io
import os
import os.path
import shutil
import tempfile
import unittest
from typing import Any, Dict, IO, cast

import fs.path
from fs.test import FSTestCases

from pdart.fs.primitives.fs_prim_adapter import FSPrimAdapter
from pdart.fs.primitives.fs_primitives import Dir, FSPrimitives, File, Node


class FSPrimitivesTestBase(object, metaclass=abc.ABCMeta):
    """
    This is not a test case, but an abstract base class for a test case.
    """

    @abc.abstractmethod
    def setUp(self) -> None:
        pass

    def test_root_node(self) -> None:
        fs = self.get_fs()
        root = fs.root_node()
        # assert that the root exists...
        self.assertTrue(root)
        # ...and is a directory
        self.assertTrue(fs.is_dir(root))
        # ...and has the right path
        self.assertEqual("/", root.path)

    def test_is_dir(self) -> None:
        fs = self.get_fs()
        root = fs.root_node()
        self.assertTrue(fs.is_dir(root))
        file = fs.add_child_file(root, "foo")
        self.assertFalse(fs.is_dir(file))

    def test_is_file(self) -> None:
        fs = self.get_fs()
        root = fs.root_node()
        self.assertFalse(fs.is_file(root))
        file = fs.add_child_file(root, "foo")
        self.assertTrue(fs.is_file(file))

    def test_get_dir_children(self) -> None:
        fs = self.get_fs()
        root = fs.root_node()
        self.assertFalse(fs.get_dir_children(root))
        file_node = fs.add_child_file(root, "file")
        dir_node = fs.add_child_dir(root, "dir")
        expected = {"file": file_node, "dir": dir_node}
        self.assertEqual(expected, fs.get_dir_children(root))

    def test_get_file_handle(self) -> None:
        fs = self.get_fs()
        root = fs.root_node()
        self.assertFalse(fs.get_children(root))
        file_node = fs.add_child_file(root, "file")
        self.assertTrue(fs.get_file_handle(file_node, "w"))

    def test_add_child_dir(self) -> None:
        fs = self.get_fs()
        root = fs.root_node()
        self.assertFalse(fs.get_children(root))
        self._assert_add_child_dir_is_correct(root, "dir")

    def _assert_add_child_dir_is_correct(self, dir: Dir, child_name: str) -> Dir:
        """
        Check that the primitive action of add_child_dir() does the
        right things.
        """
        old_children = self.get_fs().get_children(dir)
        child_dir = self.get_fs().add_child_dir(dir, child_name)
        # the result exists
        self.assertTrue(child_dir)
        # the result is a directory
        self.assertTrue(self.get_fs().is_dir(child_dir))
        # the result is in the parent's directory
        self.assertEqual(child_dir, self.get_fs().get_dir_child(dir, child_name))
        # the result has the right path
        self.assertEqual(child_dir.path, fs.path.join(dir.path, child_name))
        # assert the children changed, and differ only by the new entry
        new_children = self.get_fs().get_children(dir)
        self.assertEqual(new_children[child_name], child_dir)
        del new_children[child_name]
        self.assertEqual(old_children, new_children)
        return child_dir

    def test_add_child_file(self) -> None:
        fs = self.get_fs()
        root = fs.root_node()
        self.assertFalse(fs.get_children(root))
        self._assert_add_child_file_is_correct(root, "file")
        file = fs.get_dir_child(root, "file")
        self.assertTrue(fs.is_file(file))
        self.assertEqual(file, fs.get_dir_child(root, "file"))
        self.assertEqual("/file", file.path)

    def _assert_add_child_file_is_correct(self, dir: Dir, child_name: str) -> File:
        """
        Check that the primitive action of add_child_file() does the
        right things.
        """
        old_children = self.get_fs().get_children(dir)
        child_file = self.get_fs().add_child_file(dir, child_name)
        # the result exists
        self.assertTrue(child_file)
        # the result is a file
        self.assertTrue(self.get_fs().is_file(child_file))
        # the result is in the parent's directory
        self.assertEqual(child_file, self.get_fs().get_dir_child(dir, child_name))
        # the result has the right path
        self.assertEqual(child_file.path, fs.path.join(dir.path, child_name))
        # assert the children changed, and differ only by the new entry
        new_children = self.get_fs().get_children(dir)
        self.assertEqual(new_children[child_name], child_file)
        del new_children[child_name]
        self.assertEqual(old_children, new_children)
        return child_file

    def test_remove_child(self) -> None:
        fs = self.get_fs()
        root = fs.root_node()
        dir = fs.add_child_dir(root, "dir")
        file = fs.add_child_file(root, "file")
        self._assert_remove_child_is_correct(root, "dir")
        self._assert_remove_child_is_correct(root, "file")

    def _assert_remove_child_is_correct(self, dir: Dir, child_name: str) -> None:
        """
        Check that the primitive action of remove_child() does the
        right things.
        """
        fs = self.get_fs()
        old_children = fs.get_children(dir)
        old_child = fs.get_dir_child(dir, child_name)
        fs.remove_child(dir, child_name)
        # assert the children changed, and differ only by the entry
        new_children = fs.get_children(dir)
        self.assertTrue(child_name not in new_children)
        del old_children[child_name]
        self.assertEqual(old_children, new_children)

    # The following are defined as abstract.  Their implementations
    # come from mixing with unittest.TestCase.  I can't inherit from
    # TestCase here because then py.test will try to construct and run
    # this abstract class.

    @abc.abstractmethod
    def assertTrue(self, cond: Any, msg: object = None) -> None:
        pass

    @abc.abstractmethod
    def assertFalse(self, cond: Any, msg: object = None) -> None:
        pass

    @abc.abstractmethod
    def assertEqual(self, lhs: Any, rhs: Any, msg: object = None) -> None:
        pass

    # This is also defined in the real TestCase.

    @abc.abstractmethod
    def get_fs(self) -> FSPrimitives:
        pass


class OSFSPrimitives(FSPrimitives):
    def __init__(self, root: str) -> None:
        FSPrimitives.__init__(self)
        self.root = root

    def __str__(self) -> str:
        return f"OSFSPrimitives({self.root})"

    def __repr__(self) -> str:
        return f"OSFSPrimitives({self.root!r})"

    def _to_sys_path(self, path: str) -> str:
        path_part = path.lstrip("/")
        return os.path.join(self.root, path_part)

    def add_child_dir(self, parent_node: Dir, filename: str) -> Dir:
        path = fs.path.join(parent_node.path, filename)
        sys_path = self._to_sys_path(path)
        os.mkdir(sys_path)
        return Dir(self, path)

    def add_child_file(self, parent_node: Dir, filename: str) -> File:
        path = fs.path.join(parent_node.path, filename)
        sys_path = self._to_sys_path(path)
        with open(sys_path, "w"):
            pass
        return File(self, path)

    def get_dir_children(self, node: Dir) -> Dict[str, Node]:
        dir_path = node.path
        dir_sys_path = self._to_sys_path(dir_path)
        res = dict()
        for filename in os.listdir(dir_sys_path):
            child_path = fs.path.join(dir_path, filename)
            if os.path.isfile(child_path):
                child_node: Node = File(self, child_path)
            else:
                child_node = Dir(self, child_path)
            res[str(filename)] = child_node
        return res

    def get_file_handle(self, node: File, mode: str) -> IO[Any]:
        sys_path = self._to_sys_path(node.path)
        return io.open(sys_path, fs.mode.Mode(mode).to_platform_bin())

    def is_file(self, node: Node) -> bool:
        sys_path = self._to_sys_path(node.path)
        return os.path.isfile(sys_path)

    def remove_child(self, parent_node: Dir, filename: str) -> None:
        child = self.get_dir_child(parent_node, filename)
        sys_path = self._to_sys_path(child.path)
        if self.is_file(child):
            os.remove(sys_path)
        else:
            os.rmdir(sys_path)

    def root_node(self) -> Dir:
        return Dir(self, "/")


_TMP_DIR: str = os.path.abspath("tmp_osfs_prims")


class TestOSFSPrimitives(unittest.TestCase, FSPrimitivesTestBase):
    def setUp(self) -> None:
        self.tmpdir = tempfile.mkdtemp()
        self.fs = OSFSPrimitives(self.tmpdir)

    def get_fs(self) -> FSPrimitives:
        return self.fs

    def tearDown(self) -> None:
        shutil.rmtree(self.tmpdir)


class OSFSPrimAdapter(FSPrimAdapter):
    def __init__(self, root_dir: str) -> None:
        FSPrimAdapter.__init__(self, OSFSPrimitives(root_dir))

    def getsyspath(self, path: str) -> str:
        return cast(OSFSPrimitives, self.prims)._to_sys_path(path)


class TestOSFSPrimAdapter(FSTestCases, unittest.TestCase):
    def make_fs(self) -> FSPrimAdapter:
        self.tmpdir = tempfile.mkdtemp()
        return OSFSPrimAdapter(self.tmpdir)

    def test_t(self) -> None:
        with self.fs.openbin("test_openbin_exclusive", "x") as f:
            f.write(b"bananas")
