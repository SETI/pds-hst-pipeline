"""Representation of a PDS4 VID."""
import re


class VID(object):
    """Representation of a PDS4 VID."""

    def __init__(self, str):
        # type: (unicode) -> None
        """
        Create a VID object from a string, raising an exception if
        the VID string is malformed.
        """
        vs = str.split('.')

        # Check requirements
        assert len(str) <= 255, 'VID is too long'
        assert len(vs) in [1, 2], 'VID has too many components'
        for v in vs:
            assert re.match('\\A(0|[1-9][0-9]*)\\Z', v), \
                'VID is non-numeric: %s' % v

        self._VID = str
        self._major = int(vs[0])

        if len(vs) == 2:
            self._minor = int(vs[1])
        else:
            self._minor = None

    def major(self):
        """Return the major version number."""
        # type: () -> int
        return self._major

    def minor(self):
        """Return the minor version number."""
        # type: () -> int
        return self._minor or 0

    def next_major_vid(self):
        """Return the next major VID."""
        # type: () -> VID
        return VID('%d' % (self.major() + 1))

    def next_minor_vid(self):
        """Return the next minor VID."""
        # type: () -> VID
        return VID('%d.%d' % (self.major(), self.minor() + 1))

    def __cmp__(self, other):
        res = self.major() - other.major()
        if res == 0:
            res = self.minor() - other.minor()
        return res

    def __str__(self):
        return self._VID

    def __repr__(self):
        return 'VID(%r)' % self._VID
