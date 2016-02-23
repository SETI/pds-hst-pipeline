import os.path


class ArchiveFile(object):
    """A file belonging to an ArchiveComponent."""

    def __init__(self, comp, basename):
        """
        Create an ArchiveFile given the component it belongs to and
        the basename (that is, filepath without directory part) for
        the file.

        Note that this assumes that products don't contain
        subdirectories.  That won't always be true.
        """

        assert comp
        self.component = comp
        assert basename
        self.basename = basename

    def __eq__(self, other):
        return self.component == other.component and \
            self.basename == other.basename

    def __str__(self):
        return '%s in %s' % (self.basename, repr(self.component))

    def __repr__(self):
        return 'ArchiveFile(%s, %s)' % (repr(self.basename),
                                        repr(self.component))

    def fullFilepath(self):
        """Return the full, absolute filepath to the file."""
        return os.path.join(self.component.directoryFilepath(),
                            self.basename)
