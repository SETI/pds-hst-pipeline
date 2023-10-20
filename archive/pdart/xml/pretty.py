"""
Pretty-printing functionality.
"""
from pdart.xml.schema import run_subprocess, verify_label_or_raise


def pretty_print(str: bytes) -> bytes:
    """Reformat XML using xmllint --format."""
    (exit_code, stderr, stdout) = run_subprocess(["xmllint", "--format", "-"], str)
    if exit_code == 0:
        return stdout
    else:
        # ignore stdout
        raise Exception("pretty_print failed")


def pretty_and_verify(label: bytes, verify: bool) -> bytes:
    if label[:6] != b"<?xml ":
        raise ValueError("label is Not XML.")
    label = pretty_print(label)

    if label[:6] != b"<?xml ":
        raise ValueError("label is Not XML.")
    if verify:
        verify_label_or_raise(label)
    return label
