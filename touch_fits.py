"""
A utility to change the bits of a FITS file without changing its
meaning.  To do that, we simply add a HISTORY record.  The name comes
from the Unix utility "touch".
"""
import datetime
import sys
import astropy.io.fits


def touch_fits(filepath: str) -> None:
    hdulist = astropy.io.fits.open(filepath, mode="update")
    hdulist[0].header["history"] = f"touch_fits: {datetime.datetime.now().isoformat()}"
    hdulist.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: python3 touch_fits.py <filepath>", file=sys.stderr)
        sys.exit(1)
    filepath = sys.argv[1]
    touch_fits(filepath)
