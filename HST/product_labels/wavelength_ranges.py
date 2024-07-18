################################################################################
# wavelength_ranges() returns a list of wavelength ranges to use within the
# Science_Facets of HST product labels.
################################################################################

import re

# This is so we can use abbreviations below
FULL_NAMES = {
    "UV": "Ultraviolet",
    "VIS": "Visible",
    "NIR": "Near Infrared",
}

def wavelength_ranges(instrument_id, detector_ids, filter_name):
    """A list of wavelength ranges, based on HST parameters. Each item in the
    list is one of these strings:
        "'Ultraviolet", "Visible", "Near Infrared", "Infrared"

    Inputs are the values obtained from functions defined in HstParameters.py:
        instrument_id   the string returned by get_instrument_id()
        detector_ids    the list of strings returned by get_detector_ids()
        filter_name     the string returned by get_filter_name()
    """

    # Get the list as abbreviations (a subset of "UV", "VIS", "NIR")
    abbrevs = wavelength_abbrevs(instrument_id, detector_ids, filter_name)

    # Expand the abbreviations
    ranges = [FULL_NAMES[abbrev] for abbrev in abbrevs]

    # Append 'Infrared' if necessary
    if "NIR" in abbrevs:
        return ranges + ["Infrared"]
    else:
        return ranges

# This dictionary contains values that override the filter-based algorithm
# The key is one of:
#   instrument_id
#   (instrument_id, detector_id)
#   (instrument_id, detector_id, filter_name)

EASY_TRANSLATIONS = {
    "NICMOS": ["NIR"],
    "COS"   : ["UV"],
    "GHRS"  : ["UV"],
    "FGS"   : ["VIS"],
    ("ACS", "SBN"): ["UV"],
    ("HSP", "UV1"): ["UV"],
    ("HSP", "UV2"): ["UV"],
    ("HSP", "VIS"): ["UV", "VIS", "NIR"],
    ("HSP", "POL"): ["UV", "VIS", "NIR"],
    ("HSP", "PMT"): ["UV", "VIS", "NIR"],
    ("STIS", "FUV-MAMA"): ["UV"],
    ("STIS", "NUV-MAMA"): ["UV"],
    ("STIS", "CCD", "G230L" ): ["UV"],
    ("STIS", "CCD", "G230LB"): ["UV"],
    ("STIS", "CCD", "G230MB"): ["UV"],
    ("STIS", "CCD", "G430L" ): ["UV" , "VIS"],
    ("STIS", "CCD", "G750L" ): ["VIS", "NIR"],
    ("STIS", "CCD", "G750M" ): ["VIS", "NIR"],
    ("STIS", "CCD", "MIRVIS"): ["VIS", "NIR"],
    ("WFC3", "IR"): ["NIR"],
    ("FOS", "AMBER", "MIRROR"): ["UV", "VIS", "NIR"],
    ("FOS", "AMBER", "PRISM" ): ["UV", "VIS", "NIR"],
    ("FOS", "BLUE" , "MIRROR"): ["UV", "VIS"],
    ("FOS", "BLUE" , "PRISM" ): ["UV", "VIS"],
}

def wavelength_abbrevs(instrument_id, detector_ids, filter_name):
    """A list of abbreviated wavelength ranges: "UV", "VIS", and/or "NIR"."""

    # Handle easy cases, where wavelength range only depends on instrument or
    # instrument/detector. This also handles the special filter "PRISM" on FOS,
    # whose range depends on the detector.
    if instrument_id in EASY_TRANSLATIONS:
        return EASY_TRANSLATIONS[instrument_id]

    results = []
    for det in detector_ids:
        key = (instrument_id, det)
        if key in EASY_TRANSLATIONS:
            results += EASY_TRANSLATIONS[key]

        key = (instrument_id, det, filter_name)
        if key in EASY_TRANSLATIONS:
            results += EASY_TRANSLATIONS[key]

        # Strip a STIS suffix, e.g., "+ND_3" or "+LONG_PASS", and try again
        if instrument_id == "STIS" and "+" in filter_name:
            key = (instrument_id, det, filter_name.partition('+')[0])
            if key in EASY_TRANSLATIONS:
                results += EASY_TRANSLATIONS[key]

    if results:
        return ranges_union(results)

    # Return ranges based on filter name(s)
    return ranges_from_filter(filter_name)

# Boundaries between the definitions of UV, VIS, and NIR, with tiny adjustments
# for known filters
UV_MAX  = 400
VIS_MIN = 388
VIS_MAX = 702
NIR_MIN = 647

# This regular expression matches the number inside a filter name, e.g.,
# 'F606W'. This is generally the center wavelength of the bandpass in nm.
FILTER_REGEX = re.compile(r"[FG][QR]?(\d+)[A-Z].*")

# This is a list of filter/grating names whose ranges would not be correctly
# inferred from the name based on the usual algorithm. This is needed for two
# reasons:
#   - the filter name is non-standard, e.g., "CLEAR" or "POL60V"
#   - the filter is wide enough to encompass a range that is outside the one
#     expected based on the center wavelength alone.

# Instrument names in comments after each line indicate that MST has verified
# the entry from that instrument's handbook in Spring 2022.

FILTER_EXCEPTIONS = {
    "CLEAR" : ["UV", "VIS", "NIR"],
    "F130LP": ["UV", "VIS", "NIR"],     # FOC, WFPC2
    "F165LP": ["UV", "VIS", "NIR"],     # ACS
    "F180LP": ["UV", "VIS"],            # FOC
    "F200LP": ["UV", "VIS", "NIR"],     # WFC3
    "F305LP": ["UV", "VIS"],            # FOC
    "F320LP": ["UV", "VIS"],            # FOC
    "F350LP": ["UV", "VIS", "NIR"],     # WFC3
    "F370LP": ["UV", "VIS"],            # FOC
    "F372M" : ["UV", "VIS"],            # FOC
    "F380W" : ["UV", "VIS"],            # WFPC2
    "F410M" : ["UV", "VIS"],            # WFPC2
    "F430W" : ["UV", "VIS"],            # FOC
    "F435W" : ["UV", "VIS"],            # ACS
    "F600LP": ["VIS", "NIR"],           # WFC3
    "F606W" : ["VIS", "NIR"],           # WFPC2
    "F622W" : ["VIS", "NIR"],           # WFPC2
    "F625W" : ["VIS", "NIR"],           # ACS
    "F675W" : ["VIS", "NIR"],           # WFPC2
    "F702W" : ["VIS", "NIR"],           # WFPC2
    "F718M" : ["VIS", "NIR"],           # WF/PC
    "FQCH4N"   : ["VIS", "NIR"],        # WFPC2
    "FQCH4N15" : ["VIS", "NIR"],        # WFPC2
    "FQCH4N33" : ["VIS"],               # WFPC2
    "FQCH4P15" : ["VIS", "NIR"],        # WFPC2
    "FQUVN"    : ["UV", "VIS"],         # WFPC2
    "FQUVN33"  : ["UV", "VIS"],         # WFPC2
    "FR418N"   : ["UV", "VIS"],         # WFPC2
    "FR418N18" : ["UV", "VIS"],         # WFPC2
    "FR418N33" : ["UV", "VIS"],         # WFPC2
    "FR418P15" : ["UV", "VIS"],         # WFPC2
    "FR459M"  : ["UV", "VIS"],          # ACS
    "FR680N"  : ["VIS", "NIR"],         # WFPC2
    "FR680N18": ["VIS", "NIR"],         # WFPC2
    "FR680N33": ["NIR"],                # WFPC2
    "FR680P15": ["VIS", "NIR"],         # WFPC2
    "G430L" : ["UV", "VIS"],            # STIS
    "G430M" : ["UV", "VIS"],            # STIS
    "G570H" : ["VIS", "NIR"],           # STIS
    "G650L" : ["UV", "VIS", "NIR"],     # STIS
    "G750L" : ["VIS", "NIR"],           # STIS
    "G750M" : ["VIS", "NIR"],           # STIS
    "G780H" : ["VIS", "NIR"],           # STIS
    "G800L" : ["VIS", "NIR"],           # STIS
    "POL0UV": ["UV", "VIS"],            # ACS
    "POL0V" : ["VIS", "NIR"],           # ACS
    "POL120UV": ["UV", "VIS"],          # ACS
    "POL120V" : ["VIS", "NIR"],         # ACS
    "POL60UV" : ["UV", "VIS"],          # ACS
    "POL60V"  : ["VIS", "NIR"],         # ACS
    "PRISM1"  : ["UV", "VIS"],          # FOC
    "PRISM2"  : ["UV", "VIS"],          # FOC
    "PRISM3"  : ["UV", "VIS"],          # FOC
}

def filter_number(filter_name):
    """The number embedded within a filter name; otherwise, zero."""

    match = FILTER_REGEX.fullmatch(filter_name)
    if not match:
        return 0

    return int(match.group(2))

def ranges_from_one_filter(filter_name):
    """The list of wavelength ranges associated with a single filter name."""

    # Check the list of exceptions first
    if filter_name in FILTER_EXCEPTIONS:
        return FILTER_EXCEPTIONS[filter_name]

    # Otherwise, derive the wavelength ranges from the center wavelength
    # embedded within the filter name.

    # Get the embedded wavelength value
    value = filter_number(filter_name)
    if value == 0:
        raise ValueError("Filter name not handled: " + repr(filter_name))

    # Associate the wavelength with one or more ranges
    results = []
    if value <= UV_MAX:
        results.append("UV")

    if VIS_MIN <= value <= VIS_MAX:
        results.append("VIS")

    if NIR_MIN <= value:
        results.append("NIR")

    if not results:
        raise ValueError("Wavelength range not found: " + filter_name)

    return results

def ranges_from_filter(filter_name):
    """The list of wavelength ranges associated with a filter name, or else a
    sequence of names concatenated with "+"."""

    names = filter_name.split("+")
    ranges_list = [ranges_from_one_filter(name) for name in names]

    # If there are multiple overlapping filters, the returned list is the
    # mathematical intersection of the individual ranges
    return ranges_intersection(ranges_list)

def ranges_intersection(ranges_list):
    """The intersection of multiple ranges."""

    # This is easy if there's only one item
    if len(ranges_list) == 1:
        return ranges_list[0]

    # Determine the intersection
    combined = set(ranges_list[0])
    for ranges in ranges_list[1:]:
        combined = combined.intersection(set(ranges))

    return sorted_abbrevs(combined)

def ranges_union(ranges_list):
    """The union of multiple ranges."""

    # This is easy if there's only item
    if len(ranges_list) == 1:
        return [ranges_list[0]]

    # Determine the union
    combined = ranges_list[0]
    for ranges in ranges_list[1:]:
        combined += ranges

    return sorted_abbrevs(combined)  # also removes duplicates

def sorted_abbrevs(ranges):
    """Range abbreviations sorted to increasing wavelength; duplicates removed."""

    sorted = []
    if "UV" in ranges:
        sorted.append("UV")
    if "VIS" in ranges:
        sorted.append("VIS")
    if "NIR" in ranges:
        sorted.append("NIR")

    return sorted
