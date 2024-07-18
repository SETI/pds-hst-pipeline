##########################################################################################
# comets module
#
# To import:
#   from comets import comet_identifications(key)
#
# To use:
#   target_identification_list = comet_identifications(keys)
# returns a list of Target_Identification tuples. Each tuple contains:
#   (name, alt_designations, type, description, lid)
# where:
#   name                the preferred name
#   alt_designations    a list of strings indicating alternative names.
#   type                "Comet", "Centaur", etc.
#   description         a list of strings, to be separated by newlines inside
#                       the XML Target_Identification object.
#   lid                 the LID of the object, omitting "urn:...:target:"
#
# The input argument can be a single identification or else a list of one or more such
# identifications. If multiple identifications are provided, they must be internally
# consistent, all referring to the same comet. This can be used as an internal consistency
# check. For example,
#    comet_identifications(['1P', 'Halley'])
# will return the info for Halley's comet, but
#    comet_identifications(['2P', 'Halley'])
# will raise a ValueError, indicating that the two identifications are inconsistent.
#
# The reason comet_identifications returns a list is for possible cases of associated
# objects, such as fragments of a single comet, that might all simultaneously be targets
# of a given observation.
#
# For additional internal checks, you can also call the function with explicit values for
# orbital elements a, e and i. If provided, an error will arise if any specified orbital
# element differs from that in the MPC by the specified threshhold of accuracy.
#
# For comets that have been re-designated to be minor planets, comet_identifications
# returns the minor planet information, but with cometary names included in the list of
# alt_designations.
#
# NOTE: This module uses information extracted manually from two sources:
#  https://pds-smallbodies.astro.umd.edu/data_sb/resources/periodic_comets.shtml
#  https://ssd.jpl.nasa.gov/sbdb_query.cgi#x
# Information from these pages is saved in files PDS_COMETS_TXT.py and SSD_COMETS_CSV.py.
# It will occasionally be necessary to update these two files as new comets are
# discovered, renamed, reclassified, or are assigned numeric NAIF IDs.
##########################################################################################

import sys
import math
import re

from .PDS_COMETS_TXT import PDS_COMETS_TXT
from .SSD_COMETS_CSV import SSD_COMETS_CSV
from .ICQ_COMETS_TXT import ICQ_COMETS_TXT

from .. import minor_planets
from .. import lids

NUMBERED_COMET_PATTERN = r'(\d+)[A-Z]'
NUMBERED_COMET_REGEX = re.compile(NUMBERED_COMET_PATTERN)

COMET_DESIGNATION_PATTERN = r'[A-Z]/-?\d{1,4} [A-Z][A-Z]?[0-9]*'
COMET_DESIGNATION_REGEX = re.compile(COMET_DESIGNATION_PATTERN)

OLD_REGEX1 = re.compile(r'\d{4}[a-z]1?')            # E.g., 1987g1, 1869a
OLD_REGEX2 = re.compile(r'-?\d{1,4}(| [IVX]+)')     # E.g., -239, 1948 VII

NUMBER_SUFFIX_PATTERN = r'(.+?)[- ](\d+)'
NUMBER_SUFFIX_REGEX = re.compile(NUMBER_SUFFIX_PATTERN)

# Indicates whether to print detailed information during intialization
VERBOSE = False

##########################################################################################
# INFORMATION REPAIRS
#
# These definitions may need to be modified as new comets are added or comets are
# reclassified.
##########################################################################################

# A list of comets that have been re-disgnated as minor planets
COMETS_NOW_MINOR_PLANETS = [
    # (number, letter, designation, name, index, MP number)
    ( 95, 'P', '1977 UB'   , 'Chiron'           ,  0,   2060),
    (107, 'P', '1949 W1'   , 'Wilson-Harrington',  1,   4015),
    (133, 'P', '1996 N2'   , 'Elst-Pizarro'     ,  0,   7968),
    (174, 'P', '2000 EC98' , 'Echeclus'         ,  0,  60558),
    (176, 'P', '1999 RE70' , 'LINEAR'           , 52, 118401),
    (288, 'P', '2006 VW139', ''                 ,  0, 300163),
    (282, 'P', '2003 BM80' , ''                 ,  0, 323137),
    (362, 'P', '2008 GO98' , ''                 ,  0, 457175),
]

# Repairs known errors in the SSD comet list
def REPAIR_SSD_COMETS(ssd_comets):

    # SSD comet list does not include "2I" designation for Borisov
    borisov = [c for c in ssd_comets if c.designation == 'C/2019 Q4'][0]
    borisov.number = 2
    borisov.letter = 'I'

    # SSD comet list is missing 'Oumuamua, with SSD refers to as a "hyperbolic
    # asteroid
    oumuamua = CometInfo(1, 'I', 'A/2017 U1', '', "'Oumuamua", 0, 3788040)
    ssd_comets.append(oumuamua)

# Repairs known errors in the PDS comet list
def REPAIR_PDS_COMETS(pds_comets):

    return

##########################################################################################
# Class for comet information
##########################################################################################

class CometInfo(object):

    ######## CONSTRUCTORS

    def __init__(self, number, letter, designation, suffix, name,
                 index=0, naif_id=0, alt_designations=[]):

        def clean_str(string):
            """Removes surrounding whitespace and optional quotes."""
            string = string.strip()
            if string.startswith('"') and string.endswith('"'):
                string = string[1:-1].strip()

            # There should never be a space after a dash. This shows up for some
            # negative years. (This removes up to three spaces after a dash!)
            string = string.replace('-  ', '-')
            string = string.replace('- ', '-')

            return string

        def clean_int(value):
            """Converts from string if necessary; could be blank."""
            if isinstance(value, str):
                value = clean_str(value).strip()
                if value:
                    return int(value)
                else:
                    return 0

            return value

        self.number      = clean_int(number)
        self.letter      = clean_str(letter)
        self.designation = clean_str(designation)
        self.suffix      = clean_str(suffix)
        self.name        = clean_str(name)
        self.index       = clean_int(index)
        self.naif_id     = clean_int(naif_id)
        self.alt_designations = [clean_str(d) for d in alt_designations]

        self.fragments    = []      # List of fragments, if any
        self.parent       = None    # Parent if this is a fragment
        self.index_needed = False   # True if name without index is ambiguous
        self.mp_number    = 0       # Minor planet number if any, used for
                                    # re-classified bodies

        # Tests and further cleanup
        assert len(self.letter) == 1, 'Invalid letter: ' + letter

        # Remove a suffix appended to the designation
        if not self.suffix:
            # Dash separates suffix from designation, but not if the year is negative!
            parts = self.designation.rpartition('-')
            if parts[0] and not parts[0].endswith('/'):
                self.suffix = '-' + parts[2]
                self.designation = parts[0]

        # If the designation is a comet number, fill in the comet number
        if self.designation:
            match = NUMBERED_COMET_REGEX.fullmatch(self.designation)
            if match:
                self.number = int(self.designation[:-1])
                assert self.designation.endswith(self.letter), 'Letter mismatch'
                self.designation = ''

            # Otherwise, designation starts with a letter and slash
            else:
                if '/' not in self.designation:
                    self.designation = self.letter + '/' + self.designation

        # Make sure the designation has a proper format
        if self.designation:
            match = COMET_DESIGNATION_REGEX.fullmatch(self.designation)
            assert match, 'Invalid designation: ' + self.designation

        for designation in self.alt_designations:
            match = (COMET_DESIGNATION_REGEX.fullmatch(designation)
                     or OLD_REGEX1.fullmatch(designation)
                     or OLD_REGEX2.fullmatch(designation))
            assert match, 'Invalid designation: ' + designation

        # If the name ends with a number suffix, set the index
        if self.name and not self.index:
            match = NUMBER_SUFFIX_REGEX.fullmatch(self.name)
            if match:
                self.name = match.group(1)
                self.index = int(match.group(2))

    PDS_PATTERN = r"( *\d*)([A-Z])/(\d{4} \w+)(-\w+|) +([A-Za-z '-]+) ?(\d*)\s*"
    PDS_REGEX = re.compile(PDS_PATTERN)

    def from_pds(rec):
        match = CometInfo.PDS_REGEX.fullmatch(rec)
        if not match:
            raise ValueError('Illegal row in PDS-comets.txt: ' + rec.rstrip())

        (number, letter, designation, suffix, name, index) = match.groups()

        # The PDS list has some mystery suffixes. Ignore.
        suffix = ''

        return CometInfo(number, letter, designation, suffix, name, index)

    SSD_FULLNAME_INDEX_PATTERN = r'.* (\d+)(|-[A-Z]\w*)\)?"'
    SSD_FULLNAME_INDEX_REGEX = re.compile(SSD_FULLNAME_INDEX_PATTERN)

    def from_ssd(rec):
        values = rec.rstrip().split(',')
        if len(values) != 5:
            raise ValueError('Illegal column count in SSD-comets.csv: ' + rec.rstrip())

        (naif_id, fullname, designation, name, letter) = values

        # Get the index value from the fullname
        match = CometInfo.SSD_FULLNAME_INDEX_REGEX.fullmatch(fullname)
        if match:
            index = int(match.group(1))
        else:
            index = 0

        return CometInfo(0, letter, designation, '', name, index, naif_id)

    ICQ_PATTERN1 = (r' *(\d*)([A-Z])/? *(-? *\d+ [A-Z]\w+|)(|-[A-Z])' +
                    r' +\((.*?) *\)\s*(.*)')
    ICQ_PATTERN2 = r'= (|\w+) += ?(|-? *\d+ ?[IVX]*)\s*'

    ICQ_REGEX1 = re.compile(ICQ_PATTERN1)
    ICQ_REGEX2 = re.compile(ICQ_PATTERN2)

    def from_icq(rec):
        match = CometInfo.ICQ_REGEX1.fullmatch(rec)
        if not match:
            raise ValueError('Invalid ICQ record: ' + rec)

        (number, letter, designation, suffix, name, more) = match.groups()
        alt_designations = []

        if more.strip():
            match = CometInfo.ICQ_REGEX2.fullmatch(more)
            old_names = match.groups()
            alt_designations = [n.strip() for n in old_names if n]
        else:
            alt_designations = []

        return CometInfo(number, letter, designation, suffix, name, 0, 0,
                         alt_designations)

    ######## UTILITIES

    def __str__(self):
        return ('CometInfo(' + str(self.number)  + '|' +
                               self.letter       + '|' +
                               self.designation  + '|' +
                               self.suffix       + '|' +
                               self.name         + '|' +
                               str(self.index)   + '|' +
                               str(self.naif_id) + '|' +
                               str(self.alt_designations) + ')')

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other, ignore_suffix=False):
        # ignore_suffix=True identifies comets as equal even if they refer to different
        # fragments

        # Fragment suffixes must match unless otherwise specified
        if not ignore_suffix and self.suffix != other.suffix:
            return False

        # Matching numbers confirm equality
        if self.number and other.number:
            return self.number == other.number

        # One matching designation confirms equality
        self_designations  = set([self.designation]  + self.alt_designations)
        other_designations = set([other.designation] + other.alt_designations)
        overlap = self_designations.intersection(other_designations)
        if overlap and overlap != {''}:
            return True

        return False

    def merge(self, other):

        self.number      = self.number    or other.number
        self.letter      = self.letter    or other.letter
        self.suffix      = self.suffix    or other.suffix
        self.name        = self.name      or other.name
        self.index       = self.index     or other.index
        self.naif_id     = self.naif_id   or other.naif_id
        self.mp_number   = self.mp_number or other.mp_number

        # Merge lists of designations, keeping preferred designation in front
        if self.designation:
            new_designations = [self.designation]
        elif other.designation:
            new_designations = [other.designation]
        else:
            new_designations = ['']

        for d in (self.alt_designations +
                  [other.designation] + other.alt_designations):
            if d and d not in new_designations:
                new_designations.append(d)

        self.designation = new_designations[0]
        self.alt_designations = new_designations[1:]

    def copy(self, deep=False):

        clone = CometInfo(self.number, self.letter, self.designation,
                          self.suffix, self.name, self.index, self.naif_id,
                          self.alt_designations)
        clone.index_needed = self.index_needed
        clone.fragments = clone.fragments

        if deep:
            clone.alt_designations = list(clone.alt_designations)
            clone.fragments = list(clone.fragments)

        return clone

    def full_names(self, ignore_suffix=False):

        # Interpret the name
        if self.name.startswith('Great') and self.name.endswith('comet'):
            great_name = self.name + ' of ' + self.designation[2:6]
            name_alone = ''
            name_index = ''
        elif self.index:
            great_name = ''
            name_alone = self.name
            name_index = self.name + ' ' + str(self.index)
        else:
            great_name = ''
            name_alone = self.name
            name_index = ''

        # Interpret the suffix
        if ignore_suffix:
            suffix = ''
        else:
            suffix = self.suffix

        names = []
        if self.number:
          for n in (name_alone, name_index):
            if not n: continue

            names += [str(self.number) + self.letter + '/' + n + suffix]

          names += [str(self.number) + self.letter + suffix]

        if self.designation:
          for n in (name_index, name_alone):
            if not n: continue

            names += [self.designation + suffix + ' (' + n + ')',
                      self.designation + suffix + ' ' + n]

            if suffix:
                names += [self.designation + ' ' + n + suffix]

        for designation in [self.designation] + self.alt_designations:
          if designation:
            names += [designation + suffix]

        if great_name:
            names += [great_name + suffix]
        elif name_index:
            names += [name_index + suffix]
        elif name_alone and not self.index_needed:
            names += [name_alone + suffix]

        if self.naif_id:
            names += ['NAIF ID ' + str(self.naif_id)]

        return names

    def lid(self):

        if self.number:
            string = f'{self.number}{self.letter}_{self.name}'
        elif self.name.startswith('Great'):
            string = self.designation
        else:
            string = f'{self.designation}_{self.name}'

        if self.index:
            string += f'_{self.index}'

        string += self.suffix

        # Cleanup for LID formation rule compliance
        return lids.clean('comet.' + string)

    def target_identifications(self):

        names = [n for n in self.full_names() if not n.isdigit() and
                                                 not n.startswith('-')]
        targets = [(names[0], names[1:], 'Comet', [], self.lid())]

        if self.fragments:
            targets += [f.target_identifications()[0] for f in self.fragments]

        return targets

    def mpc_key(self):

        if self.number:
            return str(self.number) + self.letter
        else:
            return self.designation

########################################
# Load PDS comets
########################################

recs = PDS_COMETS_TXT.split('\n')
recs = [rec for rec in recs if not rec.rstrip() == '']
pds_comets = [CometInfo.from_pds(rec) for rec in recs]

# Treat Chiron as an asteroid (Centaur), in accordance with MPC
pds_comets = [c for c in pds_comets if c.name != 'Chiron']
REPAIR_PDS_COMETS(pds_comets)

########################################
# Load SSD comets
########################################

SSD_COMETS_CSV = SSD_COMETS_CSV.replace('PANSTARRS', 'PanSTARRS')

recs = SSD_COMETS_CSV.split('\n')
recs = [rec for rec in recs if not rec.rstrip() == '']
ssd_comets = [CometInfo.from_ssd(rec) for rec in recs]
REPAIR_SSD_COMETS(ssd_comets)

########################################
# Load ICQ comets
########################################

recs = ICQ_COMETS_TXT.split('\n')
recs = [rec for rec in recs if not rec.rstrip() == '']
icq_comets = [CometInfo.from_icq(rec) for rec in recs]

# Merge duplicates (of which there can be many, one per apparition)
merged_comets = []
for comet in icq_comets:
    try:
        k = merged_comets.index(comet)
    except ValueError:      # not found
        merged_comets.append(comet)
        continue

    merged_comets[k].merge(comet)

icq_comets = merged_comets

########################################
# Load comets now re-classified
########################################

mp_comets = []
for info in COMETS_NOW_MINOR_PLANETS:
    (number, letter, designation, name, index, mp_number) = info
    comet = CometInfo(number, letter, designation, '', name, index)
    comet.mp_number = mp_number
    mp_comets.append(comet)

########################################
# Merge lists
########################################

for comet in pds_comets:
    try:
        k = ssd_comets.index(comet)
    except ValueError:
        # This happens if a PDS comet does not have a NAIF ID.
        if VERBOSE:
            print('# PDS comet not found in SSD list: ' + str(comet))
        continue

    comet.merge(ssd_comets[k])
    del ssd_comets[k]

# Output in 12/20...
# PDS comet not found in SSD list: CometInfo(0|P|P/2019 T5||ATLAS|9|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2014 J1||Catalina|14|0|[])
# PDS comet not found in SSD list: CometInfo(174|P|P/2000 EC98||Echeclus|0|0|[])
# PDS comet not found in SSD list: CometInfo(133|P|P/1996 N2||Elst-Pizarro|1|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2019 Y2||Fuls|2|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2015 D5||Kowalski|10|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2016 Q4||Kowalski|11|0|[])
# PDS comet not found in SSD list: CometInfo(0|C|C/2017 Y3||Leonard|2|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2005 SD||LINEAR|49|0|[])
# PDS comet not found in SSD list: CometInfo(176|P|P/1999 RE70||LINEAR|52|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2004 V5||LINEAR-Hill|1|0|[])
# PDS comet not found in SSD list: CometInfo(0|C|C/2012 C3||PanSTARRS|7|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2015 R1||PanSTARRS|46|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2016 J1||PanSTARRS|59|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2016 S1||PanSTARRS|61|0|[])
# PDS comet not found in SSD list: CometInfo(0|C|C/2018 A4||PanSTARRS|79|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2015 Q1||Scotti|10|0|[])
# PDS comet not found in SSD list: CometInfo(0|D|D/1993 F2||Shoemaker-Levy|9|0|[])
# PDS comet not found in SSD list: CometInfo(107|P|P/1949 W1||Wilson-Harrington|1|0|[])
# PDS comet not found in SSD list: CometInfo(0|P|P/2010 L5||WISE|7|0|[])

COMETS = pds_comets + ssd_comets

for comet in COMETS:
    try:
        k = icq_comets.index(comet)
    except ValueError:
        continue

    comet.merge(icq_comets[k])
    del icq_comets[k]

COMETS += icq_comets

for comet in mp_comets:
    try:
        k = COMETS.index(comet)
    except ValueError:
        # This happens if a minor planet comet is not already in the list
        if VERBOSE:
            print('# Minor planet comet not found in list: ' + str(comet))
        COMETS.append(comet)
        continue

    COMETS[k].merge(comet)

# Minor planet comet not found in list: CometInfo(288|P|P/2006 VW139|||0|0|[])
# Minor planet comet not found in list: CometInfo(282|P|P/2003 BM80|||0|0|[])
# Minor planet comet not found in list: CometInfo(362|P|P/2008 GO98|||0|0|[])

# Drop index values where the name is otherwise unique and the index is unity.
# Set index_needed to True for reused names.

# COMETS_BY_NAME is a dictionary indexed by the name of the discoverer(s). Each entry is
# a list of comets with the same discoverer. Each unique comet in this list is
# represented by a sub-list, with one entry per fragment. For most comets, which do not
# have multiple fragments, the sub-list has unit length.

comet_list_by_name = {}     # mixture of multiple comets and multiple fragments
for comet in COMETS:
    if comet.name in comet_list_by_name:
        comet_list_by_name[comet.name].append(comet)
    else:
        comet_list_by_name[comet.name] = [comet]

COMETS_BY_NAME = {}         # grouped by unique comet
for (name, comet_list) in comet_list_by_name.items():
    unique_comets = []
    for comet in comet_list:
        appended = False
        for test_list in unique_comets:
            if test_list[0].__eq__(comet, ignore_suffix=True):
                test_list.append(comet)
                appended = True
                break

        if not appended:
            unique_comets.append([comet])

    COMETS_BY_NAME[name] = unique_comets

    if len(unique_comets) > 1:
        for comet in comet_list:
            comet.index_needed = True

    elif comet.index == 1:
        comet.index = 0
        # This happens if a PDS comet has an index of one, but the name is actually
        # unique. It occurs because the PDS list unnecessariily assigns an index to every
        # comet, even if its name is unique.
        if VERBOSE:
            print('# Index removed:', comet)

# Create a big dictionary indexed by every possible name
COMET_LOOKUP = {}           # indexed by any name
BEST_NAME_LOOKUP = {}       # indexed by preferred name

for comet in COMETS:
    full_names = comet.full_names()
    if not full_names:
        raise ValueError('No names: ' + str(comet))

    BEST_NAME_LOOKUP[full_names[0]] = comet

    for name in full_names:
        if name in COMET_LOOKUP:
            # This indicates that one of the names for two different CometInfo objects
            # is not unique. The ideal response is to update the algorithm used by
            # full_names() so that it does not return this name. However, beyond that,
            # this is not a big deal; it just means that two different
            # Target_Identification objects will share the same alt_designation.
            print('WARNING: Duplicated name:', name, COMET_LOOKUP[name], comet)
        else:
            COMET_LOOKUP[name] = comet
            COMET_LOOKUP[name.upper()] = comet
            COMET_LOOKUP[name.upper().replace('-',' ')] = comet

    lid = comet.lid()
    parts = lid.partition('comet.')

    for key in (lid, parts[1] + parts[2], parts[2]):
        COMET_LOOKUP[key] = comet
        COMET_LOOKUP[key.upper()] = comet

# Link up fragments to their parent comets
for comet in COMETS:
    if not comet.suffix: continue

    keys = comet.full_names(ignore_suffix=True)
    best = keys[0]

    # Identify or create a parent comet
    if best in BEST_NAME_LOOKUP:
        parent = BEST_NAME_LOOKUP[best]
    else:
        parent = comet.copy()
        parent.suffix    = ''
        parent.fragments = []
        parent.naif_id   = 0

        if VERBOSE:
            print('# Creating parent comet:', parent)

        # Update the dictionaries
        BEST_NAME_LOOKUP[best] = parent
        for key in keys:
            COMET_LOOKUP[key] = parent
            COMET_LOOKUP[key.upper()] = parent
            COMET_LOOKUP[key.upper().replace('-',' ')] = parent

    # Insert this fragment into the parent's list
    parent.fragments.append(comet)
    comet.parent = parent

    # Fill in any missing designations
    if parent.designation and not comet.designation:
        comet.designation = parent.designation
    if parent.alt_designations and not comet.alt_designations:
        comet.alt_designations = parent.alt_designations

# Creating parent comet: CometInfo(57|P|||duToit-Neujmin-Delporte|0|0|[])
# Creating parent comet: CometInfo(101|P|||Chernykh|0|0|[])
# Creating parent comet: CometInfo(205|P|||Giacobini|0|0|[])
# Creating parent comet: CometInfo(213|P|||Van Ness|0|0|[])
# Creating parent comet: CometInfo(332|P|||Ikeya-Murakami|0|0|[])
# Creating parent comet: CometInfo(0|C|C/1860 D1||Liais|0|0|[])
# Creating parent comet: CometInfo(0|C|C/1882 R1||Great September comet|0|0|[])
# Creating parent comet: CometInfo(0|C|C/1947 X1||Southern comet|0|0|[])
# Creating parent comet: CometInfo(0|C|C/1956 F1||Wirtanen|0|0|[])
# Creating parent comet: CometInfo(0|C|C/1965 S1||Ikeya-Seki|0|0|[])
# Creating parent comet: CometInfo(0|C|C/1969 O1||Kohoutek|0|0|[])
# Creating parent comet: CometInfo(0|C|C/1975 V1||West|0|0|[])
# Creating parent comet: CometInfo(0|C|C/1994 G1||Takamizawa-Levy|0|0|[])
# Creating parent comet: CometInfo(0|C|C/1996 J1||Evans-Drinkwater|0|0|[])
# Creating parent comet: CometInfo(0|C|C/2001 A2||LINEAR|0|0|[])
# Creating parent comet: CometInfo(0|C|C/2003 S4||LINEAR|0|0|[])
# Creating parent comet: CometInfo(0|P|P/2004 V5||LINEAR-Hill|0|0|[])
# Creating parent comet: CometInfo(0|C|C/2005 A1||LINEAR|0|0|[])
# Creating parent comet: CometInfo(0|P|P/2013 R3||Catalina-PanSTARRS|0|0|[])
# Creating parent comet: CometInfo(0|C|C/2015 E61||PanSTARRS|0|0|[])
# Creating parent comet: CometInfo(0|P|P/2016 J1||PanSTARRS|0|0|[])
# Creating parent comet: CometInfo(0|C|C/2020 P4|||0|0|[])

# A name followed by "+" returns a list of all the comets associated with that discoverer
# or team
for (name, comet_lists) in COMETS_BY_NAME.items():
    new_list = []
    for comet_list in comet_lists:
        first = comet_list[0]
        if first.parent:
            new_list.append(first.parent)
        else:
            new_list.append(first)

    BEST_NAME_LOOKUP[name + '+'] = new_list
    COMET_LOOKUP[name + '+'] = new_list
    COMET_LOOKUP[name.upper() + '+'] = new_list
    COMET_LOOKUP[name.upper().replace('-',' ') + '+'] = new_list

##########################################################################################

LID_REGEX       = re.compile(r'.*(comet\..*)', re.I)
PAREN_REGEX     = re.compile(r'(.*)\((.*)\)(.*)')
NO_SUFFIX_REGEX = re.compile(r'(.*?)-[A-Z]\d?')

def identify_comet(keys, warnings=[], ignore_suffix=False):
    """A CometInfo object, based on the identifications given."""

    # Find comet keys in big dictionary
    # Raise an error if multiple comets match
    if isinstance(keys, str):
        keys = [keys]

    original_keys = keys
    if comet_identifications.DEBUG:
        print('original comet keys:', keys)

    new_keys = []
    for key in keys:

        # This deals with a LID that includes something before "comet."
        match = LID_REGEX.fullmatch(key)
        if match:
            new_keys += [match.group(1)]
            continue

        # This deals with random ways of combining things inside parentheses
        match = PAREN_REGEX.fullmatch(key)
        if match:
            new_keys += [match.group(1).strip(),
                         match.group(2).strip(),
                         match.group(3).strip()]
            continue

        new_keys += [key]

    # Remove duplicates and blanks
    keys = []
    for key in new_keys:
        key = key.upper()

        if not key: continue
        if key in keys: continue

        keys.append(key)

    if ignore_suffix:
        stripped_keys = []
        for key in keys:
            match = NO_SUFFIX_REGEX.fullmatch(key)
            if match:
                key = match.group(1)

            stripped_keys.append(key)

        keys = stripped_keys

    if comet_identifications.DEBUG:
        print('processed comet keys:', keys)

    comet = None
    name = None     # A name alone could be ambiguous if missing an index
    for key in keys:
        key = key.upper()
        if key in COMET_LOOKUP:
            test_comet = COMET_LOOKUP[key]

            # Save the first match
            if not comet:
                comet = test_comet
                if comet_identifications.DEBUG:
                    print('comet found:', comet)

            # Subsequent matches must be consistent
            elif comet != test_comet:
              raise ValueError('Inconsistent comets designations: {comet}, {test_comet}')

        # Handle an ambiguous discoverer name string
        elif key + '+' in COMET_LOOKUP:
            test_name = COMET_LOOKUP[key + '+']

            # Save the first
            if not name:
                name = key
                if comet_identifications.DEBUG:
                    print('ambiguous comet name found:', name)

            # Subsequent matches must be consistent
            if name != key:
                raise ValueError('Inconsistent comet designations: {name}, {key}')

        # Handle an unrecognized key
        elif comet_identifications.IGNORE_EXTRA_NAMES:
            warnings.append(f'WARNING: Ignored comet identifier: {key}')
            if comet_identifications.DEBUG:
                print('comet identifier ignored:', key)

        else:
            raise ValueError(f'Unrecognized comet identifier: {key}')

    if not comet:
        if name:
            raise ValueError(f'Ambiguous comet name: {name}')
        else:
            raise KeyError(f'No comet found: {original_keys}')

    # If an ambiguous name was found (due to missing index), check it now
    if name:
        for c in COMET_LOOKUP[name + '+']:
            if c == comet:
                break

        if c != comet:
            raise ValueError(f'Name and comet do not match: {name}, {comet}')

    return comet

def append_comet_designations(info, comet):
    """Append the name and designations of this comet to a given list of identification
    tuples.

    Needed  to include cometary designations for hybrid objects that are also minor
    planets.

    The comet can be identified by a CometInfo object, a key, or a list of keys.
    """

    def merge_lists(list1, list2):

        # Find any PDS-specific designations, so we can keep them at the end
        early_names = []
        late_names = []
        for name in list1:
            if name.startswith('Minor Planet'):
                late_names.append(name)
            elif name.startswith('NAIF ID'):
                late_names.append(name)
            else:
                early_names.append(name)

        # Insert unique names into the list
        for name in list2:
            if name not in list1:
                if name.startswith('NAIF ID'):
                    late_names.append(name)
                else:
                    early_names.append(name)

        return early_names + late_names

    #### Begin active code

    if not isinstance(comet, CometInfo):
        comet = identify_comet(comet)

    comet_names = comet.full_names()

    # Append the comet names to each tuple in the list.
    # Right now, every identification tuple from the minor_planets module has unit length,
    # so this loop might not be strictly necessary.
    new_info = []
    for (mp_name, alt_designations, mp_type, desc, lid) in info:
        alt_designations = merge_lists(alt_designations, comet_names)

        desc.append('NOTE: This body is designated as both a minor planet and a comet.')

        new_info.append((mp_name, alt_designations, mp_type, desc, lid))

    return new_info

def comet_identifications(keys, a=0., e=0., i=0., q=0., warnings=[], ignore_suffix=False):

    # Identify the comet
    comet = identify_comet(keys, warnings=warnings, ignore_suffix=ignore_suffix)

    # If this is also a minor planet, return minor planet info with additional
    # designations. The minor_planet_identifications function will test a,e,i if
    # necessary.
    if comet.mp_number:
        # Set check_comets to False here because the comet has already been identified
        # and we are already prepared to append the additional designations.
        info = minor_planets.minor_planet_identifications(comet.mp_number,
                                                          a, e, i, q,
                                                          warnings,
                                                          check_comets=False)
        return append_comet_designations(info, comet)

    # Check orbital elements if necessary
    if a or e or i or q:
        (_, a0, e0, i0, q0) = minor_planets.get_mpc_info(comet.mpc_key())
        minor_planets.check_elements(a, e, i, q, a0, e0, i0, q0, warnings,
                                     comet.full_names()[0])

    return comet.target_identifications()

# These attributes of the function control its behavior...

# Indicates whether extraneous, unrecognized names should be ignored
comet_identifications.IGNORE_EXTRA_NAMES = False

# Prints out useful info for debugging
comet_identifications.DEBUG = False

##########################################################################################
# Export a few items to the minor_planets module now that comets is initialized
##########################################################################################

minor_planets.append_comet_designations = append_comet_designations
minor_planets.COMETS_NOW_MINOR_PLANETS = COMETS_NOW_MINOR_PLANETS

minor_planets.COMET_DESIGNATIONS = {}
for (number, letter, _, _, _, mp_number) in COMETS_NOW_MINOR_PLANETS:
    minor_planets.COMET_DESIGNATIONS[mp_number] = str(number) + letter

##########################################################################################
