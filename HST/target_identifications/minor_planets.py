##########################################################################################
# minor_planets module
#
# To import:
#   from minor_planets import minor_planet_identifications
#
# To use:
#   minor_planet_identifications(keys)
# returns a list of Target_Identification tuples. Each tuple contains:
#   (name, alt_designations, body_type, description, lid)
# where:
#   name                the preferred name
#   alt_designations    a list of strings indicating alternative names.
#   body_type           "Asteroid", "Centaur", etc.
#   description         a list of strings, to be separated by newlines inside
#                       the XML Target_Identification object.
#   lid                 the LID of the object, omitting "urn:...:target:"
#
# The input argument can be a single identification or else a list of one or more such
# identifications. If multiple identifications are provided, they must be internally
# consistent, all referring to the same minor planet. This can be used as an internal
# consistency check. For example,
#    minor_planet_identifications([2060, 'Chiron'])
# will return the info for the centaur Chiron, but
#    minor_planet_identifications([2061, 'Chiron'])
# will raise a ValueError, indicating that the two identifications are inconsistent.
#
# The reason it returns a list is for possible cases of associated objects, such as
# fragments of a single comet, that might all simultaneously be targets of a given
# observation. However, this capability is not yet used for minor planets.
#
# For additional internal checks, you can also call the function with explicit values for
# orbital elements a, e and i. If provided, an error will arise if any specified orbital
# element differs from that in the MPC by the specified threshhold of accuracy.
#
# For minor planets that were also classified as comets, the returned tuple includes the
# comet names as alternative designations. However, you cannot use a key that follows
# cometary nomenclature; use comets.comet_identifications() for that purpose.
#
# NOTE: This module performs real-time queries of the Minor Planet Center. As a result,
# expect a time delay of a few seconds for each identification.
##########################################################################################

import bs4
import math
import os
import re
import requests

from . import lids
from .OLD_STYLE_KBO_IDS import OLD_STYLE_KBO_IDS, NEW_STYLE_KBO_IDS

MINOR_PLANET_TYPES = ('Asteroid', 'Centaur', 'Trans-Neptunian Object', 'Dwarf Planet')

# The set of names also used for a moon or comet
DUPLICATED_NAMES = set([
    'Metis', 'Adrastea', 'Amalthea', 'Io', 'Europa', 'Leda', 'Pan', 'Prometheus',
    'Pandora', 'Epimetheus', 'Dione', 'Rhea', 'Cordelia', 'Ophelia', 'Bianca',
    'Desdemona', 'Titania', 'Galatea', 'Larissa', 'Halley',
])

# Manually update this with a list of minor planet numbers for dwarf planets
DWARF_PLANETS = [
         1,     # Ceres
    134340,     # Pluto
    136108,     # Haumea
    136199,     # Eris
    136472,     # Makemake
]

# Placeholders for the comets module
#
# This module cannot import comets, because comets imports minor_planets. The solution is
# that, when the comets module initializes, it defines these attributes of the
# minor_planets module.

append_comet_designations = None
COMETS_NOW_MINOR_PLANETS = None
COMET_DESIGNATIONS = None

##########################################################################################
# Classifications of minor planets
##########################################################################################

def mp_type_from_a(a):
    """Return the string "Asteroid", "Trans-Neptunian Object", or "Centaur"
    based on the semimajor axis in AU.
    """

    if a < 5.4:
        return 'Asteroid'

    if a > 30.1:
        return 'Trans-Neptunian Object'

    return 'Centaur'

##########################################################################################
# MPC access utilities
##########################################################################################

WEBCACHING = True
WEBCACHE = os.path.join(os.path.dirname(__file__), 'WEBCACHE')

URL_PREFIX = 'https://minorplanetcenter.net/db_search/show_object?object_id='

def get_mpc_info(key):
    """Get key information about this body from the Minor Planet Center. If the
    request fails, return None."""

    # Retrieve from cache if available
    filepath = os.path.join(WEBCACHE, key.upper().replace('/','-')) + '.html'
    if os.path.exists(filepath):
        with open(filepath, 'rb') as f:
            html = f.read()

    # Otherwise, retrieve from MPC
    else:
        url = URL_PREFIX + requests.utils.quote(str(key), safe='/')
        request = requests.get(url, allow_redirects=True)
        if request.status_code == 200:
            html = request.content
        else:
            html = b''

        if WEBCACHING:
            # Delete observation table because it can be huge
            parts = html.partition(b'<h2>Observations</h2>')
            if parts[2]:
                before_table = parts[2].partition(b'<table>')[0]
                after_table  = parts[2].partition(b'</table>')[2]
                html = parts[0] + parts[1] + before_table + after_table
            else:
                html = b''

            with open(filepath, 'wb') as f:
                f.write(html)

    if not html:
        return None

    soup = bs4.BeautifulSoup(html, 'html.parser')
    divs = soup.find_all('div')
    divs = [d for d in divs if 'id' in d.attrs and d.attrs['id'] == 'main']

    # Mal-formed pages indicate an unknown error
    if len(divs) == 0:
        raise ValueError('No main <div> for ' + key)

    if len(divs) > 1:
        raise ValueError('Multiple main <div>s for ' + key)

    # One sign of failure
    try:
        info = divs[0].h3.text
    except AttributeError:
        return None

    # Another sign of failure
    if info.strip().startswith('Data about'):
        return None

    #### Get names
    divs = soup.find_all('div')
    divs = [d for d in divs if 'id' in d.attrs and d.attrs['id'] == 'main']

    info = divs[0].h3.text

    # Clean up white space
    parts = info.split()
    info = ' '.join(parts)

    # Break into parts
    parts = info.split('=')
    names = [p.strip() for p in parts]

    #### Get orbital elements (or zero if not found)

    a = 0.
    e = 0.
    i = 0.
    q = 0.

    try:
        trs = soup.table.find_all('tr')
    except AttributeError:
        pass            # No orbit found
    else:
        for tr in trs[1:]:
            parts = tr.text.split('semimajor axis (AU)')
            if len(parts) > 1:
                try:
                    a = float(parts[1])
                except ValueError:
                    # For hyperbolic comets, semimajor axis is blank
                    pass

            parts = tr.text.split('eccentricity')
            if len(parts) > 1:
                e = float(parts[1])

            parts = tr.text.split('inclination (°)')
            if len(parts) > 1:
                i = float(parts[1])

            parts = tr.text.split('perihelion distance (AU)')
            if len(parts) > 1:
                q = float(parts[1])

    return (names, a, e, i, q)

##########################################################################################
# Class for minor planet information
##########################################################################################

class MinorPlanetInfo(object):

    ######## CONSTRUCTOR

    def __init__(self, body_type, number, designation, name,
                       alt_designations=[]):

        self.body_type = body_type
        self.number = int(number)
        self.designation = designation
        self.name = name
        self.alt_designations = alt_designations
        self.naif_id = 2000000 + self.number

        assert body_type in MINOR_PLANET_TYPES, \
                             'Invalid minor planet type: ' + body_type

        if not designation and alt_designations:
            self.designation = alt_designations[0]
            self.alt_designations = alt_designations[1:]

        assert self.number or self.designation or self.name, 'Unnamed body'

    ######## UTILITIES

    def __str__(self):
        return ('MinorPlanetInfo(' + self.body_type    + '|' +
                                     str(self.number)  + '|' +
                                     self.designation  + '|' +
                                     self.name         + '|' +
                                     str(self.alt_designations) + ')')

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):

        # Matching numbers confirm equality
        if self.number and other.number:
            return self.number == other.number

        # Matching names confirm equality
        if self.name and other.name:
            return self.name.lower() == other.name.lower()

        # One matching designation confirms equality
        self_designations  = set([self.designation]  + self.alt_designations)
        other_designations = set([other.designation] + other.alt_designations)
        overlap = self_designations.intersection(other_designations)
        if overlap:
            return True

        return False

    def merge(self, other):

        self.number = self.number or other.number
        self.name   = self.name   or other.name

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

        clone = MinorPlanetInfo(self.number, self.designation, self.name,
                                self.alt_designations)

        if deep:
            clone.alt_designations = list(clone.alt_designations)

        return clone

    def full_names(self):

        names = []
        if self.number and self.name:
            names += [f'{self.number} {self.name}',
                      f'({self.number}) {self.name}']

        designations = [self.designation] + self.alt_designations
        old_style_designations = []
        for designation in designations:
            if not designation:
                continue

            if self.number:
                names += [f'({self.number}) {designation}']

            names += [designation]

            if designation in OLD_STYLE_KBO_IDS:
                old_style_designations.append(OLD_STYLE_KBO_IDS[designation])

        names += old_style_designations

        if self.name and self.name not in DUPLICATED_NAMES:
            names += [self.name]

        if self.number:
            names += [f'Minor Planet {self.number}',
                      f'NAIF ID {self.naif_id}']

        return names

    def lid(self):

        if self.number:
            if self.name:
                string = f'{self.number}_{self.name}'
            else:
                string = f'{self.number}_{self.designation}'
        else:
            string = self.designation

        return lids.clean(self.body_type + '.' + string)

    def target_identifications(self):

        names = self.full_names()
        return [(names[0], names[1:], self.body_type, [], self.lid())]

    def mpc_key(self):

        if self.number:
            return str(self.number) + self.letter
        else:
            return self.designation

##########################################################################################

NUMBER      = r'\(?(\d+)\)?'
NAME        = r"([A-Z`'][A-Z `'\.\|!-]*[A-Z])"
DESIGNATION = r'(\d{4} *[A-Z]{1,2}\d*)'

NUMBER_REGEX             = re.compile(NUMBER)
NAME_REGEX               = re.compile(NAME, re.I)
NUMBER_NAME_REGEX        = re.compile(NUMBER + r' *' + NAME, re.I)
NUMBER_DESIGNATION_REGEX = re.compile(NUMBER + r' *' + DESIGNATION, re.I)
DESIGNATION_REGEX        = re.compile(DESIGNATION)

LID_REGEX = re.compile(r'.*(asteroid|centaur|trans-neptunian_object)\.(.*)')
NAIF_ID_REGEX = re.compile(r'NAIF ID (2\d{6})', re.I)
MINOR_PLANET_REGEX = re.compile(r'Minor Planet (\d+)', re.I)

def minor_planet_identifications(keys, a=0., e=0., i=0., q=0.,
                                 warnings=[], check_comets=True):

    # Note that check_comets is just used for communication with the comets
    # module when identifying a hybrid body. Set to False to prevent the comets
    # module from identifying the body a second time.

    if isinstance(keys, (str,int)):
        keys = [keys]

    original_keys = keys
    if minor_planet_identifications.DEBUG:
        print('original minor planet keys:', keys)

    # Handle cases of keys that are longer than necessary, and/or match
    # PDS-specific alternative designations
    new_keys = []
    for key in keys:
        if isinstance(key, int):
            new_keys += [str(key)]
            continue

        # If this is in the format of a LID, clean it up and continue trying...
        match = LID_REGEX.fullmatch(key)
        if match:
            key = match.group(2).upper().replace('_', ' ')

        match = DESIGNATION_REGEX.fullmatch(key)
        if match:
            new_keys += [match.group(1)]
            continue

        match = NUMBER_REGEX.fullmatch(key)
        if match:
            new_keys += [match.group(1)]
            continue

        match = NUMBER_NAME_REGEX.fullmatch(key)
        if match:
            new_keys += [match.group(1), match.group(2)]
            continue

        match = NUMBER_DESIGNATION_REGEX.fullmatch(key)
        if match:
            new_keys += [match.group(1), match.group(2)]
            continue

        match = NAIF_ID_REGEX.fullmatch(key)
        if match:
            new_keys += [match.group(1)]
            continue

        match = MINOR_PLANET_REGEX.fullmatch(key)
        if match:
            new_keys += [match.group(1)]
            continue

        match = NAME_REGEX.fullmatch(key)
        if match:
            new_keys += [match.group(1)]
            continue

        if key in NEW_STYLE_KBO_IDS:
            new_keys += [NEW_STYLE_KBO_IDS[key]]
            continue

        # Handle an unrecognized key
        if minor_planet_identifications.IGNORE_EXTRA_NAMES:
            warnings.append(f'Malformed minor planet identifier, ignored: {key}')
            if minor_planet_identifications.DEBUG:
                print('malformed minor planet identifier, ignored:', key)
        else:
            raise ValueError(f'Unrecognized minor planet identifier: {key}')

    if not new_keys:
        raise KeyError(f'No matching minor planets: {keys}')

    # Remove any duplicates
    keys = []
    for key in new_keys:
        key = key.upper()
        if key not in keys:
            keys.append(key)

    if minor_planet_identifications.DEBUG:
        print('processed minor planet keys:', keys)

    # Look up at MPC, gather names
    names = None
    for key in keys:
        info = get_mpc_info(key)

        # Name not found at MPC
        if not info:
            if minor_planet_identifications.IGNORE_EXTRA_NAMES:
                warnings.append(f'Unrecognized MPC identifier, ignored: {key}')
                if minor_planet_identifications.DEBUG:
                    print('unrecognized MPC identifier, ignored:', key)
                continue
            else:
                raise ValueError(f'Unrecognized MPC identifier: {key}')

        (test_names, a0, e0, i0, q0) = info
        if info and minor_planet_identifications.DEBUG:
            print('MPC record found:', key, test_names)

        # Save the first match
        if not names:
            names = test_names
            if minor_planet_identifications.DEBUG:
                print('minor planet identified:', key)

        # Subsequent names must match the first
        elif names != test_names:
            raise ValueError('Inconsistent minor planet designations: ' +
                             str(original_keys))

    if not names:
        raise KeyError(f'No matching minor planets: {keys}')

    # "(number) Name"?
    match = NUMBER_NAME_REGEX.fullmatch(names[0])
    if match:
        number = int(match.group(1))
        name   = match.group(2)
        designations = names[1:]
    else:
      match = NUMBER_REGEX.fullmatch(names[0])
      if match:
        # "(number)"?
        number = int(match.group(1))
        name   = ''
        designations = names[1:]
      else:
        # yyyy xxnn
        number = 0
        name   = ''
        designations = names

    if minor_planet_identifications.DEBUG:
        print('minor planet number:', number)
        print('minor planet name:', name)
        print('minor planet designations:', designations)

# Skip this; some "official" MPC identifiers are not in the stardard form
#     for designation in designations:
#         match = DESIGNATION_REGEX.fullmatch(designation)
#         if not match:
#             raise ValueError(f'Invalid designation: {designation}')

    # Define the type
    if number and number in DWARF_PLANETS:
        body_type = 'Dwarf Planet'
    else:
        body_type = mp_type_from_a(a0)

    if minor_planet_identifications.DEBUG:
        print(f'minor planet body type: {body_type}, (a = {a0})')

    # Define the MinorPlanetInfo; include fixes for Pluto
    if name == 'Pluto':
        mp = MinorPlanetInfo(body_type, number, '', name, [])
        mp.naif_id = 999
    elif len(designations) == 0:
        mp = MinorPlanetInfo(body_type, number, '', name, [])
    elif len(designations) == 1:
        mp = MinorPlanetInfo(body_type, number, designations[0], name, [])
    else:
        mp = MinorPlanetInfo(body_type, number, designations[0], name,
                             designations[1:])

    # Check orbital elements if necessary
    check_elements(a, e, i, q, a0, e0, i0, q0, warnings, (mp.name or mp.designation
                                                                  or str(mp.number)))

    # Get the Target Identification tuples
    info = mp.target_identifications()

    # For minor planets once designated as comets, fill in the complete list of
    # alt_designations
    if check_comets and COMET_DESIGNATIONS and mp.number in COMET_DESIGNATIONS:
        info = append_comet_designations(info, COMET_DESIGNATIONS[mp.number])

    return info

# These attributes of the function control its behavior...

# Indicates whether extraneous, unrecognized names should be ignored
minor_planet_identifications.IGNORE_EXTRA_NAMES = False

# Prints out useful info for debugging
minor_planet_identifications.DEBUG = False

##########################################################################################

THRESHOLD = 0.05

def check_elements(a, e, i, q, a0, e0, i0, q0, warnings=[], name=''):

    # Check orbital elements only if necessary
    if (a,e,i,q) == (0.,0.,0.,0.):
        return

    # Orbital elements not available
    if (a0,e0,i0,q0) == (0.,0.,0.,0.):
        if name:
            warnings.append(f'MPC orbital elements are not available for "{name}"')
        else:
            warnings.append('MPC orbital elements are not available')
        return

    values = []
    errors = []
    if a:
        errors.append(abs(a - a0) / max(a, a0))     # relative error
        values += [(a, a0)]

    if q:
        errors.append(abs(q - q0) / max(q, q0))     # relative error
        values += [(q, q0)]

    if e:
        errors.append(abs(e - e0))                  # absolute error
        values += [(e, e0)]

    if i:
        errors.append(abs(i - i0) * math.pi/180.)   # absolute, in radians
        values += [(i, i0)]

    errors.sort()
    if errors[-2] > THRESHOLD:
        if name:
            raise ValueError(f'Orbital element mismatch for "{name}": {values}')
        else:
            raise ValueError(f'Orbital element mismatch: {values}')

    if errors[-1] > THRESHOLD:
        if name:
            warnings.append(f'One orbital element mismatch for "{name}": {values}')
        else:
            warnings.append(f'One orbital element mismatch: {values}')

##########################################################################################
