##########################################################################################
# hdu_data_descriptions.py
##########################################################################################

import pdslogger
from . import suffix_info
from .hst_dictionary_support import hst_detector_ids_for_file

# Description text format string based on EXTNAME.
# Not included are the EXTNAMEs used by the first data object (hdulist[1]), e.g., "SCI",
# "CGR", "DGR", "JIT", "SHH", "SHP", "TRL"; these descriptions are defined in the
# SUFFIX_INFO dictionary.

DESCRIPTIONS = {
    'ERR'     : 'Array containing the uncertainty (sigma) for each corresponding sample '+
                'in the {name} {noun}, HDU[{hdu}].',
    'DQ'      : 'Array containing the data quality mask as a bit-encoded array of flags '+
                'representing known status or problem conditions of the {name} {noun}, ' +
                'HDU[{hdu}].',
    'SAMP'    : 'Array containing the number of input samples used to compute each '     +
                'sample of the {name} {noun}, HDU[{hdu}].',
    'TIME'    : 'Array containing the effective integration time of each sample in the ' +
                'science {name} {noun}, HDU[{hdu}].',
    'EVENTS'  : 'Photon event list.',
    'GTI'     : 'Good time intervals table, indicating time intervals when events could '+
                'have been detected in the {name} {noun}, HDU[{hdu}].',
    'TIMELINE': 'Second by second values for spacecraft position, solar and target '     +
                'altitude, and count rates for background and the most prominent '       +
                'airglow lines, to use for interpreting the {name} {noun}, HDU[{hdu}].',
    'WHT'     : 'Weight image, giving the relative weight of pixels in the {name} '      +
                '{noun}, HDU[{hdu}].',
    'CTX'     : 'Context image, stored as a bit-encoded array that indicates which of '  +
                'the input images contributed to each pixel in the {name} {noun}, '      +
                'HDU[{hdu}].',
    'D2IMARR' : 'Array describing the filter-independent corrections for the CCD '       +
                'pixel-grid irregularities of this image.',
    'WCSDVARR': 'Array describing small-scale, filter-dependent distortions of this '    +
                'image.',
    'WCSCORR' : 'Table of the history of WCS changes for this image.',
    'HDRTAB'  : 'Complilation of header keywords that have a unique value for each '     +
                'each input image contributing to this combined image.',
    'HDRLET'  : 'List of independent WCS solutions, stored as a compact embedded FITS '  +
                'file embedded into a 1-D byte array.',
    'SNAP1'   : 'Instrument and telescope parameters collected immediately before the '  +
                'exposure began.',
    'SNAP2'   : 'Instrument and telescope parameters collected during the course of the '+
                'exposure.',

# These are just for NICMOS and I'm not sure what they mean
    'MNCIRSPD': 'NICMOS Cooling System (NCS) telemetry: CIrc Rot SPeeD.',
    'MNCORSPD': 'NICMOS Cooling System (NCS) telemetry: COmp Rot SPeeD.',
    'MNCONTRL': 'NICMOS Cooling System (NCS) telemetry: CONTRoL point temperature.',
    'MNPDSTPT': 'NICMOS Cooling System (NCS) telemetry: PiD SeT Point Temperature.',
    'MNCRADAT': 'NICMOS Cooling System (NCS) telemetry: Cpl RAD A Temperature.',
    'MNHTREJT': 'NICMOS Cooling System (NCS) telemetry: HeaT REJect Temperature.',
    'MNPNCOLT': 'NICMOS Cooling System (NCS) telemetry: Pri NiC OutLet Temperature.',
    'MNRNCILT': 'NICMOS Cooling System (NCS) telemetry: Red NiC InLet Temperature.',
    'MNRNCOLT': 'NICMOS Cooling System (NCS) telemetry: Red NiC OutLet Temperature.',
    'NDWTMP11': 'NICMOS telemetry: NIC1 mounting cup temperature.',
    'NDWTMP13': 'NICMOS telemetry: NIC3 mounting cup temperature.',
    'NDWTMP14': 'NICMOS telemetry: Cold well temperature.',
    'NDWTMP21': 'NICMOS telemetry: NIC2 cold mask temperature.',
    'NDWTMP22': 'NICMOS telemetry: NIC3 colds mast temperature.',
}

# These are needed for the waivered instruments where HDU[0] contains data.
EXTRA_DESCRIPTIONS = {
    'WFPC': {
        'C0H': 'Table of information about each of the 2-D calibrated images '    +
               'stacked into the 3-D data array of HDU[0].',
        'C1H': 'Table of information about each of the 2-D data quality masks '   +
               'stacked into the 3-D data array of HDU[0].',
        'C2H': 'Table of information about each of the 2-D histograms '           +
               'stacked into the 3-D data array of HDU[0].',
        'C3H': 'Table of information about each of the 2-D saturated pixel maps ' +
               'stacked into the 3-D data array of HDU[0].',
        'D0H': 'Table of information about each of the 2-D raw images '           +
               'stacked into the 3-D data array of HDU[0].',
        'Q0H': 'Table of information about each of the 2-D data quality masks '   +
               'stacked into the 3-D data array of HDU[0].',
        'Q1H': 'Table of information about each of the 2-D data quality masks '   +
               'stacked into the 3-D data array of HDU[0].',
        'X0H': 'Table of information about each of the 2-D engineering arrays '   +
               'stacked into the 3-D data array of HDU[0].',
        'SHH': 'Table of information about the header packets of HDU[0].',
    },
    'FGS': {
        'A1H': 'Table of information about the data in HDU[0].',
        'A2H': 'Table of information about the data in HDU[0].',
        'A3H': 'Table of information about the data in HDU[0].',
    },
    'FOC': {
        'C0H': 'Table of information about the data in HDU[0].',
        'C1H': 'Table of information about the data in HDU[0].',
        'D0H': 'Table of information about the data in HDU[0].',
        'Q0H': 'Table of information about the data in HDU[0].',
        'SHH': 'Table of information about the data in HDU[0].',
        'ULH': 'Table of information about the data in HDU[0].',
    },
    'FOS': {
        'C0H': 'Table of information about the data in HDU[0].',
        'C1H': 'Table of information about the data in HDU[0].',
        'C2H': 'Table of information about the data in HDU[0].',
        'C3H': 'Table of information about the data in HDU[0].',
        'C4H': 'Table of information about the data in HDU[0].',
        'C5H': 'Table of information about the data in HDU[0].',
        'C7H': 'Table of information about the data in HDU[0].',
        'CQH': 'Table of information about the data in HDU[0].',
        'D0H': 'Table of information about the data in HDU[0].',
        'Q0H': 'Table of information about the data in HDU[0].',
        'SHH': 'Table of information about the data in HDU[0].',
        'ULH': 'Table of information about the data in HDU[0].',
        'X0H': 'Table of information about the data in HDU[0].',
        'XQH': 'Table of information about the data in HDU[0].',
    },
    'HSP': {
        'C0H': 'Table of information about the data in HDU[0].',
        'C1H': 'Table of information about the data in HDU[0].',
        'C2H': 'Table of information about the data in HDU[0].',
        'C3H': 'Table of information about the data in HDU[0].',
        'D0H': 'Table of information about the data in HDU[0].',
        'D1H': 'Table of information about the data in HDU[0].',
        'D2H': 'Table of information about the data in HDU[0].',
        'D3H': 'Table of information about the data in HDU[0].',
        'Q0H': 'Table of information about the data quality mask in HDU[0].',
        'Q1H': 'Table of information about the data quality mask in HDU[0].',
        'Q2H': 'Table of information about the data quality mask in HDU[0].',
        'Q3H': 'Table of information about the data quality mask in HDU[0].',
        'SHH': 'Table of information about the header packets of HDU[0].',
        'ULH': 'Table of information about the unique data log in HDU[0].',
    },
    'GHRS': {
        'C0H': 'Table of information about the data in HDU[0].',
        'C1H': 'Table of information about the data in HDU[0].',
        'C2H': 'Table of information about the data in HDU[0].',
        'C3H': 'Table of information about the data in HDU[0].',
        'C4H': 'Table of information about the data quality mask in HDU[0].',
        'C5H': 'Table of information about the data in HDU[0].',
        'C5H': 'Table of information about the data in HDU[0].',
        'CQH': 'Table of information about the data quality mask in HDU[0].',
        'D0H': 'Table of information about the data in HDU[0].',
        'Q0H': 'Table of information about the data quality mask in HDU[0].',
        'SHH': 'Table of information about the header packets of HDU[0].',
        'ULH': 'Table of information about the unique data log in HDU[0].',
        'X0H': 'Table of information about the engineering data in HDU[0].',
        'XQH': 'Table of information about the data quality mask in HDU[0].',
    },
}

DATA_CLASS_TO_NOUN = {
    'Array_1D'         : ('data array', 'data arrays'),
    'Array_2D'         : ('data array', 'data arrays'),
    'Array_2D_Image'   : ('image', 'images'),
    'Array_3D_Image'   : ('image', 'images'),
    'Array_1D_Spectrum': ('spectrum', 'spectra'),
    'Array_2D_Spectrum': ('spectral image', 'spectral images'),
    'Table_Binary'     : ('table', 'tables'),
    ''                 : ('pseudo-array', 'pseudo-arrays'),
    'UNKNOWN'          : ('UNKNOWN DATA CLASS', 'UNKNOWN DATA CLASSES'),
}

def fill_hdu_data_descriptions(ipppssoot, ipppssoot_dict, suffix, log_text, logger):
    """Fill in the description fields for all of the data objects in one file based on its
    suffix.

    Input:
        ipppssoot       the file's IPPPSSOOT.
        ipppssoot_dict  the dictionary describing every file having the same IPPPSSOOT.
        suffix          the long suffix for one file.
        log_text        True to write all description strings to the log.
        logger          pdslogger to use; None to suppress logging.

    Return:             the set of all descriptions generated. Useful to review the text
                        being generated to confirm that it is correct.
    """

    logger = logger or pdslogger.NullLogger()

    hst_dictionary  = ipppssoot_dict['hst_dictionary']
    instrument_id   = hst_dictionary['instrument_id']
    channel_id      = hst_dictionary['channel_id']
    hst_proposal_id = hst_dictionary['hst_proposal_id']

    suffix_dict = ipppssoot_dict[suffix]
    filepath  = suffix_dict['fullpath']
    hdu_dicts = suffix_dict['hdu_dictionaries']

    description_set = set()

    # Get the list of detector_ids associated with this file. For COS, HSP, and FGS, this
    # might be a shorter list than the list associated with the IPPPSSOOT overall.
    detector_ids = hst_detector_ids_for_file(hst_dictionary, filepath)

    # Create a list indicating the detector IDs associated with each HDU.
    detector_ids_list = []      # HDU index -> list of detector IDs

    # For ACS/WFC, WFC3/UVIS, and WFPC2, the detector ID varies by HDU
    if 'detector_id_vs_extver' in hst_dictionary:
        detector_id_vs_extver = hst_dictionary['detector_id_vs_extver']
        for k, hdu_dict in enumerate(hdu_dicts):
            extver = hdu_dict['extver']
            if extver and extver in detector_id_vs_extver:
                detector_ids_list.append([detector_id_vs_extver[extver]])
            else:
                detector_ids_list.append(detector_ids)

    # Otherwise, each HDU applies to all detectors
    else:
        detector_ids_list = len(hdu_dicts) * [detector_ids]

    # Create a mapping from extver to the index of the relevant science HDU.
    science_hdu_index = {}      # extver -> index of first HDU with this extver
    for k, hdu_dict in enumerate(hdu_dicts):
        extver = hdu_dict['extver']
        if extver not in science_hdu_index:
            science_hdu_index[extver] = k

    # Some files have no EXTVER. When this happens, science_hdu_index[0] should refer to
    # the first HDU with an EXTNAME.
    if set(science_hdu_index.keys()) == {0}:
        for k, hdu_dict in enumerate(hdu_dicts):
            if hdu_dict['extname']:
                science_hdu_index[0] = k
                break

    # However, if the first HDU has data, it always has extver == 0.
    if hdu_dicts[0]['data']['data_class']:
        science_hdu_index[0] = 0

    # Determine the HDU index of the first data object class in this file
    if len(science_hdu_index) == 1:
        first_sci = list(science_hdu_index.values())[0]
    elif hdu_dicts[0]['data']['data_class']:
        first_sci = 0
    else:
        first_sci = science_hdu_index[1]

    sci_class = hdu_dicts[first_sci]['data']['data_class']
    sci_extname = hdu_dicts[first_sci]['extname']

    # Get the science data object description text for this file
    sci_fmt = suffix_info.get_hdu_description_fmt(suffix, instrument_id, channel_id)
    associated_suffix = suffix_info.get_associated_suffix(suffix, instrument_id,
                                                                  channel_id)

    # Identify the science objects in the associated file
    associated_hdu_index = {}   # extver -> index of first HDU with this extver
    associated_hdu_dicts = hdu_dicts
    associated_suffix_dict = suffix_dict
    if not associated_suffix:
        associated_hdu_index = science_hdu_index
        associated_sci_class = sci_class
    else:
        if isinstance(associated_suffix, str):
            selected_pairs = [pair for pair in suffix_dict['reference_list']
                              if pair[1] == associated_suffix]
        else:       # it's a set of candidate suffixes
            selected_pairs = [pair for pair in suffix_dict['reference_list']
                              if pair[1] in associated_suffix]

        if not selected_pairs:
            logger.warn(f'Missing associated suffix _{associated_suffix} for file',
                        filepath)
            associated_hdu_index = science_hdu_index
            associated_sci_class = sci_class
        else:
            (associated_ipppssoot, associated_suffix) = selected_pairs[0]
            associated_suffix_dict = (ipppssoot_dict['by_ipppssoot']
                                                    [associated_ipppssoot]
                                                    [associated_suffix])
            associated_hdu_dicts = associated_suffix_dict['hdu_dictionaries']
            for k, hdu_dict in enumerate(associated_hdu_dicts):
                extver = hdu_dict['extver']
                if extver not in associated_hdu_index:
                    associated_hdu_index[extver] = k

            try:
                first_sci = associated_hdu_index[1]
            except KeyError:
                first_sci = 0
            associated_sci_class = associated_hdu_dicts[first_sci]['data']['data_class']

    # Create dictionaries to fill in the description strings
    # These values are common across the entire file.
    (noun, nouns) = DATA_CLASS_TO_NOUN[sci_class]
    ic = instrument_id + ('/' + channel_id if channel_id != instrument_id else '')
    local_words = {
        'I'     : instrument_id,
        'P'     : hst_proposal_id,
        'IC'    : ic,
        'F'     : suffix_dict['basename'],
        'C'     : suffix_dict['collection_lid'],
        'lidvid': suffix_dict['product_lidvid'],
        'noun'  : noun,
        'nouns' : nouns,
    }

    (noun, nouns) = DATA_CLASS_TO_NOUN[associated_sci_class]
    associated_words = {
        'I'     : instrument_id,
        'P'     : hst_proposal_id,
        'IC'    : ic,
        'F'     : associated_suffix_dict['basename'],
        'C'     : associated_suffix_dict['collection_lid'],
        'lidvid': associated_suffix_dict['product_lidvid'],
        'noun'  : noun,
        'nouns' : nouns,
    }

    # Fill in the product title
    product_title_fmt = suffix_info.get_product_title_fmt(suffix, instrument_id,
                                                                  channel_id)
    suffix_dict['product_title'] = product_title_fmt.format(**associated_words)

    # Fill in the data description for each HDU.
    # If the data object is empty, update the header description instead.
    descriptions_for_log = []
    for k, hdu_dict in enumerate(hdu_dicts):
        extname = hdu_dict['extname']
        extver  = hdu_dict['extver']

        # Don't worry about the primary header if it has no data
        if not extname and not hdu_dicts[0]['data']['data_class']:
            continue            # use default description, already filled in

        # Select the dictionary to use for updating the format string
        if extname == sci_extname:
            description_fmt = sci_fmt
            hdu_words = associated_words
            science_hdu_dict = associated_hdu_dicts[associated_hdu_index[extver]]
        else:
            try:
                description_fmt = DESCRIPTIONS[extname]
            except KeyError:
                try:
                    description_fmt = EXTRA_DESCRIPTIONS[instrument_id][extname]
                except KeyError:
                    logger.error(f'Unrecognized EXTNAME "{extname}" in HDU[{k}]',
                                 filepath)
                    description_fmt = ''

            hdu_words = local_words
            science_hdu_dict = hdu_dicts[science_hdu_index[extver]]

        hdu_words['hdu'] = science_hdu_dict['index']
        hdu_words['name'] = science_hdu_dict['name']

        hdu_detector_ids = detector_ids_list[k]
        if len(hdu_detector_ids) == 1:
            hdu_words['D'] = hdu_detector_ids[0]
        else:
            hdu_words['D'] = ''

        # Create the data description
        description = description_fmt.format(**hdu_words)
        description = ' '.join(description.split())     # suppress repeated spaces
        hdu_dict['data']['description'] = description

        # If the data object is empty, update the header description instead.
        if hdu_dict['data']['is_empty']:
            parts = hdu_dict['header']['description'].split('\n')
            parts[0] = parts[0].rstrip('.') + ': ' + description
            hdu_dict['header']['description'] = '\n'.join(parts)

        # Prepare info for the log
        hdu_name = hdu_dict["name"]
        description_set.add(hdu_name + ': ' + description)
        descriptions_for_log.append(f'    HDU[{k}] ("{hdu_name}"): {description}')

    # One very long debug message, if requested
    if log_text and descriptions_for_log:
        logger.debug('New data descriptions', suffix_dict['fullpath']
                     + '\n' + '\n'.join(descriptions_for_log))

    return description_set

##########################################################################################
