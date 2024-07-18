##########################################################################################
# product_labels/__init__.py
##########################################################################################

import datetime
import fnmatch
import os
import shutil
import sys
from collections import defaultdict

import astropy.io.fits as pyfits
import pdslogger

from . import suffix_info
from .date_support            import (get_header_date,
                                      get_trl_timetags,
                                      merge_trl_timetags,
                                      get_label_retrieval_date,
                                      get_file_creation_date,
                                      set_file_timestamp)
from .hdu_data_descriptions   import fill_hdu_data_descriptions
from .hdu_dictionary_support  import fill_hdu_dictionary, repair_hdu_dictionaries
from .hst_dictionary_support  import fill_hst_dictionary
from .get_time_coordinates    import get_time_coordinates
from .nan_support             import cmp_ignoring_nans, has_nans, rewrite_wo_nans
from .wavelength_ranges       import wavelength_ranges
from .xml_support             import get_modification_history, get_target_identifications

from target_identifications import hst_target_identifications
from pdstemplate import PdsTemplate

LABEL_SUFFIX = '.xml'
DEBUG_DESCRIPTIONS = False

this_dir = os.path.split(suffix_info.__file__)[0]
template = this_dir + '/../templates/PRODUCT_LABEL.xml'
TEMPLATE = PdsTemplate(template)

# From https://archive.stsci.edu/hlsp/ipppssoot.html, valid last chars of the IPPPSSOOT
STANDARD_TRANSMISSION_TAILS = {'b', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't'}

# This needs to be updated the program is revised
LABEL_VERSION = '1.0'
timetag = max(os.path.getmtime(__file__), os.path.getmtime(template))
LABEL_DATE = datetime.datetime.fromtimestamp(timetag).strftime("%Y-%m-%d")

def label_hst_fits_directories(directories, root='', *,
                               match_pattern = '',
                               old_directories = [],
                               old_root = '',
                               retrieval_date = '',
                               logger = None,
                               reset_dates = True,
                               replace_nans = False):
    """Process one or more directories of HST FITS files, returning the information needed
    for all of their PDS4 labels as a dictionary keyed by the basenames.

    Input:
        directories         directory path or a list of directory paths.
        root                optional path to prepend to each directory path.
        match_pattern       optional fnmatch pattern to use to filter the files processed.
        old_directories     directory path containing the previous versions of the same
                            FITS files.
        old_root            optional path to prepend to each old directory path.
        retrieval_date      date the file was retrieved from MAST, in yyyy-mm-dd format.
                            If blank, the date will be retrieved from a pre-existing
                            label, or else it will be set to the file's creation date.
        logger              pdslogger to use; None for default EasyLogger.
        reset_dates         True to reset the modification date of each file to the date
                            found in the FITS header.
        replace_nans        True to rewrite each file without NaNs if NaNs are found.
    """

    filepaths = get_filepaths(directories, root, match_pattern)

    if old_directories:
        old_filepaths = get_filepaths(old_directories, old_root, match_pattern)
    else:
        old_filepaths = []

    label_hst_fits_filepaths(filepaths, root,
                             old_filepaths = old_filepaths,
                             old_root = old_root,
                             retrieval_date = '',
                             logger = logger,
                             reset_dates = reset_dates,
                             replace_nans = replace_nans)

############################################

def label_hst_fits_filepaths(filepaths, root='', *,
                             old_filepaths = [],
                             old_root = '',
                             retrieval_date = '',
                             logger = None,
                             reset_dates = True,
                             replace_nans = False):
    """Process a list of filepaths, returning the information needed for all of their
    PDS4 labels as a dictionary keyed by the basenames.

    Input:
        filepaths           a list of file paths.
        root                an optional directory path to prepend to each file path.
        old_filepaths       a list of file paths to old versions of these files.
        old_root            optional path to prepend to each old file path.
        retrieval_date      date the file was retrieved from MAST, in yyyy-mm-dd format.
                            If blank, the date will be inferred from the file itself.
        logger              pdslogger to use; None for default EasyLogger.
        replace_nans        True to rewrite each file without NaNs if NaNs are found.
        reset_dates         True to reset the modification date of each file to the date
                            found in the FITS header.
    """

    logger = logger or pdslogger.EasyLogger()

    # Find a common root if one was not provided
    # This is not strictly necessary but creates cleaner logs by suppressing the overall
    # common directory path
    if not root:
        try:
            min_filepath = min(filepaths)
            max_filepath = max(filepaths)
        except ValueError as e:
            logger.error(str(e) + ','.join(filepaths))

        for k, chars in enumerate(zip(min_filepath, max_filepath)):
            if chars[0] != chars[1]:
                break

        for j in range(k-1, -1, -1):
            if min_filepath[j] == '/':
                break

        root = min_filepath[:j+1]

    if root:
        logger.info('Root of file paths: ' + root)
        logger.replace_root(root)

    PdsTemplate.set_logger(logger)

    # Make sure the retrieval date, if any, is valid
    if retrieval_date:
        try:
            _ = datetime.date.fromisoformat(retrieval_date)
        except ValueError:
            logger.exception(ValueError)
            sys.exit(1)

    # Create a mapping from basename to old filepath
    old_fullpath_vs_basename = {os.path.basename(f):os.path.join(old_root,f)
                                for f in old_filepaths}

    # Save basic info about each FITS file
    info_by_basename = {}           # info vs. basename
    associations_by_ipppssoot = {}  # association list keyed by IPPPSSOOT
    trl_timetags_by_ipppssoot = {}  # time tag dictionary, which maps date to date-time

    prev_instrument_id = ''
    accepted_suffixes = set()
    for filepath in filepaths:

        # Ignore files that are not FITS
        if not filepath.endswith('.fits'):
            continue

        logger.info('Reading', filepath)

        fullpath = os.path.join(root, filepath)
        basename = os.path.basename(fullpath)
        ipppssoot_plus_suffix = basename.partition('.')[0]
        (ipppssoot, _, suffix) = ipppssoot_plus_suffix.partition('_')
        suffix = suffix.lower()

        # Determine the instrument
        instrument_id = suffix_info.INSTRUMENT_FROM_LETTER_CODE[basename[0]]

        if instrument_id != prev_instrument_id:
            logger.info('Instrument identified', instrument_id)
            accepted_suffixes = suffix_info.ACCEPTED_SUFFIXES[instrument_id]
            prev_instrument_id = instrument_id

        # If this suffix is not accepted, skip it
        if suffix not in accepted_suffixes:
            logger.warn(f'Suffix {suffix} rejected', filepath)
            continue

        ##################################################################################
        # Initialize the dictionary of file info keyed by basename, info_by_basename.
        #
        # "basename"       : basename of FITS file.
        # "filepath"       : path to FITS file as given in input.
        # "fullpath"       : full path to FITS file, including root path.
        # "ipppssoot"      : first nine letters of basename, before first underscore.
        # "suffix"         : suffix following first underscore, excluding ".fits".
        # "short_suffix"   : suffix without a trailing "_a", "_b", etc.
        # "lid_suffix"     : text to be appended to the IPPPSSOOT in the LID, e.g., "_a"
        #                    or "_1".
        # "group_ipppssoot": the IPPPSSOOT under which this file will be grouped; the last
        #                    character may differ from the actual IPPPSSOOT.
        # "retrieval_date" : the creation date on the file.
        ##################################################################################

        # Identify the retrieval date
        if retrieval_date:
            file_retrieval_date = retrieval_date
        else:
            # Take it from a pre-existing label, if any
            file_retrieval_date = get_label_retrieval_date(fullpath, LABEL_SUFFIX)

        if not file_retrieval_date:
            # Otherwise, use the file creation date. This is a bit dangerous because the
            # pipeline can modify the creation dates of files, meaning this date might be
            # wrong if you run the pipeline on the same directory a second time.
            file_retrieval_date = get_file_creation_date(fullpath)[:10]

        basename_dict = {
            'basename' : basename,
            'filepath' : filepath,
            'fullpath' : fullpath,
            'ipppssoot': ipppssoot,
            'suffix'   : suffix,
            'collection_name': suffix_info.collection_name(suffix, instrument_id),
            'lid_suffix'     : suffix_info.lid_suffix(suffix),
            'group_ipppssoot': ipppssoot,
            'retrieval_date' : file_retrieval_date,
            'label_version'  : LABEL_VERSION,
            'label_date'     : LABEL_DATE,
        }

        info_by_basename[basename] = basename_dict

        ##################################################################################
        # Gather info about the prior version, if any, of this FITS file
        #
        # "previous_fullpath"   : path to previous version of FITS file, or "".
        # "version_id"          : version_id for this product as a two-integer tuple.
        # "previous_xml"        : content old XML label as a single string.
        # "modification_history": list of modification history attributes from old label.
        # "fits_is_identical"   : True if this FITS file is identical to the old version.
        ##################################################################################

        previous_fullpath = old_fullpath_vs_basename.get(basename, '')
        if previous_fullpath:
            with open(previous_fullpath[:-5] + LABEL_SUFFIX) as f:
                xml_content = f.read()

            modification_history = get_modification_history(xml_content)
            old_version = modification_history[-1]['version_id']

            fits_is_identical = cmp_ignoring_nans(fullpath, previous_fullpath)
            if fits_is_identical:
                logger.info('Previous data is identical', filepath)
                version_id = (old_version[0], old_version[1]+1)
            else:
                logger.info('Previous data found', filepath)
                version_id = (old_version[0]+1, 0)
        else:
            version_id = (1, 0)
            xml_content = ''
            modification_history = []
            fits_is_identical = False

        # Update the dictionary
        basename_dict['previous_fullpath'   ] = previous_fullpath
        basename_dict['version_id'          ] = version_id
        basename_dict['previous_xml'        ] = xml_content
        basename_dict['modification_history'] = modification_history
        basename_dict['fits_is_identical'   ] = fits_is_identical

        ##################################################################################
        # Read fundamental info from a few specific files
        #
        # associations_by_ipppssoot[ipppssoot] = dictionary of associated IPPPSSOOTs.
        # trl_timetags_by_ipppssoot[ipppssoot] = dictionary mapping dates to date-times.
        ##################################################################################

        # Open the (original) FITS file
        original_path = fullpath + '-original'
        try:
            if os.path.exists(original_path):
                path = original_path
                hdulist = pyfits.open(original_path)
            else:
                path = fullpath
                hdulist = pyfits.open(fullpath)
        except OSError as e:
            logger.error(str(e) + path)
        # If this is an association file, read its contents for the IPPPSSOOT
        if suffix == 'asn':
            logger.info('Reading associations', filepath)
            associations = read_associations(hdulist, basename)
            associations_by_ipppssoot.update(associations)

        # If this is a TRL or PDQ file, get the dictionary of date-times vs. date
        if suffix in ('trl', 'pdq'):
            logger.info('Reading dates', filepath)
            timetag_dict = get_trl_timetags(hdulist[1], filepath, logger)
            if ipppssoot in trl_timetags_by_ipppssoot:
                merge_trl_timetags(trl_timetags_by_ipppssoot[ipppssoot], timetag_dict)
            else:
                trl_timetags_by_ipppssoot[ipppssoot] = timetag_dict

        ##################################################################################
        # Read structure and content info from this file
        #
        # "hdu_dictionaries": a list of dictionaries describing the content of each HDU.
        # "internal_date"   : the internal date, if any, from the first FITS header.
        # "has_nans"        : True if any data array in the file contains NaN.
        ##################################################################################

        # Gather the HDU structure info
        hdu_dictionaries = []
        for k,hdu in enumerate(hdulist):
            try:
                hdu_dictionaries.append(fill_hdu_dictionary(hdu, k, instrument_id,
                                                            filepath, logger))
            except (ValueError, TypeError) as e:
                logger.error('Irrecoverable error reading FITS file', filepath)
                logger.exception(e)
                break

        # Fix known errors
        repair_hdu_dictionaries(hdu_dictionaries, filepath, logger)

        # Update the dictionary
        basename_dict['hdu_dictionaries'] = hdu_dictionaries
        basename_dict['internal_date'   ] = get_header_date(hdulist)
        basename_dict['has_nans'        ] = has_nans(hdulist)

        hdulist.close()

    ######################################################################################
    # Define an alternative way to access the basename dictionaries:
    #   info_by_ipppssoot[ipppssoot][suffix] = basename_dict[basename]
    ######################################################################################

    info_by_ipppssoot = defaultdict(dict)
    for basename, basename_dict in info_by_basename.items():
        ipppssoot = basename_dict['ipppssoot']
        suffix = basename_dict['suffix']
        info_by_ipppssoot[ipppssoot][suffix] = basename_dict

    ######################################################################################
    # Report any old file path that isn't in the new retrieval, but should be
    ######################################################################################

    for basename in old_fullpath_vs_basename:
        ipppssoot_plus_suffix = basename.partition('.')[0]
        (ipppssoot, _, suffix) = ipppssoot_plus_suffix.partition('_')
        if ipppssoot in info_by_ipppssoot:
            if suffix not in info_by_ipppssoot[ipppssoot]:
                logger.error(f'Old filename {basename} is missing from new retrieval',
                             old_fullpath_vs_basename[basename])

    ######################################################################################
    # Annoyingly, some files match only by IPPPSSOO, not by IPPPSSOOT. Specifically, for
    # the some files such as jif, jit, cmh, cmi, and cmj, the character before the
    # underscore is "j". Some details are here:
    #   https://archive.stsci.edu/hlsp/ipppssoot.html
    # where the options for the final character are "b", "m", "n", "o", "p", "q", "r",
    # "s", are "t". (No mention of "j".)
    ######################################################################################

    # Organize IPPPSSOOTS by IPPPSSOO
    ipppssoots_from_ipppssoo = defaultdict(list)
    for ipppssoot, ipppssoot_dict in info_by_ipppssoot.items():
        if ipppssoot[-1].isnumeric():       # ignore associated files ending in a digit
            continue

        ipppssoots_from_ipppssoo[ipppssoot[:-1]].append(ipppssoot)

    # Save a mapping from actual IPPPSSOOT to standardized IPPPSSOOT
    group_ipppssoot_dict = {}

    # For each IPPPSSOO...
    for ipppssoo, ipppssoots in ipppssoots_from_ipppssoo.items():

        # If there's only one transmission character, there's nothing to do
        if len(ipppssoots) == 1:
            continue

        # Identify standard and non-standard transmission characters
        tails = {ipppssoot[-1] for ipppssoot in ipppssoots}
        standard_tails = STANDARD_TRANSMISSION_TAILS & tails
        nonstandard_tails = tails - standard_tails

        standard_tails = list(standard_tails)
        standard_tails.sort()

        nonstandard_tails = list(nonstandard_tails)
        nonstandard_tails.sort()

        # Select the reference transmission character, favoring standard values
        tails = standard_tails + nonstandard_tails
        reference_tail = tails[0]

        # Log an error if something happened that we don't understand
        if len(standard_tails) == 0:
            logger.error('No standard transmission characters found for ' +
                         f'"{ipppssoo}": {nonstandard_tails}')
        elif len(standard_tails) > 1:
            logger.error('Multiple standard transmission characters found for ' +
                         f'"{ipppssoo}": {standard_tails}')

            # Better to select the most common IPPPSSOOT
            tuples = []
            for tail in standard_tails:
                files_with_ipppssoot = len(info_by_ipppssoot[ipppssoo + tail])
                tuples.append((-files_with_ipppssoot, tail))

            tuples.sort()
            reference_tail = tuples[0][1]

        # Move each file under the reference IPPPSSOOT
        reference_ipppssoot = ipppssoo + reference_tail
        reference_dict = info_by_ipppssoot[reference_ipppssoot]
        for tail in tails:
            if tail == reference_tail:
                continue

            ipppssoot = ipppssoo + tail
            for suffix, suffix_dict in info_by_ipppssoot[ipppssoot].items():
                if suffix in reference_dict:
                    logger.error('Duplicated files with the same IPPPSSOO and suffix [1]',
                                 reference_dict[suffix]['filepath'])
                    logger.error('Duplicated files with the same IPPPSSOO and suffix [2]',
                                 suffix_dict['filepath'])
                    continue

                reference_dict[suffix] = suffix_dict
                reference_dict[suffix]['group_ipppssoot'] = reference_ipppssoot
                logger.info('Moved to standardized IPPPSSOOT ' + reference_ipppssoot,
                            suffix_dict['filepath'])
                group_ipppssoot_dict[ipppssoot] = reference_ipppssoot

            # If necessary, transfer the timetag info as well
            if ipppssoot in trl_timetags_by_ipppssoot:
                trl_timetags_by_ipppssoot[reference_ipppssoot] = \
                                                    trl_timetags_by_ipppssoot[ipppssoot]
                del trl_timetags_by_ipppssoot[reference_ipppssoot]

            del info_by_ipppssoot[ipppssoot]

    ######################################################################################
    # On rare occasions involving programs with multiple instruments, the IPPPSSOOT
    # of an "_asn" file can end in "0" but all the other files end in another digit.
    ######################################################################################

    solo_ipppssoots = []
    for ipppssoot, ipppssoot_dict in info_by_ipppssoot.items():

        suffixes = set(ipppssoot_dict.keys())
        if suffixes != {'asn'}:
            continue

        for k in range(0,9):
            alt_ipppssoot = ipppssoot[:-1] + str(k)
            if alt_ipppssoot == ipppssoot:
                continue
            if alt_ipppssoot not in info_by_ipppssoot:
                continue

            suffix_dict = ipppssoot_dict['asn']
            suffix_dict['group_ipppssoot'] = alt_ipppssoot
            info_by_ipppssoot[alt_ipppssoot]['asn'] = suffix_dict
            logger.info('Moved to standardized IPPPSSOOT ' + alt_ipppssoot,
                        suffix_dict['filepath'])
            solo_ipppssoots.append(ipppssoot)
            break

    for ipppssoot in solo_ipppssoots:
        del info_by_ipppssoot[ipppssoot]

    ######################################################################################
    # For all IPPPSSOOT dictionaries:
    #
    # "all_suffixes": set of all suffixes for this IPPPSSOOT.
    # "ipppssoot"   : IPPPSSOOT.
    # "timetags"    : timetags derived from the TRL and PDQ files.
    #
    # For both IPPPSSOOT and basename dictionaries:
    # "by_basename" : overall dictionary info_by_basename.
    # "by_ipppssoot": overall dictionary info_by_ipppssoot.
    #
    # Also, for basename dictionaries:
    # "ipppssoot_dict": dictionary for this file's IPPPSSOOT.
    #####################################################################################

# This isn't right, because the absence of a TRL file doesn't mean that we should
# completely ignore the associated data files.

#     removed_ipppssoot = []
#     for ipppssoot, ipppssoot_dict in info_by_ipppssoot.items():
#         all_suffixes = set(ipppssoot_dict.keys())
#         ipppssoot_dict['all_suffixes'] = all_suffixes
#         ipppssoot_dict['ipppssoot'] = ipppssoot
#         # ipppssoot_dict['timetags'] = trl_timetags_by_ipppssoot[ipppssoot]
#         try:
#             ipppssoot_dict['timetags'] = trl_timetags_by_ipppssoot[ipppssoot]
#         except KeyError:
#             # If trl timetags is missing, we remove this ipppssoot from info_by_ipppssoot
#             logger.error(f'Missing trl timetags for {ipppssoot}')
#             removed_ipppssoot.append(ipppssoot)
#
#         ipppssoot_dict['by_basename'] = info_by_basename
#         ipppssoot_dict['by_ipppssoot'] = info_by_ipppssoot
#
#         for suffix in all_suffixes:
#             basenamed_dict = ipppssoot_dict[suffix]
#             basenamed_dict['by_basename'] = info_by_basename
#             basenamed_dict['by_ipppssoot'] = info_by_ipppssoot
#             basenamed_dict['ipppssoot_dict'] = ipppssoot_dict
#
#     # Bypass this ipppssoot by removing them from info_by_ipppssoot & info_by_basename
#     for ipppssoot in removed_ipppssoot:
#         del info_by_ipppssoot[ipppssoot]
#         for basename in list(info_by_basename):
#             if ipppssoot in basename:
#                 del info_by_basename[basename]

    for ipppssoot, ipppssoot_dict in info_by_ipppssoot.items():
        all_suffixes = set(ipppssoot_dict.keys())
        ipppssoot_dict['all_suffixes'] = all_suffixes
        ipppssoot_dict['ipppssoot'] = ipppssoot
        ipppssoot_dict['by_basename'] = info_by_basename
        ipppssoot_dict['by_ipppssoot'] = info_by_ipppssoot

        for suffix in all_suffixes:
            basenamed_dict = ipppssoot_dict[suffix]
            basenamed_dict['by_basename'] = info_by_basename
            basenamed_dict['by_ipppssoot'] = info_by_ipppssoot
            basenamed_dict['ipppssoot_dict'] = ipppssoot_dict

        # If there's no TRL file, don't use a mapping from date to date-time
        if ipppssoot in trl_timetags_by_ipppssoot:
            ipppssoot_dict['timetags'] = trl_timetags_by_ipppssoot[ipppssoot]
        else:
            ipppssoot_dict['timetags'] = {}
            logger.warn(f'Missing trl timetags for {ipppssoot}')

    ######################################################################################
    # "associates"        : list of tuples (associated ipppssoot, memtype) for this
    #                       IPPPSSOOT if the associated IPPPSSOOT was found.
    # "missing_associates": list of tuples (associated ipppssoot, memtype) for this
    #                       IPPPSSOOT if the associated IPPPSSOOT _not_ found.
    # "parent"            : the "parent" IPPPSSOOT to which this IPPPSSOOT is associated,
    #                       if any.
    ######################################################################################

    # Make sure every we will have these items for every IPPPSSOOT
    for ipppssoot, ipppssoot_dict in info_by_ipppssoot.items():
        ipppssoot_dict['associates'] = []
        ipppssoot_dict['missing'] = []
        ipppssoot_dict['parent'] = ''

    # Insert the downward associations list from the ASN file
    for ipppssoot, associations in associations_by_ipppssoot.items():

        ippsoot_dict = info_by_ipppssoot[ipppssoot]

        # Warn if an associated file is missing
        validated_associates = []
        missing_associates = []
        for (associate, memtype) in associations:
            if associate in info_by_ipppssoot:
                validated_associates.append((associate, memtype))
                logger.debug(f'Associate "{memtype}" found for ' + ipppssoot, associate)

            else:
                missing_associates.append((associate, memtype))
                logger.warn('Missing associate for ' + ipppssoot, associate)

        ippsoot_dict['associates'] = validated_associates
        ippsoot_dict['missing_associates'] = missing_associates

        # Also record "parent" associations
        for (associate, _) in validated_associates:
            info_by_ipppssoot[associate]['parent'] = ipppssoot

    ######################################################################################
    # Identify the SPT/SHM/SHF file for each IPPPSSOOT.
    #
    # "spt_suffix"  : SPT/SHM/SHF suffix.
    # "spt_fullpath": Full path to the SPT file.
    ######################################################################################

    # Identify all the SPT files; make a list of IPPPSSOOTs without one
    ipppssoots_wo_spts = []
    for ipppssoot, ipppssoot_dict in info_by_ipppssoot.items():
        spt_suffixes = suffix_info.SPT_SUFFIXES & ipppssoot_dict['all_suffixes']
        if spt_suffixes:
            spt_suffix = spt_suffixes.pop()
            ipppssoot_dict['spt_suffix'] = spt_suffix
            ipppssoot_dict['spt_fullpath'] = ipppssoot_dict[spt_suffix]['fullpath']
        else:
            ipppssoots_wo_spts.append(ipppssoot)

    # If an IPPPSSOOT has no SPT, use that of one of its associates
    # This is OK because associates all use the same target and optical elements
    for ipppssoot in ipppssoots_wo_spts:
        ipppssoot_dict = info_by_ipppssoot[ipppssoot]
        spt_candidates = []
        for (associate, _) in ippsoot_dict['associates']:
            if 'spt_suffix' in info_by_ipppssoot[associate]:
                spt_candidates.append(associate)
                break

        if not spt_candidates:
            # There's no way for the pipeline to proceed from here
            logger.fatal('No SPT/SHM/SHF file found for ' + ipppssoot)
            raise IOError('No SPT/SHM/SHF file found for ' + ipppssoot)

        associate = spt_candidates[0]
        associate_dict = info_by_ipppssoot[associate]
        spt_suffix = associate_dict['spt_suffix']
        spt_fullpath = associate_dict[spt_suffix]['fullpath']
        ipppssoot_dict['spt_suffix'] = spt_suffix
        ipppssoot_dict['spt_fullpath'] = spt_fullpath

    ######################################################################################
    # Identify the reference file for each IPPPSSOOT. This is the file that serves as the
    # primary source of HST dictionary values.
    #
    # "reference_suffix"  : first reference suffix.
    # "reference_suffixes": set of all reference_suffixes.
    #
    # Then gather the metadata for each IPPPSSOOT:
    #
    # "hst_dictionary"        : HST dictionary.
    # "hst_proposal_id"       : HST proposal ID.
    # "instrument_id"         : instrument ID.
    # "instrument_name"       : full name of instrument.
    # "channel_id"            : channel ID.
    # "time_coordinates"      : tuple (start_time, stop_time)
    # "wavelength_ranges"     : list of wavelength_range values.
    # "target_identifications": list of Target_Identification tuples.
    ######################################################################################

    # Identify all reference files and generate the needed info
    no_reference_ipppssoots = []
    for ipppssoot, ipppssoot_dict in info_by_ipppssoot.items():
        all_suffixes = ipppssoot_dict['all_suffixes']

        # Check the REF_SUFFIXES list
        reference_suffixes = suffix_info.REF_SUFFIXES & all_suffixes
        tag = 'Reference'

        # If that fails, check the ALT_REF_SUFFIXES list
        if not reference_suffixes:
            reference_suffixes = suffix_info.ALT_REF_SUFFIXES & all_suffixes
            tag = 'Alternative reference'

        ipppssoot_dict['reference_suffixes'] = reference_suffixes

        # If there is no reference data file...
        if not reference_suffixes:

            # Maybe this is supposed to happen
            spt_hdulist = pyfits.open(ipppssoot_dict['spt_fullpath'])
            try:
                scidata = spt_hdulist[0].header['SCIDATA']
            except KeyError:
                scidata = None
            spt_hdulist.close()

            if scidata:
                logger.error('Reference file missing for ' + ipppssoot)
                no_reference_ipppssoots.append(ipppssoot)
                continue

            logger.warn('No science data for ' + ipppssoot)
            reference_suffix = ''
            ipppssoot_dict['reference_suffix'] = ''

        else:
            reference_suffix_list = list(reference_suffixes)
            reference_suffix_list.sort()
            reference_suffix = reference_suffix_list[0]
            count = len(reference_suffix_list)
            if count > 1:
                for k,suffix in enumerate(reference_suffix_list):
                    logger.debug(f'Multiple {tag.lower()} files found for {ipppssoot} ' +
                                 f'({k+1}/{count})', ipppssoot_dict[suffix]['filepath'])
            else:
                logger.debug(f'{tag} file found for ' + ipppssoot,
                             ipppssoot_dict[reference_suffix]['filepath'])

            ipppssoot_dict['reference_suffix'] = reference_suffix

        # Fill the HST dictionary
        if reference_suffix:
            fullpath = ipppssoot_dict[reference_suffix]['fullpath']
            ref_hdulist = pyfits.open(fullpath)
        else:
            fullpath = ipppssoot_dict['spt_fullpath']
            ref_hdulist = None

        spt_hdulist = pyfits.open(ipppssoot_dict['spt_fullpath'])
        hst_dictionary = fill_hst_dictionary(ref_hdulist, spt_hdulist, fullpath, logger)
        instrument_id = hst_dictionary['instrument_id']
        channel_id = hst_dictionary['channel_id']

        ipppssoot_dict['hst_dictionary'] = hst_dictionary
        ipppssoot_dict['hst_proposal_id'] = hst_dictionary['hst_proposal_id']
        ipppssoot_dict['instrument_id'] = instrument_id
        ipppssoot_dict['channel_id'] = channel_id
        ipppssoot_dict['instrument_name'] = suffix_info.INSTRUMENT_NAMES[instrument_id]

        # Time coordinates
        time_coordinates = get_time_coordinates(ref_hdulist, spt_hdulist, fullpath,
                                                                          logger)
        ipppssoot_dict['time_coordinates'] = time_coordinates[:2]
        ipppssoot_dict['time_is_actual'] = time_coordinates[2]

        # Wavelength ranges
        try:
            ranges = wavelength_ranges(hst_dictionary['instrument_id'],
                                       hst_dictionary['detector_ids'],
                                       hst_dictionary['filter_name'])
        except (ValueError, KeyError):
            ranges = ['UNK']
            logger.error('Undetermined wavelength range for ' +
                         hst_dictionary['filter_name'], fullpath)

        ipppssoot_dict['wavelength_ranges'] = ranges

        # Target identifications
        spt_fullpath = ipppssoot_dict['spt_fullpath']
        if reference_suffix:
            xml_content = ipppssoot_dict[reference_suffix]['previous_xml']
        else:
            xml_content = ''

        if xml_content:                 # use the old identification if available
            target_ids = get_target_identifications(xml_content)
        else:
            try:
                target_ids = hst_target_identifications(spt_hdulist[0].header,
                                                        spt_fullpath, logger)
                logger.debug('Target ' + str([rec[0] for rec in target_ids]),
                             ipppssoot_dict['spt_fullpath'])
            except (ValueError, KeyError) as e:
                target_ids = [('UNK', [], 'UNK', [], 'UNK')]
                logger.error(str(e).strip('"'), ipppssoot_dict['spt_fullpath'])
            except Exception as e:
                target_ids = [('UNK', [], 'UNK', [], 'UNK')]
                logger.exception(e)

        ipppssoot_dict['target_identifications'] = target_ids

        if ref_hdulist:
            ref_hdulist.close()
        spt_hdulist.close()

    # Delete IPPPSSOOTs without a reference file from structures
    no_reference_ipppssoots.sort()
    for ipppssoot in no_reference_ipppssoots:
        ipppssoot_dict = info_by_ipppssoot[ipppssoot]
        all_suffixes = list(ipppssoot_dict['all_suffixes'])
        all_suffixes.sort()
        for suffix in all_suffixes:
            basename = ipppssoot_dict[suffix]['basename']
            logger.error('File ignored due to missing reference file',
                         info_by_basename[basename]['fullpath'])
            del info_by_basename[basename]

        del info_by_ipppssoot[ipppssoot]

    ######################################################################################
    # For each basename, fill in some basics from the IPPPSSOOT dictionary
    #
    # "hst_dictionary"        : HST dictionary.
    # "instrument_id"         : instrument ID.
    # "instrument_name"       : full name of instrument.
    # "channel_id"            : channel ID.
    # "time_coordinates"      : tuple (start_time, stop_time)
    # "wavelength_ranges"     : list of wavelength_range values.
    # "target_identifications": list of Target_Identification tuples.
    #
    # Also, fill in some suffix-based info:
    #
    # "processing_level": processing_level ("Raw", "Calibrated", etc.)
    # "collection_title": collection_title
    # "browse_info"     : suffix-based browse_info.
    # "collection_lid"  : LID for this product.
    # "product_lid"     : LID for this product.
    # "product_lidvid"  : LIDVID for this product.
    ######################################################################################

    for basename, basename_dict in info_by_basename.items():
        for key in ('hst_dictionary', 'instrument_id', 'instrument_name', 'channel_id',
                    'time_coordinates', 'wavelength_ranges', 'target_identifications'):
            basename_dict[key] = basename_dict['ipppssoot_dict'][key]

        suffix = basename_dict['suffix']
        instrument_id = basename_dict['instrument_id']
        channel_id = basename_dict['channel_id']
        processing_level = suffix_info.get_processing_level(suffix, instrument_id,
                                                                    channel_id)
        collection_title_fmt = suffix_info.get_collection_title_fmt(suffix, instrument_id,
                                                                            channel_id)

        hst_dictionary = basename_dict['hst_dictionary']
        ipppssoot_dict = basename_dict['ipppssoot_dict']

        ic = instrument_id + ('/' + channel_id if channel_id != instrument_id else '')
        icp_dict = {'I': instrument_id, 'IC': ic, 'P': hst_dictionary['hst_proposal_id']}
        collection_title = collection_title_fmt.format(**icp_dict)

        basename_dict['processing_level'] = processing_level
        basename_dict['collection_title'] = collection_title

        collection_lid = ('urn:nasa:pds:hst_' + str(hst_dictionary['hst_proposal_id']) +
                          ':' + basename_dict['collection_name'])
        product_lid = (collection_lid + ':' +
                       basename_dict['ipppssoot'] + basename_dict['lid_suffix'])
        basename_dict['collection_lid'] = collection_lid
        basename_dict['product_lid'] = product_lid

        version_id = basename_dict['version_id']
        basename_dict['product_lidvid'] = (product_lid + '::' +
                                           str(version_id[0]) + ':' + str(version_id[1]))

        browse_info = suffix_info.BROWSE_SUFFIX_INFO[instrument_id, suffix]
        basename_dict['browse_info'] = browse_info

    ######################################################################################
    # Identify the associated Internal_Reference files for each basename or IPPPSSOOT +
    # suffix.
    #
    # "reference_list": set of tuples (group_ipppssoot, suffix) that will appear in the
    #                   Reference_List for this product.
    # "prior_pairs"   : set of tuples (group_ipppssoot, suffix) whose time tags must be
    #                   the same as or earlier than this file
    ######################################################################################

    # Make sure every dictionary will have these items
    for basename, basename_dict in info_by_basename.items():
        basename_dict['reference_list'] = set()
        basename_dict['prior_pairs'] = set()

    for ipppssoot, ipppssoot_dict in info_by_ipppssoot.items():
        instrument_id = ipppssoot_dict['instrument_id']

        # Create a set of tuples (ipppssoot, suffix) for every suffix
        ipppssoot_suffixes = ipppssoot_dict['all_suffixes']
        pairs = {(ipppssoot_dict[suffix]['group_ipppssoot'], suffix)
                 for suffix in ipppssoot_suffixes}

        # Add to the set the (ipppssoot, suffix) pair for every associate
        for (associate, memtype) in ipppssoot_dict['associates']:
            suffixes = info_by_ipppssoot[associate]['all_suffixes']
            pairs |= {(associate, suffix) for suffix in suffixes}

        # Organize these pairs by suffix
        pairs_by_suffix = defaultdict(set)
        for pair in pairs:
            pairs_by_suffix[pair[1]].add(pair)

        all_suffixes = set(pairs_by_suffix.keys())

        # The reference suffixes for this IPPPSSOOT are associated with all pairs
        # They must also be older than all other observational pairs with the same
        # IPPPSSOOT
        reference_suffixes = ipppssoot_dict['reference_suffixes']
        for suffix in reference_suffixes:
            ipppssoot_dict[suffix]['reference_list'] |= pairs
            for pair in pairs:
                if (pair[0] == ipppssoot
                    and pair[1] not in reference_suffixes
                    and suffix_info.is_observational(pair[1], instrument_id)):
                        ipppssoot_dict[pair[1]]['prior_pairs'].add((ipppssoot, suffix))

        # The TRL and PDQ files are newer than all of the observational files except raw
        # They point to every other file
        for suffix in {'trl', 'pdq'} & ipppssoot_suffixes:
            derived = {p for p in pairs
                       if suffix_info.get_processing_level(p[1], instrument_id)
                          not in ('Raw', 'Ancillary')}
            ipppssoot_dict[suffix]['reference_list'] |= pairs
            ipppssoot_dict[suffix]['prior_pairs'] |= derived

        # Every ASN file points to the reference files for this ipppssoot and associates
        # The ASN file must be newer than any of these files
        if 'asn' in ipppssoot_dict:
            reference_pairs = ({p for p in pairs if p[1] in reference_suffixes} |
                               {p for p in pairs if p[1] in suffix_info.REF_SUFFIXES})
            ipppssoot_dict['asn']['reference_list'] |= reference_pairs
            ipppssoot_dict['asn']['prior_pairs'] |= reference_pairs

        # Insert the prior suffixes from SUFFIX_INFO
        for suffix in ipppssoot_suffixes:
            prior_suffixes = suffix_info.get_prior_suffixes(suffix, instrument_id)
            prior_suffixes &= all_suffixes      # limit to suffixes that actually exist
            prior_pairs = {p for p in pairs if p[1] in prior_suffixes}

            suffix_dict = ipppssoot_dict[suffix]
            suffix_dict['reference_list'] |= prior_pairs
            suffix_dict['prior_pairs'] |= prior_pairs

    # Finalize the reference list
    for basename, basename_dict in info_by_basename.items():

        # Remove self
        me = (basename_dict['group_ipppssoot'], basename_dict['suffix'])
        basename_dict['reference_list'] -= {me}

        # Remove files with conflicting lid_suffixes
        lid_suffix = basename_dict['lid_suffix']
        if lid_suffix:
            excluded_suffixes = suffix_info.excluded_lid_suffixes(lid_suffix)
            new_set = {pair for pair in basename_dict['reference_list']
                       if suffix_info.lid_suffix(pair[1]) not in excluded_suffixes}
            basename_dict['reference_list'] = new_set

        # Make each association bi-directional
        for (associate, suffix) in basename_dict['reference_list']:
            group_ipppssoot = group_ipppssoot_dict.get(associate, associate)
            info_by_ipppssoot[group_ipppssoot][suffix]['reference_list'].add(me)

    ######################################################################################
    # "reference_basenames": sorted list of reference basenames
    ######################################################################################

    # Convert each reference_list to a sorted list of basenames
    for basename, basename_dict in info_by_basename.items():
        reference_basenames = []
        for (ipppssoot, suffix) in basename_dict['reference_list']:
            reference_basename = info_by_ipppssoot[ipppssoot][suffix]['basename']
            reference_basenames.append(reference_basename)
        reference_basenames.sort()
        basename_dict['reference_basenames'] = reference_basenames

    ######################################################################################
    # Fill in each modification date for each basename; make sure it is later than its
    # internal date and also later than the modification date of any of its priors.
    #
    #  "modification_date": best guess at a modification date in "yyyy-mm-ddThh:mm:ss"
    #                       format.
    ######################################################################################

    # Fill in each modification date; make sure it is no earlier than its internal date
    # and also later than the modification date of any of its priors.
    for basename, basename_dict in info_by_basename.items():
        fill_modification_dates(basename_dict, info_by_ipppssoot, logger)

    # Update the modification dates on the files if necessary
    if reset_dates:
        for basename, basename_dict in info_by_basename.items():
            modification_date = basename_dict['modification_date']
            if modification_date != '0000':     # year zero means date is unknown
                fullpath = basename_dict['fullpath']
                set_file_timestamp(fullpath, modification_date)
                logger.info('Modification date set to ' + modification_date, fullpath)

    ######################################################################################
    # Fill in the description fields for all data objects, i.e.,
    #   info_by_basename[basename]["hdu_dictionaries"][k]["data"]["description"]
    #
    # Also, fill in the product tile for each basename,
    #   info_by_basename[basename]["product_title"] = product_title
    #
    # If DEBUG_DESCRIPTIONS is True, the text of descriptions will be written into the
    # log, file by file, as a series of DEBUG messages. Otherwise, a sorted list of
    # unique descriptions is written to the log as DEBUG messages.
    ######################################################################################

    descriptions = set()
    for ipppssoot, ipppssoot_dict in info_by_ipppssoot.items():
        for suffix in ipppssoot_dict['all_suffixes']:
            descriptions |= fill_hdu_data_descriptions(ipppssoot, ipppssoot_dict,
                                                       suffix, DEBUG_DESCRIPTIONS, logger)

    if not DEBUG_DESCRIPTIONS:
        shortened = set()
        for desc in descriptions:
            if '(' in desc:
                before = desc.partition('(')[0]
                after = desc.partition(')')[2]
                desc = before + '(...)' + after
            shortened.add(desc)
        shortened = list(shortened)
        shortened.sort()
        logger.debug('Descriptions:\n    ' + '\n    '.join(shortened))

    ######################################################################################
    # Rewrite the file without NaNs if necessary
    ######################################################################################

    for basename, basename_dict in info_by_basename.items():
        basename_dict['nan_replacement'] = 0.
        basename_dict['hdus_with_nans'] = []

        if not basename_dict['has_nans']:
            continue

        fullpath = basename_dict['fullpath']
        if not replace_nans:
            logger.warn('NaNs found but not replaced', fullpath)
            continue

        original_path = fullpath + '-original'
        if os.path.exists(original_path):
            (nan_replacement,
             hdus_with_nans) = rewrite_wo_nans(original_path, rewrite=False)
            logger.info(f'NaNs already replaced with {nan_replacement}', fullpath)
        else:
            shutil.copy(fullpath, original_path)
            (nan_replacement,
             hdus_with_nans) = rewrite_wo_nans(fullpath, rewrite=True)
            logger.info(f'NaNs replaced with {nan_replacement}', fullpath)

        basename_dict['nan_replacement'] = nan_replacement
        basename_dict['hdus_with_nans'] = hdus_with_nans

    ######################################################################################
    # Write the new labels
    ######################################################################################

    for basename, basename_dict in info_by_basename.items():
        label_path = basename_dict['fullpath'].replace('.fits', LABEL_SUFFIX)
        TEMPLATE.write(basename_dict, label_path)
        if TEMPLATE.ERROR_COUNT == 1:
            logger.error('1 error encountered', label_path)
        elif TEMPLATE.ERROR_COUNT > 1:
            logger.error(f'{TEMPLATE.ERROR_COUNT} errors encountered', label_path)

##########################################################################################

def get_filepaths(directories, root='', match_pattern='', extension='.fits'):
    """Generate a list of file paths for processing.

    Input:
        directories     directory path of a list of directory paths.
        root            optional path to prepend to each directory path.
        match_pattern   optional fnmatch pattern to use to filter the returned list.
    """

    # Allow directories to be a single string
    if isinstance(directories, str):
        directories = [directories]

    # Prepare to strip off the root
    if root:
        root = root.rstrip('/') + '/'
        lroot = len(root)
    else:
        lroot = 0

    # Find file paths
    filepaths = []
    for directory in directories:
        fulldir = os.path.join(root, directory)
        for local_root, dirs, files in os.walk(fulldir):

            # Apply the match pattern if provided
            if match_pattern:
                files = [f for f in files if fnmatch.fnmatch(f, match_pattern)]

            # FITS files only
            files = [file for file in files if file.endswith(extension)]

            # Apply the necessary portion of the local root
            filepaths += [os.path.join(local_root[lroot:], file) for file in files]

    filepaths.sort()
    return filepaths

##########################################################################################

# Order of columns in an association table
MEMNAME = 0
MEMTYPE = 1

def read_associations(hdulist, basename):
    """Read the specified ASN file and return the list of associated IPPPSSOOT values in
    a dictionary keyed by each combined IPPPSSnnn. Each item in the list is a tuple
    (IPPPSSOOT value, memtype), where memtype has been truncated if it contains a suffix
    (following the dash).
    """

    table = hdulist[1].data

    # Identify all the products and their keys
    products = [rec for rec in table if rec[MEMTYPE].startswith('PROD')]
    if not products:
        raise IOError(f'No product MEMTYPES found in {basename}')

    if len(products) == 1:
        product = products[0][MEMNAME].lower()
        exposures = [(rec[MEMNAME].lower(), rec[MEMTYPE]) for rec in table
                     if not rec[MEMTYPE].startswith('PROD')]
        associations = {product: exposures}

    else:
        # products should be "PROD-<key>"
        product_by_key = {rec[MEMTYPE].partition('-')[2]:rec[MEMNAME].lower()
                          for rec in products}
        if '' in product_by_key:
            raise IOError(f'Inconsistent product MEMTYPES in {basename}; ' +
                          'should all be "PROD-<key>"')

        associations_by_key = {key:list() for key in product_by_key}
        for rec in table:
            (memtype, _, key) = rec[MEMTYPE].partition('-')
            if memtype == 'PROD':
                continue

            if key not in associations_by_key:
                raise IOError(f'MEMTYPE does not match product keys in {basename}: ' +
                              rec[MEMTYPE])

            associations_by_key[key].append((rec[MEMNAME].lower(), memtype))

        associations = {}
        for key, product in product_by_key.items():
            associations[product] = associations_by_key[key]

    return associations

##########################################################################################

def fill_modification_dates(basename_dict, info_by_ipppssoot, logger):
    """Fill in the modification date associated with this basename."""

    # Return if this modification date is already filled in
    if basename_dict.get('modification_date', ''):
        return

    # Determine the earliest date for this IPPPSSOOT from the TRL if any
    ipppssoot_dict = basename_dict['ipppssoot_dict']
    trl_dates = list(ipppssoot_dict['timetags'].values())
    earliest_trl_date = min(trl_dates) if trl_dates else '0000'     # year zero

    # If this file has no internal date, make sure its date is no earlier than the
    # earliest TRL date
    internal_date = basename_dict['internal_date']
    if not internal_date:
        internal_date = earliest_trl_date

    # Derived products should not be older than the earliest TRL date
    if basename_dict['processing_level'] not in ('Ancillary', 'Raw'):
        internal_date = max(internal_date, earliest_trl_date)

    # This avoids a possible infinite loop, even if blank
    basename_dict['modification_date'] = internal_date

    # Check for the modification dates of priors, recursively
    dates_found = [internal_date]
    for (prior_ipppssoot, prior_suffix) in basename_dict['prior_pairs']:
        prior_basename_dict = info_by_ipppssoot[prior_ipppssoot][prior_suffix]
        fill_modification_dates(prior_basename_dict, info_by_ipppssoot, logger)
        dates_found.append(prior_basename_dict['modification_date'])

    # Select the latest date
    latest_date = max(dates_found)

    # Fill in hh:mm:ss from TRL timetag dictionary if necessary and possible
    latest_date = ipppssoot_dict['timetags'].get(latest_date, latest_date)

    basename_dict['modification_date'] = latest_date
    if latest_date == '0000':
        logger.warn('Modification date unknown', basename_dict['fullpath'])
    else:
        logger.debug('Modification date = ' + latest_date, basename_dict['fullpath'])

##########################################################################################
