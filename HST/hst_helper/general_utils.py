##########################################################################################
# hst_helper/general_utils.py
##########################################################################################
import csv
import datetime
import os
import pdslogger

from hst_helper import (CITATION_INFO_DICT,
                        DOCUMENT_EXT,
                        INST_ID_DICT,
                        PROGRAM_INFO_FILE,
                        TARG_ID_DICT)
from hst_helper.fs_utils import (get_formatted_proposal_id,
                                 get_program_dir_path,
                                 get_instrument_id_from_fname)

from product_labels.xml_support import (get_modification_history,
                                        get_target_identifications)
from citations import Citation_Information
from xmltemplate import XmlTemplate


def create_xml_label(template_path, label_path, data_dict, logger):
    """Create xml label with given template path, label path, and data dictionary.

    Inputs:
        template_path:    the path of the label template.
        label_path:       the path of the label to be created.
        logger:           pdslogger to use; None for default EasyLogger.
    """
    logger = logger or pdslogger.EasyLogger()
    logger.info(f'Create label using template from: {template_path}.')
    TEMPLATE = XmlTemplate(template_path)
    XmlTemplate.set_logger(logger)

    logger.info('Insert data to the label template.')
    TEMPLATE.write(data_dict, label_path)
    if TEMPLATE.ERROR_COUNT == 1:
        logger.error('1 error encountered', label_path)
    elif TEMPLATE.ERROR_COUNT > 1:
        logger.error(f'{TEMPLATE.ERROR_COUNT} errors encountered', label_path)

def create_csv(csv_path, data, logger):
    """Create csv with given csv file path and data to be written into the csv file.

    Inputs:
        csv_path:    the path of the csv file.
        data:        a list that contains row data to be written into the csv file. Each
                     row data in the list is a list of column values for the row.
        logger:      pdslogger to use; None for default EasyLogger.
    """
    logger = logger or pdslogger.EasyLogger()
    logger.info(f'Create csv: {csv_path}.')


    # open the file in the write mode
    with open(csv_path, 'w') as f:
        # create the csv writer
        writer = csv.writer(f)

        # write rows to the csv file
        writer.writerows(data)

def get_citation_info(proposal_id, logger):
    """Search for proposal files & program info file stored at pipeline directory to
    obtain the citation info for a given proposal id.

    Inputs:
        proposal_id:    a proposal id.
        logger:         pdslogger to use; None for default EasyLogger.
    """
    formatted_proposal_id = get_formatted_proposal_id(proposal_id)
    if formatted_proposal_id in CITATION_INFO_DICT:
        return CITATION_INFO_DICT[formatted_proposal_id]
    logger = logger or pdslogger.EasyLogger()
    logger.info(f'Get citation info for: {proposal_id}.')

    pipeline_dir = get_program_dir_path(proposal_id, None, root_dir='pipeline')

    for file in os.listdir(pipeline_dir):
        _, _, ext = file.rpartition('.')
        if ext in DOCUMENT_EXT or file == PROGRAM_INFO_FILE:
            fp = pipeline_dir + f'/{file}'
            if formatted_proposal_id not in CITATION_INFO_DICT:
                CITATION_INFO_DICT[formatted_proposal_id] = (
                                                Citation_Information.create_from_file(fp))
            return CITATION_INFO_DICT[formatted_proposal_id]

def get_instrument_id_set(proposal_id, logger):
    """Walk through all downloaded files, store data in the INST_ID_DICT, and return
    a set of instrument ids for a given propsal id.

    Inputs:
        proposal_id:    a proposal id.
        logger:         pdslogger to use; None for default EasyLogger.
    """
    formatted_proposal_id = get_formatted_proposal_id(proposal_id)

    # Get imstrument id
    if formatted_proposal_id not in INST_ID_DICT:
        logger = logger or pdslogger.EasyLogger()
        logger.info(f'Get instrument ids for: {proposal_id}.')
        # Walk through all the downloaded files from MAST in staging directory
        files_dir = get_program_dir_path(proposal_id, None, root_dir='staging')
        for root, dirs, files in os.walk(files_dir):
            for file in files:
                inst_id = get_instrument_id_from_fname(file)
                formatted_proposal_id = str(proposal_id).zfill(5)
                if inst_id is not None:
                    INST_ID_DICT[formatted_proposal_id].add(inst_id)

    return INST_ID_DICT[formatted_proposal_id]

def get_mod_history_from_label(prev_label_path, current_version_id):
    """Compare current version id with the one in the exisiting (old) label. Return a
    list of modification history to be used in new label.

    Inputs:
        prev_label_path:      the path of the exisiting xml label
        current_version_id:   the current version id of the new bundle
    """
    mod_history = []
    if os.path.exists(prev_label_path):
        with open(prev_label_path) as f:
            xml_content = f.read()
            modification_history = get_modification_history(xml_content)
            old_version = modification_history[-1][1]
            if old_version != current_version_id:
                mod_history = modification_history

    return mod_history

def get_target_id_form_label(proposal_id, prev_label_path):
    """Get the target identification info from the exisitng label, if there is no
    existing label, walk through all downloaded files from Mast, store data in the
    TARG_ID_DICT, and return a lsit of target ids for a given propsal id.

    Inputs:
        proposal_id:          a proposal id.
        prev_label_path:      the path of the exisiting xml label
    """
    formatted_proposal_id = get_formatted_proposal_id(proposal_id)
    if formatted_proposal_id not in TARG_ID_DICT:
        if os.path.exists(prev_label_path):
            with open(prev_label_path) as f:
                xml_content = f.read()
                # use the old identification if available
                if xml_content:
                        target_ids = get_target_identifications(xml_content)
                        for targ in target_ids:
                            if targ not in TARG_ID_DICT[formatted_proposal_id]:
                                TARG_ID_DICT[formatted_proposal_id].append(targ)
        else:
            # TODO: might need to walk through bundles dir depending on if the files have
            # been moved to the bundles dir.
            # Walk through all the downloaded files from MAST in staging directory
            files_dir = get_program_dir_path(proposal_id, None, root_dir='staging')
            for root, _, files in os.walk(files_dir):
                for file in files:
                    if file.endswith('.xml'):
                        fp = os.path.join(root, file)
                        with open(fp) as f:
                            xml_content = f.read()
                            target_ids = get_target_identifications(xml_content)
                            for targ in target_ids:
                                if targ not in TARG_ID_DICT[formatted_proposal_id]:
                                    TARG_ID_DICT[formatted_proposal_id].append(targ)

    return TARG_ID_DICT[formatted_proposal_id]
