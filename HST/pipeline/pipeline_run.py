#!/usr/bin/env python3
##########################################################################################
# pipeline/pipeline_run.py
#
# Syntax:
# pipeline_run.py [-h] [--proposal-ids PROPOSAL_IDS [PROPOSAL_IDS ...]]
#                 [--log LOG] [--quiet]
#                 [--max-subproc-cnt MAX_SUBPROC_CNT]
#                 [--max-allowed-time MAX_ALLOWED_TIME] [--get-ids]
#
# Enter the --help option to see more information.
#
# The script to start hst pipeline process for the given proposal ids.
##########################################################################################

import argparse
import datetime
import os
import pdslogger
import sys

from hst_helper import HST_DIR
from organize_files import clean_up_staging_dir
from query_hst_moving_targets import query_hst_moving_targets
from queue_manager import run_pipeline
import queue_manager

# Set up parser
parser = argparse.ArgumentParser(
    description="""pipeline_run: run pipeline with the default list of ids or the given
                proposal ids.""")

parser.add_argument('--proposal-ids', nargs='+', type=str, default='',
    help='The proposal id for the MAST query.')

parser.add_argument('--log', '-l', type=str, default='',
    help="""Path and name for the log file. The name always has the current date and time
         appended. If not specified, the file will be written to the current logs
         directory and named "run-pipeline-<date>.log".""")

parser.add_argument('--quiet', '-q', action='store_true',
    help='Do not also log to the terminal.')

parser.add_argument('--max-subproc-cnt', '--max-subproc',
    type=int, action='store', default=20,
    help='Max number of subprocesses to run at a time for one pipeline process.')

parser.add_argument('--max-allowed-time', '--max-time',
    type=int, action='store', default=1800,
    help='Max allowed subprocess time in seconds before it gets killed.')

parser.add_argument('--get-ids', '-g', action='store_true',
    help='Fetch ids from MAST to update the id list before running pipeline.')

# Default list of program ids
ids_li = ['15648', '13667', '10161', '11113', '08152', '15142', '11650', '04600',
          '10423', '15929', '09354', '11573', '10719', '08699', '07430', '07583',
          '12004', '15972', '11806', '09747', '11556', '07443', '05437', '09386',
          '16040', '05168', '09713', '11644', '16995', '11969', '11178', '09738',
          '15456', '06259', '07428', '16183', '15502', '06679', '10862', '06024',
          '15925', '17292', '13668', '05590', '15500', '13633', '12045', '08258',
          '08274', '07324', '05642', '06793', '05836', '09235', '14334', '16057',
          '13937', '16313', '13631', '15625', '11990', '16009', '07769', '15360',
          '16913', '08391', '11498', '15674', '12891', '04225', '14053', '05424',
          '07317', '13315', '16704', '12725', '09745', '06758', '02483', '16053',
          '05217', '07885', '02145', '10156', '13051', '06559', '05466', '12041',
          '08229', '11115', '15143', '09975', '05831', '03214', '06009', '14092',
          '15261', '14195', '15262', '12535', '06818', '15357', '16074', '08169',
          '10514', '10770', '15171', '15706', '16266', '16088', '05832', '06852',
          '11972', '14103', '14661', '12244', '08634', '12468', '16929', '16878',
          '15505', '03178', '12801', '07276', '11169', '06447', '07313', '08583',
          '10512', '15122', '12049', '15460', '07315', '16306', '15344', '12245',
          '11102', '11226', '08661', '15379', '13679', '13864', '17291', '13118',
          '15144', '04442', '09823', '10165', '10534', '15425', '17294', '11292',
          '14864', '09052', '11998', '14650', '13612', '05837', '09060', '11984',
          '16056', '12897', '05361', '16501', '06630', '10786', '16454', '16280',
          '08580', '11190', '07178', '15097', '10870', '07858', '10502', '11118',
          '14524', '09035', '12237', '11187', '12665', '16310', '12176', '06452',
          '15259', '16880', '13694', '09384', '13492', '06509', '14928', '13414',
          '06692', '06591', '12894', '07239', '07589', '14931', '13198', '16674',
          '12312', '17206', '09426', '08680', '10102', '06013', '07427', '15372',
          '13934', '17275', '05828', '14475', '10992', '14261', '01478', '11310',
          '05452', '14192', '13692', '14756', '17187', '12305', '14634', '16452',
          '09344', '09341', '13693', '09678', '05823', '10083', '06025', '14138',
          '05313', '04764', '04580', '09393', '14798', '10545', '13873', '15249',
          '04976', '06752', '05489', '10799', '05167', '13031', '11957', '06753',
          '11156', '14105', '15595', '15665', '15419', '05834', '13396', '10268',
          '14811', '11559', '05392', '12537', '13404', '11361', '10774', '14492',
          '10782', '11314', '14263', '04771', '05843', '06219', '11096', '06008',
          '06853', '14136', '09809', '15965', '10805', '16084', '14113', '10800',
          '06634', '06741', '07180', '06295', '06613', '05508', '16790', '13328',
          '10398', '02231', '14612', '12243', '16904', '05653', '06743', '13664',
          '04662', '16663', '08148', '06662', '14217', '16499', '11109', '16795',
          '17289', '05757', '05829', '11055', '07823', '12597', '07429', '09440',
          '04521', '13067', '17254', '03365', '16315', '16087', '17142', '06733',
          '14752', '14839', '06854', '11970', '10422', '07616', '08657', '14474',
          '10065', '12980', '15328', '15622', '06028', '17163', '07321', '08108',
          '05348', '12436', '09385', '09056', '09110', '04748', '15207', '08398',
          '08802', '15978', '13713', '10871', '06806', '16924', '11418', '11630',
          '14627', '08577', '15424', '13609', '08224', '13035', '10140', '09725',
          '09185', '05329', '10170', '08658', '14790', '12234', '07319', '02560',
          '06497', '15233', '05662', '04784', '14064', '05221', '16309', '07822',
          '10427', '13619', '10341', '08539', '17297', '06715', '09119', '16041',
          '08660', '07444', '04886', '04113', '08405', '16043', '03707', '03803',
          '09748', '03069', '06663', '07308', '03887', '10473', '06755', '07243',
          '12003', '15342', '06736', '16993', '03744', '10508', '14491', '12463',
          '05493', '04549', '06141', '09906', '07916', '04966', '10557', '13794',
          '15638', '07223', '10192', '04601', '02627', '04005', '14133', '13488',
          '05440', '14884', '06481', '01290', '16044', '05218', '03373', '13716',
          '14755', '13675', '16199', '17089', '13610', '04204', '02442', '15481',
          '03848', '11085', '03997', '02481', '12900', '13438', '16781', '03862',
          '15973', '12601', '13865', '14616', '11518', '11971', '09904', '05769',
          '08171', '05591', '13311', '17293', '14044', '06513', '10144', '02432',
          '07244', '11536', '13795', '10456', '04548', '01206', '07582', '05830',
          '05215', '12435', '15405', '13620', '05414', '05624', '17330', '12053',
          '12675', '05835', '08871', '05640', '13197', '05844', '13503', '17302',
          '16463', '15248', '06648', '04015', '05391', '04001', '13829', '13474',
          '01285', '07939', '15108', '12433', '13055', '05824', '16917', '06650',
          '16418', '16308', '14045', '04353', '09320', '11170', '16652', '01289',
          '06841', '15421', '09585', '05417', '11656', '10357', '06218', '05020',
          '16167', '10783', '14040', '15406', '16073', '03644', '03186', '13712',
          '06328', '16687', '13005', '13229', '03319', '04024', '14499', '06521',
          '06027', '12077', '08117', '16192', '16077', '16496', '10507', '11645',
          '11566', '15868', '09496', '09391', '15158', '06786', '07179', '15343',
          '10860', '06760', '09355', '07445', '14936', '03616', '10115', '14042',
          '15423', '09508', '10506', '06531', '05216', '10801', '07594', '15581',
          '14498', '01286', '09062', '14939', '13012', '16675', '03769', '07240',
          '04253', '15821', '06023', '03617', '05438', '16195', '13805', '08876',
          '06216', '11537', '02625', '12478', '06030', '09268', '14267', '09746',
          '05341', '09685', '12119', '14485', '11956', '15371', '15678', '07792',
          '16049', '05220', '05633', '04545', '07312', '13336', '09059', '07316',
          '06145', '16112', '05782', '05321', '15167', '13632', '13475', '12607',
          '04123', '15490', '08276', '03763', '13866', '02798', '15095', '06803',
          '10781', '13047', '12538', '06846', '12042', '15159', '15504', '04476',
          '13803', '10468', '11535', '07208', '09985', '13663', '06725', '15623',
          '10625', '12883', '08158', '15409', '06215', '17214', '07309', '12660',
          '05376', '08218', '03618', '03935', '07204', '03064', '05373', '07818',
          '16921', '04518', '01086', '07241', '06840', '15111', '11539', '09256',
          '06291', '16211', '16852', '14629', '06757', '06552', '07717', '16193',
          '04790', '05776', '16214', '16770', '06315', '04608', '16111', '03899',
          '14458', '06487', '06029', '08579', '06817', '05219', '13936', '03511',
          '01268', '08800', '06774', '14257', '14112', '01288', '09259', '09991',
          '09330', '15380', '02569', '14036', '06327', '04792', '15447', '17288',
          '12239', '17150', '07314', '15976', '06150', '14106', '05842', '06762',
          '05660', '13502', '13863', '14930', '05021', '09302', '12887', '13402',
          '07439', '04257', '16018', '10461', '06621', '16987', '11540', '02602',
          '05330', '07182', '07605', '06015', '16089', '05743', '08659', '06316',
          '04193', '13040', '13874', '07183', '14609', '07242', '15147', '04203',
          '06595', '03994', '16712', '05150', '12792', '12235', '05825', '09296',
          '15450', '17099', '03401', '15848', '11649', '04060', '04198', '12395',
          '12037', '08578', '11538', '14891', '04965', '07884', '05661', '05783',
          '10555', '11186', '16989', '06842', '15492', '04076', '15493', '09112',
          '01080', '05238', '11497', '14903', '05483', '15098', '05385', '05175',
          '02957', '17315', '14795', '16118', '13736', '08231', '01431', '07176',
          '02890', '07181', '06926', '10095', '02891', '16886', '01083', '13199', '03858']

# Parse and validate the command line
args = parser.parse_args()
get_ids = args.get_ids

queue_manager.MAX_ALLOWED_TIME = args.max_allowed_time
queue_manager.MAX_SUBPROCESS_CNT = args.max_subproc_cnt

LOG_DIR = HST_DIR['pipeline'] + '/logs'

logger = pdslogger.PdsLogger('pds.hst.run-pipeline')
if not args.quiet:
    logger.add_handler(pdslogger.stdout_handler)

# Define the log file
now = datetime.datetime.now().strftime('%Y-%m-%dT%H-%M-%S')
if args.log:
    if os.path.isdir(args.log):
        logpath = os.path.join(args.log, 'run-pipeline-' + now + '.log')
    else:
        parts = os.path.splitext(args.log)
        logpath = parts[0] + '-' + now + parts[1]
else:
    os.makedirs(LOG_DIR, exist_ok=True)
    logpath = LOG_DIR + '/run-pipeline-' + now + '.log'

logger.add_handler(pdslogger.file_handler(logpath))
LIMITS = {'info': -1, 'debug': -1, 'normal': -1}
logger.open('run-pipeline ' + ' '.join(sys.argv[1:]), limits=LIMITS)

# Query MAST to get the updated list of program ids with moving target flag set to True
# when --get-ids is in the command and there is no passed in proposal ids.
if get_ids and not args.proposal_ids:
    ids_li = query_hst_moving_targets(logger=logger)

# If some proposal ids are passed in, we ignore the default list of full ids. Pipeline
# will run on the passed in list of ids
proposal_ids = args.proposal_ids if args.proposal_ids else ids_li

run_pipeline(proposal_ids, logger)
# Clean up the staging directories
for id in proposal_ids:
    clean_up_staging_dir(id, logger)
logger.close()

##########################################################################################
