################################################################################
# run_spt_tests.py
#
# This stand-alone application attempts to create a list of Target Idenfication
# tuples for every t the file SPT_TESTS.py. On success, it prints out the
# SPT filename followed by one or more target LIDs. On failure, it prints out
# the error message and continues. At the end, it prints a summary.
#
# Usage:
#   python3 run_spt_tests.py
#
################################################################################

from target_identifications import hst_target_identifications
from target_identifications.TESTS.SPT_TESTS import SPT_TESTS

errors_raised = 0
missing_targets = 0

for (filename,d) in SPT_TESTS:
    try:
        targets = hst_target_identifications(d, filename)
    except Exception as e:
        print('****', filename, 'ERROR', e)
        errors_raised += 1
    else:
        short_lids = [t[-1].rpartition('target:')[-1] for t in targets]
        if short_lids:
            print(filename, *tuple(short_lids))
        else:
            print(filename, 'NO TARGETS IDENTIFIED')
            missing_targets += 1

print()
print('Number of missing targets:', missing_targets)
print('Number of errors raised:', errors_raised)

if missing_targets or errors_raised:
    print('SPT test failed')
else:
    print('SPT test was successful')
