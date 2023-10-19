##########################################################################################
# queue_manager/__init__.py
#
# Queue manager module that will queue in the next task for the hst pipeline process.
# - run_pipeline will start a hst pipeline process for the given proposal id list.
# - queue_next_task will Queue in the next task for a given proposal id to database, and
#   wait for the open subprocess slot to execute the corresponding command.
##########################################################################################

import os
import pdslogger
import psutil
import subprocess
import time

from hst_helper.fs_utils import get_formatted_proposal_id
from queue_manager.task_queue_db import (add_a_prog_id_task_queue,
                                         add_a_subprocess,
                                         create_task_queue_table,
                                         db_exists,
                                         erase_all_task_queue,
                                         get_all_subprocess_info,
                                         get_next_task_to_be_run,
                                         get_pid_by_prog_id_task_and_visit,
                                         get_total_number_of_subprocesses,
                                         init_task_queue_table,
                                         remove_a_subprocess_by_pid,
                                         update_a_prog_id_task_status)
from queue_manager.config import (DB_PATH,
                                  PYTHON_EXE,
                                  HST_SOURCE_ROOT,
                                  MAX_ALLOWED_TIME,
                                  MAX_SUBPROCESS_CNT,
                                  TASK_NUM_TO_CMD_MAPPING,
                                  TASK_NUM_TO_PREV_TASK_MAPPING,
                                  TASK_NUM_TO_PRI_MAPPING)
from sqlalchemy.exc import OperationalError

def run_pipeline(proposal_ids, logger=None):
    """With a given list of proposal ids, run pipeline for each program id.

    Inputs:
        proposal_ids    a list of proposal ids.
        logger          pdslogger to use; None for default EasyLogger.
    """
    logger = logger or pdslogger.EasyLogger()
    logger.info(f'Run pipeline with proposal ids: {proposal_ids}')

    try:
        init_task_queue_table()
    except OperationalError as e:
        if 'already exists' in e.__repr__():
            erase_all_task_queue()
        elif 'no such table' in e.__repr__():
            create_task_queue_table()
        else:
            logger.error('Failed to create task queue table!')
            raise Exception('Failed to create task queue table!') # fatal error

    # Kick start the pipeline for each proposal id
    # proc_li = []
    for prog_id in proposal_ids:
        try:
            proposal_id = int(prog_id)
        except ValueError:
            logger.warn(f'Proposal id: {prog_id} is not valid.')

        formatted_proposal_id = get_formatted_proposal_id(proposal_id)
        # Start hst pipeline for each proposal id
        logger.info(f'Starting to run pipeline for {proposal_id}')
        logger.info(f'Queue query_hst_moving_targets for {proposal_id}')
        queue_next_task(formatted_proposal_id, '', 0, logger)

def queue_next_task(proposal_id, visit_info, task_num, logger):
    """Queue in the next task for a given proposal id to database, and wait for the open
    subprocess slot to execute the corresponding command. Once there is an open slot,
    update the task queue status

    1. Update the next task for a given proposal id to database. (task status: 0)
    2. Wait for an open slot in subprocess list, once there is an open slot, update
    the task status to 1 in the database.
    3. Run the task command (spawn the subprocess)

    Inputs:
        proposal_id    the proposal if of the current task.
        visit_info     a two character visit, a list of visits or ''.
        task_num       a number represents the current task.
        logger         pdslogger to use; None for default EasyLogger.

    Returns:    the child process that executes the given task.
    """
    # if DB doesn't exist, log a warning message and return
    if not db_exists():
        logger.warn(f'Task queue db: {DB_PATH} does not exist.')
        return

    formatted_proposal_id = get_formatted_proposal_id(proposal_id)
    logger = logger or pdslogger.EasyLogger()
    logger.info(f'Queue in the next task for: {formatted_proposal_id}'
                f', task num: {task_num}')

    visit_arg = ' '.join(visit_info) if isinstance(visit_info, list) else visit_info
    visit = '' if isinstance(visit_info, list) else visit_info

    priority = TASK_NUM_TO_PRI_MAPPING[task_num]
    cmd = TASK_NUM_TO_CMD_MAPPING[task_num].replace('{P}', formatted_proposal_id)
    cmd = cmd.replace('{V}', visit_arg)
    # if the task has been queued, we don't spawn duplicated subprocess.
    spawn_subproc = add_a_prog_id_task_queue(formatted_proposal_id, visit,
                                             task_num, priority, 0, cmd)
    if spawn_subproc is False:
        return

    cmd_parts = cmd.split(' ')
    program_path = os.path.join(HST_SOURCE_ROOT, cmd_parts[0])
    args = [PYTHON_EXE, program_path] + cmd_parts[1:]
    max_allowed_time = MAX_ALLOWED_TIME
    # spawn the task subprocess
    pid = run_and_maybe_wait(task_num, args, max_allowed_time,
                             formatted_proposal_id, visit, logger)

    return pid

def run_and_maybe_wait(task_num, args, max_allowed_time, proposal_id, visit, logger):
    """Run one subprocess, waiting as necessary for a slot to open up.

    Inputs:
        task_num            a number represents the current task.
        args                the command of a subprocess to be executed.
        max_allowed_time    the max time for a subprocess to be done before killing it.
        proposal_id         the proposal id of the current task.
        visit               two character visit.
        logger              pdslogger to use; None for default EasyLogger.

    Returns:    the child process that executes the given args.
    """
    # query database and see if there is a higher priority job waiting (status: 0) to
    # be run, if so, spawn that subprocess first.
    task = get_next_task_to_be_run()
    while (task is not None
           and task.proposal_id != proposal_id
           and task.visit != visit):
        cmd_parts = task.cmd.split(' ')
        program_path = os.path.join(HST_SOURCE_ROOT, cmd_parts[0])
        sub_args = [PYTHON_EXE, program_path] + cmd_parts[1::]
        run_and_maybe_wait(task.task_num, sub_args,  max_allowed_time,
                           task.proposal_id, task.visit, logger)
        # logger.debug("Spawning subprocess %s", str(sub_args))
        task = get_next_task_to_be_run()

    # wait for an open subprocess slot
    wait_for_subprocess(task_num)

    update_a_prog_id_task_status(proposal_id, visit, 1)
    logger.debug("Spawning subprocess %s", str(args))
    pid = subprocess.Popen(args)
    add_a_subprocess(pid.pid, task_num, time.time(),
                     time.time()+max_allowed_time, visit, proposal_id)

    return pid

def wait_for_subprocess(task_num, all=False):
    """Wait for one (or all) subprocess slots to open up.

    Inputs:
        task_num    a number represents the current task.
        all         a flag to determine if we are waiting for all subprocess slots to
                    open up.
    """
    subprocess_count = MAX_SUBPROCESS_CNT

    if all:
       subprocess_count = 0

    while get_total_number_of_subprocesses() > 0:
        subproc_info = get_all_subprocess_info()
        for info in subproc_info:
            pid, _, _, proc_max_time, visit, prog_id = info
            prev_task = TASK_NUM_TO_PREV_TASK_MAPPING[task_num]
            prev_pid_li = get_pid_by_prog_id_task_and_visit(prog_id, prev_task, visit)
            # If the subprocess is to be executed at the end of previous task, open the
            # slot occupied by previous task.
            if len(prev_pid_li) > 0:
                for prev_pid in prev_pid_li:
                    remove_a_subprocess_by_pid(prev_pid)
                break

            if not psutil.pid_exists(pid):
                # The subprocess completed, make the slot available for next subprocess
                remove_a_subprocess_by_pid(pid)
                break
            elif psutil.Process(pid).status == psutil.STATUS_ZOMBIE:
                # The subprocess is not running, kill the subprocess and open the slot
                psutil.Process(pid).terminate()
                remove_a_subprocess_by_pid(pid)
                break

            cur_time = time.time()
            if cur_time > proc_max_time and pid:
                # If a subprocess has been running for too long, kill it
                try:
                    psutil.Process(pid).terminate()
                except ProcessLookupError:
                    pass
                remove_a_subprocess_by_pid(pid)
                break

        if get_total_number_of_subprocesses() <= subprocess_count:
            # A slot opened up! Or all processes finished. Depending on what we're
            # waiting for.
            break
        time.sleep(1)
