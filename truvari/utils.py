"""
Miscellaneous utilites for truvari
"""
import os
import sys
import time
import signal
import logging
import warnings
import subprocess
from datetime import timedelta
from collections import namedtuple

import progressbar

def setup_progressbar(size):
    """
    Return a formatted progress bar of size

    Parameters
    ----------
    size : int
        Number of elements in the progress bar

    Returns
    -------
    progressbar.ProgressBar
    """
    return progressbar.ProgressBar(redirect_stdout=True, max_value=size, widgets=[
        ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(size), '] ',
        progressbar.Bar(),
        ' (', progressbar.ETA(), ') ',
    ])


class LogFileStderr():
    """ Write to stderr and a file
    Useful in conjunction with `setup_logging(stream=LogFileStderr('log.txt'))`
    """

    def __init__(self, fn):
        """ keep these props """
        self.name = fn
        self.file_handler = open(fn, 'w') # pylint: disable=consider-using-with

    def write(self, *args):
        """ Write to both """
        sys.stderr.write(*args)
        self.file_handler.write(*args)

    def flush(self):
        """ Flush both """
        sys.stderr.flush()
        self.file_handler.flush()


def setup_logging(debug=False, stream=sys.stderr, log_format="%(asctime)s [%(levelname)s] %(message)s"):
    """
    Create default logger

    Parameters
    ----------
    debug : boolean, optional
        Set log-level to logging.DEBUG
    stream : file handler, optional (sys.stderr)
        Where log is written
    log_format : string, optional ("%(asctime)s [%(levelname)s] %(message)s")
        Format of log lines
    """
    logLevel = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(stream=stream, level=logLevel, format=log_format)
    logging.info("Running %s", " ".join(sys.argv))

    def sendWarningsToLog(message, category, filename, lineno):
        """
        Put warnings into logger
        """
        logging.warning('%s:%s: %s:%s', filename, lineno, category.__name__, message)

    warnings.showwarning = sendWarningsToLog


class Alarm(Exception):
    """ Alarm Class for command timeouts """
    pass # pylint: disable=unnecessary-pass


def alarm_handler(signum, frame=None):
    """ Alarm handler for command timeouts """
    raise Alarm


def cmd_exe(cmd, timeout=-1, cap_stderr=True, pipefail=False):
    """
    Executes a command through the shell and captures the output while enabling
    process handling with timeouts and keyboard interrupts.

    Parameters
    ----------
    cmd : string
        The command to be executed.
    timeout : int
        Timeout for the command in minutes. So 1440 means 24 hours. -1 means never
    cap_stderr : boolean
        If True, capture the stderr and return it as part of the returned cmd_results.
        Otherwise, stderr will be streamed through to the terminal
    pipefail : boolean
        Set to True if the cmd contains pipes `|`

    Returns
    -------
    namedtuple : ret_code, stdout, stderr, run_time
        ret_code : integer, exit code of the command
        stdout : string, captured standard output of the command
        stderr : binary string, captured standard error of the command
        run_time : datetime.timedelta, execution time

    Examples
    --------
        >>> import truvari
        >>> truvari.cmd_exe("ls")
        cmd_result(ret_code=0, stdout='anno_answers.vcf\nsegment.vcf\nx.py\n', stderr=b'',
         run_time=datetime.timedelta(microseconds=9452))
    """
    cmd_result = namedtuple("cmd_result", "ret_code stdout stderr run_time")
    t_start = time.time()
    stderr = subprocess.PIPE if cap_stderr else None
    if pipefail:
        cmd = f"set -o pipefail; {cmd}"
    # pylint: disable=consider-using-with
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stdin=sys.stdin, stderr=stderr, close_fds=True,
                            start_new_session=True, executable="/bin/bash")
    signal.signal(signal.SIGALRM, alarm_handler)
    if timeout > 0:
        signal.alarm(int(timeout * 60))
    try:
        stdoutVal, stderrVal = proc.communicate()
        signal.alarm(0)  # reset the alarm
    except Alarm:
        logging.error(("Command was taking too long. "
                       "Automatic Timeout Initiated after %d"), timeout)
        os.killpg(proc.pid, signal.SIGTERM)
        proc.kill()
        return cmd_result(214, None, None, timedelta(seconds=time.time() - t_start))
    except KeyboardInterrupt:
        logging.error("KeyboardInterrupt on cmd %s", cmd)
        os.killpg(proc.pid, signal.SIGKILL)
        proc.kill()
        try:
            sys.exit(214)
        except SystemExit:
            os._exit(214) # pylint: disable=protected-access

    stdoutVal = bytes.decode(stdoutVal)
    retCode = proc.returncode
    ret = cmd_result(retCode, stdoutVal, stderrVal, timedelta(seconds=time.time() - t_start))
    return ret
