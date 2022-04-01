"""
Miscellaneous utilites for truvari
"""
import os
import re
import sys
import time
import signal
import logging
import argparse
import warnings
import subprocess
from datetime import timedelta
from collections import namedtuple

import Levenshtein
import pysam
import progressbar

HEADERMAT = re.compile(
    r"##\w+=<ID=(?P<name>\w+),Number=(?P<num>[\.01AGR]),Type=(?P<type>\w+)")

def restricted_float(x):
    """
    Restrict float to range (0,1). Raises argparse.ArgumentTypeError if float is out of range
    Used with :class:`argparse.ArgumentParser.add_argument` type parameter

    :param `x`: number to check
    :type `x`: float

    :return: input number
    :rtype: float

    Example
        >>> import truvari
        >>> truvari.restricted_float(0.5)
        0.5
        >>> truvari.restricted_float(2)
        Traceback (most recent call last):
        argparse.ArgumentTypeError: 2.0 not in range (0, 1)
    """
    x = float(x)
    if x < 0 or x > 1:
        raise argparse.ArgumentTypeError(
            f"{x} not in range (0, 1)")
    return x


def setup_progressbar(size):
    """
    Build a formatted :class:`progressbar.ProgressBar`

    :param `size`: Number of elements in the progress bar
    :type `size`: int

    :return: Formatted progress bar
    :rtype: progressbar.ProgressBar

    Example
        >>> import truvari
        >>> import time
        >>> bar = truvari.setup_progressbar(4)
        >>> for i in range(4):
        ...     bar.update(i + 1) # The bar animation updates
        ...     time.sleep(1)
        >>> bar.finish() # Final update to the bar
    """
    return progressbar.ProgressBar(redirect_stdout=True, max_value=size, widgets=[
        ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(size), '] ',
        progressbar.Bar(),
        ' (', progressbar.ETA(), ') ',
    ])


class LogFileStderr():
    """ Write to stderr and a file

    Useful in conjunction with :meth:`setup_logging`

    Example
        >>> import truvari
        >>> truvari.setup_logging(stream=LogFileStderr('log.txt'))
    """

    def __init__(self, fn):
        """ keep these props """
        self.name = fn
        # pylint: disable=consider-using-with
        self.file_handler = open(fn, 'w')

    def write(self, *args):
        """ Write to handlers """
        sys.stderr.write(*args)
        self.file_handler.write(*args)

    def flush(self):
        """ Flush handlers """
        sys.stderr.flush()
        self.file_handler.flush()


def setup_logging(debug=False, stream=sys.stderr, log_format="%(asctime)s [%(levelname)s] %(message)s"):
    """
    Create default logger

    :param `debug`: Set log-level to logging.DEBUG
    :type `debug`: boolean, optional
    :param `stream`: Where log is written
    :type `stream`: file handler, optional
    :param `log_format`: Format of log lines
    :type `log_format`: string, optional
    """
    logLevel = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(stream=stream, level=logLevel, format=log_format)
    logging.info("Running %s", " ".join(sys.argv))

    def sendWarningsToLog(message, category, filename, lineno, *args, **kwargs):  # pylint: disable=unused-argument
        """
        Put warnings into logger
        """
        logging.warning('%s:%s: %s:%s', filename, lineno,
                        category.__name__, message)

    warnings.showwarning = sendWarningsToLog


class Alarm(Exception):
    """ Alarm Class for command timeouts """
    pass  # pylint: disable=unnecessary-pass


def alarm_handler(signum, frame=None):
    """ Alarm handler for command timeouts """
    raise Alarm


def cmd_exe(cmd, timeout=-1, cap_stderr=True, pipefail=False):
    """
    Executes a command through the shell and captures the output while enabling
    process handling with timeouts and keyboard interrupts.

    :param `cmd`: The command to be executed.
    :type `cmd`: string
    :param `timeout`: Timeout for the command in minutes. So 1440 means 24 hours. -1 means never
    :type `timeout`: int
    :param `cap_stderr`: If True, capture the stderr and return it as part of the returned cmd_results.
        Otherwise, stderr will be streamed through to the terminal
    :type `cap_stderr`: boolean
    :param `pipefail`: Set to True if the cmd contains pipes `|`
    :type `pipefail`: boolean

    :return: | namedtuple of
             | ret_code - integer, exit code of the command
             | stdout - string, captured standard output of the command
             | stderr - binary string, captured standard error of the command
             | run_time - datetime.timedelta, execution time
    :rtype: namedtuple (ret_code, stdout, stderr, run_time)

    Example
        >>> import truvari
        >>> ret = truvari.cmd_exe("echo 'hello world'")
        >>> ret.ret_code
        0
        >>> ret.stdout
        'hello world\\n'
        >>> import truvari
        >>> ret = truvari.cmd_exe("sleep 5", pipefail=True, timeout=2/60) # Error logged is caught
        >>> ret.ret_code
        214
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
            os._exit(214)  # pylint: disable=protected-access

    stdoutVal = bytes.decode(stdoutVal)
    retCode = proc.returncode
    ret = cmd_result(retCode, stdoutVal, stderrVal,
                     timedelta(seconds=time.time() - t_start))
    return ret


def ref_ranges(reference, chunk_size=10000000):
    """
    Chunk reference into pieces. Useful for multiprocessing.

    :param `reference`: Filename of reference to chunk
    :type `reference`: string
    :param `chunk_size`: Length of reference chunks
    :type `chunk_size`: int, optional

    :return: generator of tuples (ref_name, start, stop)
    :rtype: iterator

    Example
        >>> import truvari
        >>> gen = truvari.ref_ranges("repo_utils/test_files/reference.fa", 1000)
        >>> print(len([_ for _ in gen]))
        1000
    """
    ref = pysam.FastaFile(reference)
    for ref_name in ref.references:  # pylint: disable=not-an-iterable
        final_stop = ref.get_reference_length(ref_name)
        start = 0
        stop = start + chunk_size
        while stop < final_stop:
            yield ref_name, start, stop
            start = stop
            stop += chunk_size
        yield ref_name, start, final_stop


def bed_ranges(bed, chunk_size=10000000):
    """
    Chunk bed regions into pieces. Useful for multiprocessing.

    :param `bed`: Filename of bed to chunk
    :type `bed`: string
    :param `chunk_size`: Length of reference chunks
    :type `chunk_size`: int, optional

    :return: generator of tuples (ref_name, start, stop)
    :rtype: iterator

    Example
        >>> import truvari
        >>> gen = truvari.bed_ranges("repo_utils/test_files/giab.bed", 1000)
        >>> print(len([_ for _ in gen]))
        992
    """
    with open(bed, 'r') as fh:
        for line in fh:
            data = line.strip().split('\t')
            start = int(data[1])
            final_stop = int(data[2])
            stop = start + chunk_size
            while stop < final_stop:
                yield data[0], start, stop
                start = stop
                stop += chunk_size
            yield data[0], start, final_stop


def help_unknown_cmd(user_cmd, avail_cmds, threshold=.5):
    """
    Guess the command in avail_cmds that's most similar to user_cmd. If there is
    no guess below threshold, don't guess.

    :param `user_cmd`: user command
    :type `user_cmd`: string
    :param `avail_cmds`: available commands
    :type `avail_cmds`: list of strings
    :param `threshold`: max score (0-1) to give a guess
    :type `threshold`: float, optional

    :return: best guess from :param: `avail_cmds` or None
    :rtype: string

    Example
        >>> import truvari
        >>> truvari.help_unknown_cmd('banno', ["bench", "anno"])
        'anno'
        >>> truvari.help_unknown_cmd('random', ["bench", "anno"])
    """
    guesses = []
    for real in avail_cmds:
        ratio = Levenshtein.ratio(user_cmd, real)
        distance = Levenshtein.distance(user_cmd, real)
        score = distance - ratio
        score2 = 1 - distance / len(real)
        dat = (score, real)
        if ratio >= threshold and score2 >= threshold:
            guesses.append(dat)
    guesses.sort()
    if not guesses:
        return None
    return guesses[0][1]
