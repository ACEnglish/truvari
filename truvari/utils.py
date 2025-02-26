"""
Miscellaneous utilites for truvari
"""
import os
import re
import sys
import gzip
import time
import signal
import logging
import argparse
import tempfile
import warnings
import subprocess
from datetime import timedelta
from collections import namedtuple
from importlib.metadata import version

import pysam
from pysam import bcftools

import truvari

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


def restricted_int(x):
    """
    Restrict int to positive. Raises argparse.ArgumentTypeError if int is negative
    Used with :class:`argparse.ArgumentParser.add_argument` type parameter

    :param `x`: number to check
    :type `x`: int

    :return: input number
    :rtype: float

    Example
        >>> import truvari
        >>> truvari.restricted_int(5)
        5
        >>> truvari.restricted_int(-2)
        Traceback (most recent call last):
        argparse.ArgumentTypeError: -2 is < 0
    """
    x = int(x)
    if x < 0:
        raise argparse.ArgumentTypeError(f"{x} is < 0")
    return x


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


def setup_logging(debug=False, stream=sys.stderr,
                  log_format="%(asctime)s [%(levelname)s] %(message)s",
                  show_version=False):
    """
    Create default logger

    :param `debug`: Set log-level to logging.DEBUG
    :type `debug`: boolean, optional
    :param `stream`: Where log is written
    :type `stream`: file handler, optional
    :param `log_format`: Format of log lines
    :type `log_format`: string, optional
    :param `show_version`: Start log with truvari version and command line args
    :type `show_version`: bool, optional
    """
    logLevel = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(stream=stream, level=logLevel, format=log_format)
    if show_version:
        logging.info(f"Truvari v{version('truvari')}")
        logging.info("Command %s", " ".join(sys.argv))

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


def cmd_exe(cmd, stdin=None, timeout=-1, cap_stderr=True, pipefail=False):
    """
    Executes a command through the shell and captures the output while enabling
    process handling with timeouts and keyboard interrupts.

    :param `cmd`: The command to be executed.
    :type `cmd`: string
    :param `stdin`: stdinput to pipe to command
    :type `stdin`: bytes
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
    m_in = sys.stdin if stdin is None else subprocess.PIPE
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                            stdin=m_in, stderr=stderr, close_fds=True,
                            start_new_session=True, executable="/bin/bash")
    if timeout > 0:
        signal.signal(signal.SIGALRM, alarm_handler)
        signal.alarm(int(timeout * 60))
    try:
        stdoutVal, stderrVal = proc.communicate(input=stdin)
        signal.alarm(0)  # reset the alarm
    except Alarm:
        logging.error(("Command was taking too long. "
                       "Automatic Timeout Initiated after %d"), timeout)
        os.killpg(proc.pid, signal.SIGTERM)
        proc.kill()
        return cmd_result(214, "", "", timedelta(seconds=time.time() - t_start))
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
    ret = cmd_result(retCode, stdoutVal, stderrVal.decode(),
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
        >>> gen = truvari.ref_ranges("repo_utils/test_files/references/reference.fa", 1000)
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
        >>> gen = truvari.bed_ranges("repo_utils/test_files/beds/giab.bed", 1000)
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


def vcf_ranges(vcf, min_dist=1000):
    """
    Chunk vcf into discrete pieces. Useful for multiprocessing.

    :param `vcf`: Filename of vcf to find ranges
    :type `vcf`: string
    :param `min_dist`: Minimum distance between entries for new range
    :type `min_dist`: int, optional

    :return: generator of tuples (ref_name, start, stop)
    :rtype: iterator

    Example
        >>> import truvari
        >>> gen = truvari.vcf_ranges("repo_utils/test_files/variants/input1.vcf.gz")
        >>> print(len([_ for _ in gen]))
        228
    """
    in_vcf = truvari.VariantFile(vcf)

    cur_chrom = None
    min_start = None
    max_end = None
    for entry in in_vcf:
        if cur_chrom is None:
            cur_chrom = entry.chrom
            min_start = entry.start
            max_end = entry.stop
        elif entry.chrom != cur_chrom:
            yield cur_chrom, min_start, max_end
            cur_chrom = entry.chrom
            min_start = entry.start
            max_end = entry.stop
        elif entry.start >= max_end + min_dist:
            yield cur_chrom, min_start, max_end
            cur_chrom = entry.chrom
            min_start = entry.start
            max_end = entry.stop
        else:
            max_end = max(max_end, entry.stop)

    yield cur_chrom, min_start, max_end


def opt_gz_open(in_fn):
    """
    Chooses file handler for plain-text files or `*.gz` files.

    returns a generator which yields lines of the file
    """
    def gz_hdlr(fn):
        with gzip.open(fn) as fh:
            for line in fh:
                yield line.decode()

    def fh_hdlr(fn):
        with open(fn) as fh:
            yield from fh

    if in_fn.endswith('.gz'):
        return gz_hdlr(in_fn)

    return fh_hdlr(in_fn)


def make_temp_filename(tmpdir=None, suffix=""):
    """
    Get a random filename in a tmpdir with an optional extension
    """
    if tmpdir is None:
        tmpdir = tempfile.gettempdir()
    # pylint: disable=protected-access
    fn = os.path.join(tmpdir, next(tempfile._get_candidate_names())) + suffix
    # pylint: enable=protected-access
    return fn


def help_unknown_cmd(user_cmd, avail_cmds, threshold=0.8):
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
        score = truvari.seqsim(user_cmd, real)
        if score >= threshold:
            guesses.append((score, real))
    guesses.sort()
    if not guesses:
        return None
    return guesses[-1][1]


def performance_metrics(tpbase, tp, fn, fp):
    """
    Calculates precision, recall, and f1 given counts by state
    Can return None if values are zero

    :param `tpbase`: found base count
    :type `tpbase`: int
    :param `tp`: found comp count
    :type `tp`: int
    :param `fn`: missed base count
    :type `fn`: int
    :param `fp`: missed comp count
    :type `fp`: int

    :return: tuple of precision, recall, f1
    :rtype: tuple, float
    """
    base_cnt = tpbase + fn
    comp_cnt = tp + fp
    if base_cnt == 0 or comp_cnt == 0:
        return None, None, None
    precision = tp / comp_cnt
    recall = tpbase / base_cnt
    neum = recall * precision
    denom = recall + precision
    f1 = 2 * (neum / denom) if denom != 0 else None
    return precision, recall, f1

def compress_index_vcf(fn, fout=None, remove=True):
    """
    Compress and index a vcf. If no fout is provided, write to fn + '.gz'

    :param `fn`: filename of vcf to work on
    :type `fn`: string
    :param `fout`: output filename
    :type `fout`: string, optional
    :param `remove`: remove fn after fout is made
    :type `remove`: bool
    """
    if fout is None:
        fout = fn + '.gz'
    with tempfile.TemporaryDirectory() as tmpdirname:
        bcftools.sort(fn, "-o", fout, "-O", "z", "--temp-dir", tmpdirname, catch_stdout=False)
    pysam.tabix_index(fout, force=True, preset="vcf")
    if remove:
        os.remove(fn)


def check_vcf_index(vcf_path):
    """
    Return true if an index file is found for the vcf
    """
    vcf_index_ext = ['tbi', 'csi']
    return any(os.path.exists(vcf_path + '.' + x) for x in vcf_index_ext)
