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
from collections import OrderedDict, namedtuple

import progressbar

class Alarm(Exception):

    """ Alarm Class for command timeouts """
    pass # pylint: disable=unnecessary-pass

def alarm_handler(signum, frame=None):
    """ Alarm handler for command timeouts """
    raise Alarm

class StatsBox(OrderedDict):
    """
    Make a blank stats box for counting TP/FP/FN and calculating performance
    """

    def __init__(self):
        super().__init__()
        self["TP-base"] = 0
        self["TP-call"] = 0
        self["FP"] = 0
        self["FN"] = 0
        self["precision"] = 0
        self["recall"] = 0
        self["f1"] = 0
        self["base cnt"] = 0
        self["call cnt"] = 0
        self["TP-call_TP-gt"] = 0
        self["TP-call_FP-gt"] = 0
        self["TP-base_TP-gt"] = 0
        self["TP-base_FP-gt"] = 0
        self["gt_concordance"] = 0

    def calc_performance(self, peek=False):
        """
        Calculate the precision/recall
        """
        do_stats_math = True
        if self["TP-base"] == 0 and self["FN"] == 0:
            logging.warning("No TP or FN calls in base!")
            do_stats_math = False
        elif self["TP-call"] == 0 and self["FP"] == 0:
            logging.warning("No TP or FP calls in comp!")
            do_stats_math = False
        elif peek:
            logging.info("Results peek: %d TP-base %d FN %.2f%% Recall", self["TP-base"], self["FN"],
                         100 * (float(self["TP-base"]) / (self["TP-base"] + self["FN"])))
        if peek:
            return

        # Final calculations
        if do_stats_math:
            self["precision"] = float(self["TP-call"]) / (self["TP-call"] + self["FP"])
            self["recall"] = float(self["TP-base"]) / (self["TP-base"] + self["FN"])
            if self["TP-call_TP-gt"] + self["TP-call_FP-gt"] != 0:
                self["gt_concordance"] = float(self["TP-call_TP-gt"]) / (self["TP-call_TP-gt"] +
                                                                         self["TP-call_FP-gt"])

        # f-measure
        neum = self["recall"] * self["precision"]
        denom = self["recall"] + self["precision"]
        if denom != 0:
            self["f1"] = 2 * (neum / denom)
        else:
            self["f1"] = "NaN"

def setup_progressbar(size):
    """
    Return a formatted progress bar of size
    """
    return progressbar.ProgressBar(redirect_stdout=True, max_value=size, widgets=[
        ' [', progressbar.Timer(), ' ', progressbar.Counter(), '/', str(size), '] ',
        progressbar.Bar(),
        ' (', progressbar.ETA(), ') ',
    ])


class LogFileStderr():

    """ Write to stderr and a file"""

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


def setup_logging(debug=False, stream=sys.stderr, log_format=None):
    """
    Default logger
    """
    logLevel = logging.DEBUG if debug else logging.INFO
    if log_format is None:
        log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(stream=stream, level=logLevel, format=log_format)
    logging.info("Running %s", " ".join(sys.argv))

    def sendWarningsToLog(message, category, filename, lineno):
        """
        Put warnings into logger
        """
        logging.warning('%s:%s: %s:%s', filename, lineno, category.__name__, message)

    # pylint: disable=unused-variable
    old_showwarning = warnings.showwarning
    warnings.showwarning = sendWarningsToLog

def cmd_exe(cmd, timeout=-1, cap_stderr=True, pipefail=False):
    """
    Executes a command through the shell.
    timeout in minutes! so 1440 mean is 24 hours.
    -1 means never
    returns namedtuple(ret_code, stdout, stderr, run_time)
    where ret_code is the exit code for the command executed
    stdout/err is the Standard Output Error from the command
    and runtime is hh:mm:ss of the execution time
    cap_stderr will capture the stderr and return it as part of the
    returned cmd_result. Otherwise, stderr will be streamed through.
    set pipefail=True if you're using pipes
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
