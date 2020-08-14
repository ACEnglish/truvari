"""
Miscellaneous utilites for truvari
"""
import sys
import logging
import warnings
from collections import OrderedDict
import progressbar


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
        self["gt_precision"] = 0
        self["gt_recall"] = 0
        self["gt_f1"] = 0

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
            if self["TP-call_TP-gt"] != 0:
                self["gt_precision"] = float(self["TP-call_TP-gt"]) / (self["TP-call_TP-gt"] +
                                                                       self["FP"] + self["TP-call_FP-gt"])
                self["gt_recall"] = float(self["TP-base_TP-gt"]) / (self["TP-base_TP-gt"] + self["FN"])

        # f-measure
        neum = self["recall"] * self["precision"]
        denom = self["recall"] + self["precision"]
        if denom != 0:
            self["f1"] = 2 * (neum / denom)
        else:
            self["f1"] = "NaN"

        neum = self["gt_recall"] * self["gt_precision"]
        denom = self["gt_recall"] + self["gt_precision"]
        if denom != 0:
            self["gt_f1"] = 2 * (neum / denom)
        else:
            self["gt_f1"] = "NaN"


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
        self.file_handler = open(fn, 'w')

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
