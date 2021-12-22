"""
Utility to run pylint and update the badge
"""
import re
import sys
from io import StringIO
import anybadge
from pylint.lint import Run
from pylint.reporters.text import TextReporter

# run pylint
output = StringIO()
Run(['--rcfile=.pylintrc', 'truvari'], reporter=TextReporter(output), exit=False)
output.seek(0)
output = output.read()

# Get the score
search = re.search("Your code has been rated at (?P<score>.*)/10", output)
pylint_score = float(search.groupdict()['score'])
if pylint_score == 10:
    pylint_score = int(pylint_score)

# Define thresholds: <2=red, <4=orange <8=yellow <10=green
thresholds = {2: 'red',
              4: 'orange',
              6: 'yellow',
              10: 'green'}

badge = anybadge.Badge('pylint', pylint_score, thresholds=thresholds)
badge.write_badge('imgs/pylint.svg', overwrite=True)

sys.stdout.write(output)

# failunder
if pylint_score != 10:
    exit(1)
