"""
Utility to run pylint and update the badge
"""
import re
import sys
from io import StringIO
import anybadge
from pylint.lint import Run
from pylint.reporters.text import TextReporter

fail_under = 9.33

# run pylint
output = StringIO()
Run(['--rcfile=.pylintrc', f'--fail-under={fail_under}', 'truvari'], reporter=TextReporter(output), exit=False)
output.seek(0)
output = output.read()

# Get the score
search = re.search("Your code has been rated at (?P<score>.*)/10", output)
pylint_score = search.groupdict()['score']

# Define thresholds: <2=red, <4=orange <8=yellow <10=green
thresholds = {2: 'red',
              4: 'orange',
              6: 'yellow',
              fail_under: 'green'}

badge = anybadge.Badge('pylint', pylint_score, thresholds=thresholds)

badge.write_badge('imgs/pylint.svg', overwrite=True)

sys.stdout.write(output)

