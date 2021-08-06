"""
Utility to run pylint and update the badge
"""
import anybadge
from pylint.lint import Run
from pylint.reporters.text import TextReporter
from io import StringIO

fail_under = 9.33
output = StringIO()
Run(['--rcfile=.pylintrc', f'--fail-under={fail_under}', 'truvari'], reporter=TextReporter(output), exit=False)
output.seek(0)

# Get the score
pylint_score = float(output.read().strip().split('\n')[-1].split(' ')[-1].split('/')[0])


# Define thresholds: <2=red, <4=orange <8=yellow <10=green
thresholds = {2: 'red',
              4: 'orange',
              6: 'yellow',
              fail_under: 'green'}

badge = anybadge.Badge('pylint', output, thresholds=thresholds)

badge.write_badge('imgs/pylint.svg')

