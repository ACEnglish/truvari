"""
Utility to run pylint and update the badge
"""
import re
import sys
import json
from io import StringIO
import anybadge

# Get the score
data = json.load(open(sys.argv[1], 'r'))
coverage_pct = int(data['totals']['percent_covered'])

# Define thresholds: <2=red, <4=orange <8=yellow <10=green
thresholds = {20: 'red',
              40: 'orange',
              70: 'yellow',
              90: 'green'}

badge = anybadge.Badge('coverage', coverage_pct,
                       thresholds=thresholds, value_suffix="%")
badge.write_badge('imgs/coverage.svg', overwrite=True)

# failunder
if coverage_pct < 85:
    exit(1)
