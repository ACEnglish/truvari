DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
python3 -m pip install sphinx sphinx_rtd_theme
# will overwrite your current environment's truvari
cd ../../
python3 setup.py install
cd -
sphinx-build $DIR $DIR/output $DIR/truvari.* 
