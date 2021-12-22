#sphinx-apidoc -f -F ../../truvari/ -o /data/docs/api/ --ext-autodoc -M 
# requires pip install  sphinx & sphinx_rtd_theme 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# will overwrite your python
python3 $DIR/../../setup.py install
sphinx-build $DIR $DIR/output $DIR/truvari.* 
