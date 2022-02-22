#sphinx-apidoc -f -F ../../truvari/ -o /data/docs/api/ --ext-autodoc -M 
# requires pip install  sphinx & sphinx_rtd_theme 
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# will overwrite your python
cd ../../
python3 setup.py install
cd -
sphinx-build $DIR $DIR/output $DIR/truvari.* 
