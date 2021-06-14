import codecs
import os.path
from setuptools import setup

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")


setup(
    name='Truvari',
    version=get_version('truvari/__init__.py'),
    author="ACEnglish",
    author_email="acenglish@gmail.com",
    url="https://github.com/spiralgenetics/truvari",
    packages=['truvari', 'truvari/annos'],
    license='MIT',
    scripts=["bin/truvari"],
    description="Structural variant comparison tool for VCFs",
    long_description=open('README.md', encoding='UTF-8').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        "ACEBinf>=1.0.2",
        "python-Levenshtein==0.12.2",
        "progressbar2>=3.41.0",
        "pysam>=0.15.2",
        "pyfaidx>=0.5.5.2",
        "intervaltree>=3.0.2",
        "joblib>=0.14.1",
        "numpy>=1.18.1",
        "pyfaidx>=0.5.8",
        "pytabix>=0.1",
        "bwapy>=0.1.4",
        "pandas>=1.0.5",
    ],
)
