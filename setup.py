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
    packages=['truvari', 'truvari/annotations'],
    license='MIT',
    description="Structural variant comparison tool for VCFs",
    long_description=open('README.md', encoding='UTF-8').read(),
    long_description_content_type='text/markdown',
    entry_points={
      'console_scripts': [
         'truvari = truvari.__main__:main'
      ]
    },
    install_requires=[
        "rich==12.5.1",
        "python-Levenshtein==0.12.2",
        "edlib>=1.3.8.post2",
        "progressbar2>=3.41.0",
        "pysam>=0.15.2",
        "intervaltree>=3.0.2",
        "joblib>=1.0.1",
        "numpy>=1.21.2",
        "pytabix>=0.1",
        "bwapy>=0.1.4",
        "pandas>=1.3.3"
    ],
)
