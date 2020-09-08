from setuptools import setup

setup(
    name='Truvari',
    version='2.0.2',
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
        "python-Levenshtein>=0.12.0",
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
