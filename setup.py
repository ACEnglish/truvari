from setuptools import setup

setup(
    name='Truvari',
    version='1.1',
    author="ACEnglish",
    author_email="acenglish@gmail.com",
    url="https://github.com/spiralgenetics/truvari",
    packages=['truvari',],
    license='MIT',
    scripts=["truvari/truvari", "truvari/consistency_report"],
    description="Structural variant comparison tool for VCFs",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        "PyVCF==0.6.8",
        "python-Levenshtein==0.12.0",
        "progressbar2==3.41.0",
        "pysam==0.15.2",
        "pyfaidx==0.5.5.2",
        "intervaltree==3.0.2",
    ],
)
