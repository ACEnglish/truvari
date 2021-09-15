```
████████╗██████╗ ██╗   ██╗██╗   ██╗ █████╗ ██████╗ ██╗
╚══██╔══╝██╔══██╗██║   ██║██║   ██║██╔══██╗██╔══██╗██║
   ██║   ██████╔╝██║   ██║██║   ██║███████║██████╔╝██║
   ██║   ██╔══██╗██║   ██║╚██╗ ██╔╝██╔══██║██╔══██╗██║
   ██║   ██║  ██║╚██████╔╝ ╚████╔╝ ██║  ██║██║  ██║██║
   ╚═╝   ╚═╝  ╚═╝ ╚═════╝   ╚═══╝  ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝
```
[![PyPI version](https://badge.fury.io/py/Truvari.svg)](https://badge.fury.io/py/Truvari)
[![pylint](imgs/pylint.svg)](https://github.com/spiralgenetics/truvari/actions/workflows/pylint.yml)
[![FuncTests](https://github.com/spiralgenetics/truvari/actions/workflows/func_tests.yml/badge.svg?branch=develop&event=push)](https://github.com/spiralgenetics/truvari/actions/workflows/func_tests.yml)
[![coverage](imgs/coverage.svg)](https://github.com/spiralgenetics/truvari/actions/workflows/func_tests.yml)
[![develop](https://img.shields.io/github/commits-since/spiralgenetics/truvari/v3.0.0)](https://github.com/spiralgenetics/truvari/commits/develop)

Structural variant toolkit for benchmarking, annotating and more for VCFs

[WIKI page](https://github.com/spiralgenetics/truvari/wiki) has detailed documentation.
See [Updates](https://github.com/spiralgenetics/truvari/wiki/Updates) on new versions.

## Installation
Truvari uses Python 3.6+ and can be installed with pip:
```
  python3 -m pip install Truvari 
```
[PyPi](https://pypi.org/project/Truvari/#history) has a history of all versions available. Pip installs all requirements EXCEPT external tools needed for running some annotations. See [anno](https://github.com/spiralgenetics/truvari/wiki/anno) for details. 

To build and install Truvari from scratch:
```
  python3 setup.py install
```
 
See [tags/](https://github.com/spiralgenetics/truvari/tags/) for a list of all available versions.

See [Development/docker](https://github.com/spiralgenetics/truvari/wiki/Development#docker) for instructions on how to
build and use a Truvari Docker container.

## Quick Start

Each sub-command contains help documentation. Start with `truvari -h` to see available commands.

The current most common Truvari use case is for structural variation benchmarking:
```
  truvari bench -b base.vcf.gz -c comp.vcf.gz -f reference.fasta -o output_dir/
```
## Truvari Commands

 - [bench](https://github.com/spiralgenetics/truvari/wiki/bench) - Performance metrics from comparison of two VCFs
 - [collapse](https://github.com/spiralgenetics/truvari/wiki/collapse) - Collapse possibly redundant VCF entries
 - [anno](https://github.com/spiralgenetics/truvari/wiki/anno) - Add SV annotations to a VCF
 - [vcf2df](https://github.com/spiralgenetics/truvari/wiki/vcf2df) - Turn a VCF into a pandas DataFrame
 - [consistency](https://github.com/spiralgenetics/truvari/wiki/consistency) - Consistency report between multiple VCFs

## More Information

Find more details and discussions about Truvari on the [WIKI page](https://github.com/spiralgenetics/truvari/wiki).

[https://www.spiralgenetics.com](https://www.spiralgenetics.com)
