```
████████╗██████╗ ██╗   ██╗██╗   ██╗ █████╗ ██████╗ ██╗
╚══██╔══╝██╔══██╗██║   ██║██║   ██║██╔══██╗██╔══██╗██║
   ██║   ██████╔╝██║   ██║██║   ██║███████║██████╔╝██║
   ██║   ██╔══██╗██║   ██║╚██╗ ██╔╝██╔══██║██╔══██╗██║
   ██║   ██║  ██║╚██████╔╝ ╚████╔╝ ██║  ██║██║  ██║██║
   ╚═╝   ╚═╝  ╚═╝ ╚═════╝   ╚═══╝  ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝
```
[![PyPI version](https://badge.fury.io/py/Truvari.svg)](https://badge.fury.io/py/Truvari)
[![pylint](imgs/pylint.svg)](https://github.com/acenglish/truvari/actions/workflows/pylint.yml)
[![FuncTests](https://github.com/acenglish/truvari/actions/workflows/func_tests.yml/badge.svg?branch=develop&event=push)](https://github.com/acenglish/truvari/actions/workflows/func_tests.yml)
[![coverage](imgs/coverage.svg)](https://github.com/acenglish/truvari/actions/workflows/func_tests.yml)
[![develop](https://img.shields.io/github/commits-since/acenglish/truvari/v3.4.0)](https://github.com/ACEnglish/truvari/compare/v3.4.0...develop)
[![Downloads](https://pepy.tech/badge/truvari)](https://pepy.tech/project/truvari)

Toolkit for benchmarking, merging, and annotating Structrual Variants

[WIKI page](https://github.com/acenglish/truvari/wiki) has detailed documentation.  
See [Updates](https://github.com/acenglish/truvari/wiki/Updates) on new versions.  
Read our [Paper](https://doi.org/10.1101/2022.02.21.481353) for more details.

## Installation
Truvari uses Python 3.6+ and can be installed with pip:
```
  python3 -m pip install Truvari 
```
For details and more installation options, see [Installation](https://github.com/acenglish/truvari/wiki/Installation) on the wiki.

## Quick Start

Each sub-command contains help documentation. Start with `truvari -h` to see available commands.

The current most common Truvari use case is for structural variation benchmarking:
```
  truvari bench -b base.vcf.gz -c comp.vcf.gz -f reference.fasta -o output_dir/
```
## Truvari Commands

 - [bench](https://github.com/acenglish/truvari/wiki/bench) - Performance metrics from comparison of two VCFs
 - [collapse](https://github.com/acenglish/truvari/wiki/collapse) - Collapse possibly redundant VCF entries
 - [anno](https://github.com/acenglish/truvari/wiki/anno) - Add SV annotations to a VCF
 - [vcf2df](https://github.com/acenglish/truvari/wiki/vcf2df) - Turn a VCF into a pandas DataFrame
 - [consistency](https://github.com/acenglish/truvari/wiki/consistency) - Consistency report between multiple VCFs
 - [divide](https://github.com/ACEnglish/truvari/wiki/divide) - Divide a VCF into independent parts
 - [segment](https://github.com/acenglish/truvari/wiki/segment) - Normalization of SVs into disjointed genomic regions

## More Information

Find more details and discussions about Truvari on the [WIKI page](https://github.com/acenglish/truvari/wiki).
