[![PyPI version](https://badge.fury.io/py/Truvari.svg)](https://badge.fury.io/py/Truvari)
[![pylint](imgs/pylint.svg)](https://github.com/acenglish/truvari/actions/workflows/pylint.yml)
[![FuncTests](https://github.com/acenglish/truvari/actions/workflows/func_tests.yml/badge.svg?branch=develop&event=push)](https://github.com/acenglish/truvari/actions/workflows/func_tests.yml)
[![coverage](imgs/coverage.svg)](https://github.com/acenglish/truvari/actions/workflows/func_tests.yml)
[![develop](https://img.shields.io/github/commits-since/acenglish/truvari/v4.2.0)](https://github.com/ACEnglish/truvari/compare/v4.2.0...develop)
[![Downloads](https://static.pepy.tech/badge/truvari)](https://pepy.tech/project/truvari)

![Logo](https://raw.githubusercontent.com/ACEnglish/truvari/develop/imgs/BoxScale1_DarkBG.png)  
Toolkit for benchmarking, merging, and annotating Structural Variants

📚 [WIKI page](https://github.com/acenglish/truvari/wiki) has detailed documentation.  
📈 See [Updates](https://github.com/acenglish/truvari/wiki/Updates) on new versions.  
📝 Read our Papers ([#1](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02840-6), [#2](https://www.biorxiv.org/content/10.1101/2023.10.29.564632v1)) to learn more.

## 💻 Installation
Truvari uses Python 3.6+ and can be installed with pip:
```
  python3 -m pip install Truvari 
```
For details and more installation options, see [Installation](https://github.com/acenglish/truvari/wiki/Installation) on the wiki.

## ⏩ Quick Start

Each sub-command contains help documentation. Start with `truvari -h` to see available commands.

The current most common Truvari use case is for structural variation benchmarking:
```
  truvari bench -b base.vcf.gz -c comp.vcf.gz -o output_dir/
```
## 🧬 Truvari Commands

 - [bench](https://github.com/acenglish/truvari/wiki/bench) - Performance metrics from comparison of two VCFs
 - [collapse](https://github.com/acenglish/truvari/wiki/collapse) - Collapse possibly redundant VCF entries
 - [anno](https://github.com/acenglish/truvari/wiki/anno) - Add SV annotations to a VCF
 - [consistency](https://github.com/acenglish/truvari/wiki/consistency) - Consistency report between multiple VCFs
 - [vcf2df](https://github.com/acenglish/truvari/wiki/vcf2df) - Turn a VCF into a pandas DataFrame
 - [segment](https://github.com/acenglish/truvari/wiki/segment) - Normalization of SVs into disjointed genomic regions
 - [stratify](https://github.com/acenglish/truvari/wiki/stratify) - Count variants per-region in vcf
 - [divide](https://github.com/ACEnglish/truvari/wiki/divide) - Divide a VCF into independent shards
 - [phab](https://github.com/ACEnglish/truvari/wiki/phab) - Harmonize variant representations using MSA
 - [refine](https://github.com/ACEnglish/truvari/wiki/refine) - Automated bench result refinement with phab

## 🔎 More Information

All documentation about Truvari is on the [WIKI](https://github.com/acenglish/truvari/wiki). Additional information about using Truvari can be found in [Discussions](https://github.com/ACEnglish/truvari/discussions)
