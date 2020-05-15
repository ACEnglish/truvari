
```
████████╗██████╗ ██╗   ██╗██╗   ██╗ █████╗ ██████╗ ██╗
╚══██╔══╝██╔══██╗██║   ██║██║   ██║██╔══██╗██╔══██╗██║
   ██║   ██████╔╝██║   ██║██║   ██║███████║██████╔╝██║
   ██║   ██╔══██╗██║   ██║╚██╗ ██╔╝██╔══██║██╔══██╗██║
   ██║   ██║  ██║╚██████╔╝ ╚████╔╝ ██║  ██║██║  ██║██║
   ╚═╝   ╚═╝  ╚═╝ ╚═════╝   ╚═══╝  ╚═╝  ╚═╝╚═╝  ╚═╝╚═╝
```

Structural variant toolkit for benchmarking, annotating and more for VCFs

[WIKI page](https://github.com/spiralgenetics/truvari/wiki) has detailed documentation.
See [Updates](https://github.com/spiralgenetics/truvari/wiki/Updates) on new versions.

## Installation
Truvari uses Python 3.7 and can be installed with pip:
```
  pip install Truvari 
```
[PyPi](https://pypi.org/project/Truvari/#history) has a history of all versions available. Pip installs all requirements EXCEPT external tools needed for running some annotations. See XYZ for details. 

To build and install Truvari from scratch:
```
  python setup.py sdist bdist_wheel
  pip install dist/Truvari-<version>.tar.gz
```
Where `<version>` is which version you just built.
 
 See [tags/](https://github.com/spiralgenetics/truvari/tags/) for a list of all available version.
 
## Quick Start

Each sub-command contains help documentation. Start with `truvari -h` to see available commands.

The current most common Truvari use case is for structural variation benchmarking:
```
  truvari bench -b base.vcf.gz -c comp.vcf.gz -r reference.fasta -o output_dir/
```
## Truvari Commands

 - bench - Performance metrics from comparison of two VCFs
 - stats - Basic SV relevant VCF stats
 - consistency - Consistency report between multiple VCFs
 - anno - Annotate a VCF

## More Information

Find more details and discussions about Truvari on the [WIKI page](https://github.com/spiralgenetics/truvari/wiki).

[https://www.spiralgenetics.com](https://www.spiralgenetics.com)
