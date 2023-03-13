Recommended
===========
For stable versions of Truvari, use pip
```
python3 -m pip install truvari
```
Specific versions can be installed via
```
python3 -m pip install truvari==3.2.0
```
See [pypi](https://pypi.org/project/Truvari/#history) for a history of all distributed releases.

Manual Installation
===================
To build Truvari directly, clone the repository and switch to a specific tag.
```
git clone https://github.com/spiralgenetics/truvari.git
git checkout tags/v3.0.0
python3 -m pip install .
```

To see a list of all available tags, run:
```
git tag -l
```

If you have an older clone of the repository and don't see the version you're looking for in tags, make sure to pull the latest changes:
```
git pull
git fetch --all --tags
```

Mamba / Conda
=============
NOTE!! There is a very old version of Truvari on bioconda that - for unknown reasons - supersedes the newer, supported versions. Users may need to specify to conda which release to build. See [this ticket](https://github.com/ACEnglish/truvari/issues/130#issuecomment-1196607866) for details.

Truvari releases are automatically deployed to bioconda. 
Users can follow instructions here (https://mamba.readthedocs.io/en/latest/installation.html) to install mamba. (A faster alternative conda compatible package manager.)

Creating an environment with Truvari and its dependencies.
```
mamba create -c conda-forge -c bioconda -n truvari truvari
```

Alternatively, see the [conda page](https://anaconda.org/bioconda/truvari) for details
```
conda install -c bioconda truvari
```

Building from develop
=====================
The default branch is `develop`, which holds in-development changes. This is for developers or those wishing to try experimental features and is not recommended for production. Development is versioned higher than the most recent stable release with an added suffix (e.g. Current stable release is `3.0.0`, develop holds `3.1.0-dev`). If you'd like to install develop, repeat the steps above but without `git checkout tags/v3.0.0`. See [wiki](https://github.com/spiralgenetics/truvari/wiki/Development#git-flow) for details on how branching is handled.

Docker
======
See [Development](https://github.com/spiralgenetics/truvari/wiki/Development#docker) for details on building a docker container.
