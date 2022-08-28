# Truvari API
Many of the helper methods/objects are documented such that developers can reuse truvari in their own code. To see developer documentation, visit [readthedocs](https://truvari.readthedocs.io/en/latest/).

Documentation can also be seen using
```python
import truvari
help(truvari)
```

# docker

A Dockerfile exists to build an image of Truvari. To make a Docker image, clone the repository and run
```bash
docker build -t truvari .
```

You can then run Truvari through docker using
```bash
docker run -v `pwd`:/data -it truvari
```
Where `pwd` can be whatever directory you'd like to mount in the docker to the path `/data/`, which is the working directory for the Truvari run. You can provide parameters directly to the entry point.
```bash
docker run -v `pwd`:/data -it truvari anno svinfo -i example.vcf.gz
```

If you'd like to interact within the docker container for things like running the CI/CD scripts
```bash
docker run -v `pwd`:/data --entrypoint /bin/bash -it truvari
```
You'll now be inside the container and can run FuncTests or run Truvari directly
```bash
bash repo_utils/truvari_ssshtests.sh
truvari anno svinfo -i example.vcf.gz
```

# CI/CD

Scripts that help ensure the tool's quality. Extra dependencies need to be installed in order to run Truvari's CI/CD scripts. 

```bash
pip install pylint anybadge coverage
```

Check code formatting with 
```bash
python repo_utils/pylint_maker.py
```
We use [autopep8](https://pypi.org/project/autopep8/) (via [vim-autopep8](https://github.com/tell-k/vim-autopep8)) for formatting.

Test the code and generate a coverage report with 
```bash
bash repo_utils/truvari_ssshtests.sh
```

Truvari leverages github actions to perform these checks when new code is pushed to the repository. We've noticed that the actions sometimes hangs through no fault of the code. If this happens, cancel and resubmit the job. Once FuncTests are successful, it uploads an artifact of the `coverage html` report which you can download to see a line-by-line accounting of test coverage.

# git flow

To organize the commits for the repository, we use [git-flow](https://danielkummer.github.io/git-flow-cheatsheet/). Therefore, `develop` is the default branch, the latest tagged release is on `master`, and new, in-development features are within `feature/<name>`

When contributing to the code, be sure you're working off of develop and have run `git flow init`.
 
# versioning

Truvari uses [Semantic Versioning](https://semver.org/). As of v3.0.0, a single version is kept in the code under `truvari/__init__.__version__`. We try to keep the suffix `-dev` on the version in the develop branch. When cutting a new release, we may replace the suffix with `-rc` if we've built a release candidate that may need more testing/development. Once we've committed to a full release that will be pushed to PyPi, no suffix is placed on the version.

# docs

The github wiki serves the documentation most relevant to the `develop/` branch. When cutting a new release, we freeze and version the wiki's documentation with the helper utility `docs/freeze_wiki.sh`.

# Creating a release
Follow these steps to create a release

0) Bump release version
1) Run tests locally
2) Update API Docs
3) Freeze the Wiki
4) Ensure all code is checked in
5) Do a [git-flow release](https://danielkummer.github.io/git-flow-cheatsheet/)
6) Use github action to make a testpypi release
7) Check test release
```bash
python3 -m venv test_truvari
python3 -m pip install --index-url https://test.pypi.org/simple --extra-index-url https://pypi.org/simple/ truvari
```
8) Use GitHub action to make a pypi release
9) Change Updates Wiki
10) Download release-tarball.zip from step #8’s action
11) Create release (include #9) from the tag
12) Checkout develop and Bump to dev version and README ‘commits since’ badge