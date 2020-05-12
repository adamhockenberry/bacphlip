# BACPHLIP - a bacteriophage lifestyle prediction tool

*Adam J. Hockenberry and Claus O. Wilke*

[![PyPI version](https://badge.fury.io/py/bacphlip.svg)](https://badge.fury.io/py/bacphlip)
![PyPI - Downloads](https://img.shields.io/pypi/dm/bacphlip)
![PyPI - License](https://img.shields.io/pypi/l/bacphlip)
[![Build Status](https://travis-ci.com/adamhockenberry/bacphlip.svg?branch=master)](https://travis-ci.com/adamhockenberry/bacphlip)
[![Coverage Status](https://img.shields.io/codecov/c/github/adamhockenberry/bacphlip/master.svg)](https://codecov.io/github/adamhockenberry/bacphlip?branch=master)

**Reference:**

(eventual manuscript reference here)

## Installation

You can install BACPHLIP with pip:
```
pip install bacphlip
```

Alternatively, users can pull/download the latest github repository, navigate to the directory where BACPHLIP was downloaded and run:
```
pip install .
```

BACPHLIP has several required dependencies outside of the standard library: [biopython](https://pypi.org/project/biopython/), [pandas](https://pypi.org/project/pandas/), [joblib](https://pypi.org/project/joblib/), and [scikit-learn](https://pypi.org/project/scikit-learn/).

Additionally, users are required to install the [HMMER3 software suite](http://hmmer.org/) (in addition to the installation routes listed on the HMMER3 website we note that this tool can also be installed via [conda](https://anaconda.org/bioconda/hmmer)). By default, BACPHLIP assumes that HMMER3 is installed in the system path, but local paths may be provided as run-time flags (see below). 

## Examples

The most straightforwad usage of BACPHLIP is as a command line tool. Assuming that `/valid/path/to/a/genome.fasta` exists, you can call BACPHLIP with the command:
```
python -m bacphlip -i /valid/path/to/a/genome.fasta
```

This command should create 4 seperate files in the path of the target `genome.fasta` with `genome.fasta.bacphlip` containing the final model predictions (tab-separated format) in terms of probability of the input phage being either "Virulent" or "Temperate" (the other files append `.6frame`, `.hmmsearch`, and `.hmmsearch.tsv` to the genome file). Attempting to run this command a second time, assuming the first worked, should create an error since the output files already exist. This behavior can be altered with a flag to force overwrite the files:
```
python -m bacphlip -i /valid/path/to/a/genome.fasta -f 
```

Finally, a path to a local HMMER3 install (specifically, the `hmmsearch` tool) can be specified in the command line:
```
python -m bacphlip -i /valid/path/to/a/genome.fasta --local_hmmsearch /valid/path/to/hmmsearch
```


## Next steps

We have several planned next steps, including:
1. adding a tutorial for library usage as a jupyter notebook in a forthcoming `examples` folder. 
2. adding the ability to run the pipeline in a "quiet" mode
3. adding a flag for batch input of sequences. 

## Misc

The software is provided to you under the MIT license (see file `LICENSE.txt`).
The most up-to-date version of this software is available at
https://github.com/adamhockenberry/bacphlip.

The development of `BACPHLIP` is provided in a separate repository for transparency. See [bacphlip-model-dev](https://github.com/adamhockenberry/bacphlip-model-dev).

## Contributing

Pull requests addressing errors or adding new functionalities are welcome on GitHub. However, to be accepted, contributions must pass the `pytest` unit tests. 
