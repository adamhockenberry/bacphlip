# BACPHLIP - a bacteriophage lifestyle prediction tool

*Adam J. Hockenberry and Claus O. Wilke*

[![PyPI version](https://badge.fury.io/py/bacphlip.svg)](https://badge.fury.io/py/bacphlip)
![PyPI - Downloads](https://img.shields.io/pypi/dm/bacphlip)
![PyPI - License](https://img.shields.io/pypi/l/bacphlip)
[![Build Status](https://travis-ci.com/adamhockenberry/bacphlip.svg?branch=master)](https://travis-ci.com/adamhockenberry/bacphlip)
[![Coverage Status](https://img.shields.io/codecov/c/github/adamhockenberry/bacphlip/master.svg)](https://codecov.io/github/adamhockenberry/bacphlip?branch=master)

**Reference:**

Pre-print available at: <https://www.biorxiv.org/content/10.1101/2020.05.13.094805v1>


## Overview and important caveats

The BACPHLIP software is designed to test whether a given phage genome (`.fasta` formatted) is likely to be either temperate (lysogenic) or virulent (lytic). The software makes this determination by searching for a particular set of what are hypothesized to be "temperate-specific" protein domains. As such, the default assumption is that any given input file is a virulent (lytic) phage. Depending on the number and identity of various proteins that are found, this default assumption may be updated to indicate that the sequence is in fact temperate. **BACPHLIP does not perform any checks on whether the input sequence is even a phage.** Thus, random stretches of DNA will be called virulent phages (assuming that no relevant domains are found within the random sequence) not because there are any idications of the sequence being a virulent phage, but rather because no data overturns the starting assumption (that you provided the program with a phage). Similarly strange results will occur if you provide BACPHLIP with whole bacterial chromosomes, these will likely be called temperate phages simply because it's likely that several of the relevant "temperate domains" will be found somewhere within the chromosome. 

Finally, we recommend that users read through all documentation here as well as the manuscript (referenced above). Our software was trained on a dataset consisting almost entirely of phages from the order *Caudovirales*, most of which infect hosts in the orders *Actinobacteria*, *Gammaproteobacteria*, and *Bacilli*. We urge caution when using the software on species outside of these orders, but this fact may change as we update and expand training set data in future releases. 

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
bacphlip -i /valid/path/to/a/genome.fasta
```

This command should create 4 seperate files in the path of the target `genome.fasta` with `genome.fasta.bacphlip` containing the final model predictions (tab-separated format) in terms of probability of the input phage being either "Virulent" or "Temperate" (the other files append `.6frame`, `.hmmsearch`, and `.hmmsearch.tsv` to the genome file). Attempting to run this command a second time, assuming the first worked, should create an error since the output files already exist. This behavior can be altered with a flag to force overwrite the files:
```
bacphlip -i /valid/path/to/a/genome.fasta -f 
```

Finally, a path to a local HMMER3 install (specifically, the `hmmsearch` tool) can be specified in the command line:
```
bacphlip -i /valid/path/to/a/genome.fasta --local_hmmsearch /valid/path/to/hmmsearch
```

However, BACPHLIP can also be accessed and used as a python library. From a python interpreter simply type:
```
import bacphlip
bacphlip.run_pipeline('/valid/path/to/a/genome.fasta')
```

At present this is probably the easiest way to run BACPHLIP on a batch of input files:
```
import bacphlip
import glob
for infile_loc in glob.glob('/valid/path/to/a/set/of/files/*.fasta'):
    bacphlip.run_pipeline(infile_loc)
```

Finally, using BACPHLIP as a library makes individual functions available to the user in order to run and possibly troubleshoot single steps. I.e.:
```
import bacphlip
bacphlip.six_frame_translate( ... )
bacphlip.hmmsearch_py( ... )
bacphlip.process_hmmsearch( ... )
bacphlip.predict_lifestyle( ... )
```
Each function has a relevant set of arguments that should be clear from the docs. Running in this manner will give more flexibility with regard to file names and may prove useful to some users.

## Next steps

We have several planned next steps, including:
1. adding a tutorial for library usage as a jupyter notebook in a forthcoming `examples` folder. 
2. adding the ability to run the pipeline in a "quiet" mode
3. adding a flag for batch input of sequences. 
4. (insert your suggestion here)

## Misc

The software is provided to you under the MIT license (see file `LICENSE.txt`).
The most up-to-date version of this software is available at
https://github.com/adamhockenberry/bacphlip.

The development of `BACPHLIP` is provided in a separate repository for transparency. See [bacphlip-model-dev](https://github.com/adamhockenberry/bacphlip-model-dev).

## Contributing

Pull requests addressing errors or adding new functionalities are welcome on GitHub. However, to be accepted, contributions must pass the `pytest` unit tests. 
