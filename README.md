# Overview

Peter Sadowski and Michael Zeller (Team IGB) submission and scripts for sbv IMPROVER 2013 subchallenge 2

# Setup

Clone the git repository and also run the following to get the
git submodules `pylearn2` and `cybtpy`.

```
git submodule init
git submodule update
```

You will also need to extract the competition provided data from the data directory, for example:

```
cd data
unzip *.zip
```

After extracting, copy the files from the `SBV_STC_subchallenge1_gold_standard` folder
into the `SBV_STC_subchallenge2` folder.

# Requirements

Code has been tested on `CentOS 6` using `Python 2.7.6` and `R 3.0.1`.

Required R packages are `gplots` and `VGAM`.

# Relevant scripts

The neural network training is performed using the Pylearn2 submodule, with the script and readme located in
`improver2013/opt/pylearn2/pylearn2/scripts/improver2013/`

The pre- and post-analysis scripts for the competition are in the `scripts` dir, and
they work with files from the `data` and `out` directories.

## Pre-training analysis

`process-microarray.r`: Processes the GEx and P data and includes the
batch number to construct the training data as a single file.

`transform-data.r`: Implements the transformations of the GEx data to
be clipped from 0 to 4.

## Post-training analysis

`convert-to-significance-subchallenge2.r`: Creates the submission file
and figures using combinations of t-tests (bayesian, standardard),
transformations (fisher z, logistic), and combinations of training
data (using GEx, rat P only, rat and human P without GEx).
