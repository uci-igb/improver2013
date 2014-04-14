# Overview

Peter Sadowski and Michael Zeller (Team IGB) submission and scripts for sbv IMPROVER 2013 subchallenge 2

# Usage

Clone the git repository and also run the following to get the
submodules `pylearn2` and `cybtpy`.

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

# Relevant scripts

The neural network training is performed using the Pylearn2 submodule, with the script and readme located in
improver2013/opt/pylearn2/pylearn2/scripts/improver2013/

The analysis scripts for the competition are in the `scripts` dir, and
they work with files from the `data` and `out` directories.


