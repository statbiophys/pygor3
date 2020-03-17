# pygor3
pygor3 is a python3 library to manipulate IGoR inputs/outputs and easyly make plots and scripts tasks.
pygor3 repository
synopsis

# Installation
0. (Optional) Install conda https://docs.conda.io/en/latest/ (or anaconda https://www.anaconda.com/) and/or create a virtual environment
1. Install requirements:

$ pip3 install -r requirements.txt

2. Install pygor3 package from this directory

$ pip3 install -e .

# simple usage.

igor-model_export : Human readable tab-separated models in different files (http://physics.princeton.edu/~ccallan/TCRPaper/results/event_distributions.xls,) 
igor-pgen_sequences: Simple script to get the pgen of input sequences given an output filename, no batch required but optional.
igor-infer: A script to run any model given gene templates Is incomplete only works for VDJ for the moment. I found a bug for VJ I'm solving it :|
igor-bs_export: Given an IGoR batchname, species and chain it captures the model and convert the values of best scenarios to explicit values and not indexes.

