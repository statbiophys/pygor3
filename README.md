# pygor3
pygor3 is a python3 library to manipulate IGoR inputs/outputs and easyly make plots and scripts tasks.
pygor3 repository
synopsis

# Installation
0. (Optional) Install conda https://docs.conda.io/en/latest/ (or anaconda https://www.anaconda.com/) and create (or use ) a virtual environment
  $ conda create --name pygor3 python=3.7
  $ conda activate pygor3

1. Install requirements:

$ pip3 install -r requirements.txt

2. Install pygor3 package from this directory

$ pip3 install -e .

# simple usage.
## pygor3-load_database 

Script to load genomes references and alignments in a single file

usage: pygor3-load_database [-h] [-s SPECIES] [-c CHAIN] [-M MODEL_PATH]
                            [-g PATH_REF_GENOME] [-w WORKING_DIRECTORY]
                            [-b BATCH]

optional arguments:
  -h, --help            show this help message and exit
  -M MODEL_PATH, --model_path MODEL_PATH
                        IGoR model directory path, this path include
                        ref_genomes and model_parms
  -g PATH_REF_GENOME, --path_ref_genome PATH_REF_GENOME
                        Directory where genome references are stored:
                        genomicDs.fasta, genomicJs.fasta, genomicVs.fasta,
                        J_gene_CDR3_anchors.csv, V_gene_CDR3_anchors.csv
  -w WORKING_DIRECTORY, --WD WORKING_DIRECTORY
                        Path where files gonna be created.
  -b BATCH, --batch BATCH
                        Batchname to identify run. If not set random name is
                        generated

IGoR default models:
  -s SPECIES, --species SPECIES
                        Igor species
  -c CHAIN, --chain CHAIN
                        Igor chain


# Scripts
pygor3-model_export : Human readable tab-separated models in different files (http://physics.princeton.edu/~ccallan/TCRPaper/results/event_distributions.xls,) 
Usage: pygor3-model_export [options]

Options:
  -h, --help            show this help message and exit
  -s SPECIES, --species=SPECIES
                        Igor species
  -c CHAIN, --chain=CHAIN
                        Igor chain
  -b BATCH, --batch=BATCH
                        Batchname to identify run. If not set random name is
                        generated
  -o OUTPUT, --output=OUTPUT
                        filename of csv file to export data
  -p MODEL_PARAMS, --model_params=MODEL_PARAMS
                        IGoR model_params.txt
  -m MODEL_MARGINALS, --model_marginals=MODEL_MARGINALS
                        IGoR model_marginals.txt


pygor3-plot_marginals : Export csv file with real probability marginals from IGoR models.
Usage: pygor3-plot_marginals [options]

Options:
  -h, --help            show this help message and exit
  -s SPECIES, --species=SPECIES
                        Igor species
  -c CHAIN, --chain=CHAIN
                        Igor chain
  -o OUTPUT, --output=OUTPUT
                        filename of csv file to export data
  -p MODEL_PARAMS, --model_params=MODEL_PARAMS
                        IGoR model_params.txt
  -m MODEL_MARGINALS, --model_marginals=MODEL_MARGINALS
                        IGoR model_marginals.txt



[dev]
pygor3-pgen_sequences: Simple script to get the pgen of input sequences given an output filename, no batch required but optional.
pygor3-infer: A script to run any model given gene templates Is incomplete only works for VDJ for the moment. I found a bug for VJ I'm solving it :|
pygor3-bs_export: Given an IGoR batchname, species and chain it captures the model and convert the values of best scenarios to explicit values and not indexes.

