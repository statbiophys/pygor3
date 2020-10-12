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

1. Download new gene templates from IGMT website:
$ pygor3-cli imgt get-ref-genome -t VDJ --imgt-species Homo+sapiens --imgt-chain TRB

2. Make a new model from a ref_genome directory

$ pygor3-cli -M models/Homo+sapiens/TRB/ model create -t VDJ

3. Explore default model

$ pygor3-cli -M models/Homo+sapiens/TRB/ model plot

4. Infer model
$ pygor3-cli -M HumanTRB/ igor-infer -i test_seqs.csv -o myfile.db

5. Evaluate model
$ pygor3-cli -D myfile.db igor-evaluate -i test_seqs.csv -o mewdatabase.db

# if -o option is None then use the same myfile.db, so the sequences to use should be 

Show the tables in database:
$ pygor3-cli db-ls -D myfile.db

Delete elements in database:
$ pygor3-cli db-rm -D myfile.db --sequences --model --alignments

Extract elements in database:
$ pygor3-cli db-extract -D myfile.db --sequences --model --alignments -o newfile.db

Attach elements in database:
$ pygor3-cli db-attacht -D myfile.db --sequences --model --alignments -i another.db






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


## pygor3-load_database 

usage: pygor3-naive_align [-h] [-o OUTPUT] [-D DATABASE]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        filename of csv file to export data
  -D DATABASE, --database DATABASE
                        Igor database created with database script.


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







# pygor command line

## Common options
Usage: pygor [OPTIONS] COMMAND [ARGS]...

Options:
  -s, --set_igor_species TEXT     Species in IGoR's format: human, mouse
  -c, --set_igor_chain TEXT       Chain in IGoR's format, e.g. alpha, beta,
                                  TRB, TRA

  -m, --set_igor_model <model_parms.txt> <model_marginals.txt>
                                  IGoR model_params.txt and
                                  model_marginals.txt filenames.

  -M, --set_model_path <model_directory_path>
                                  IGoR model directory path, this path include
                                  ref_genomes and model_parms

  -g, --set_path_ref_genome TEXT  Directory where genome references are
                                  stored: genomicDs.fasta,  genomicJs.fasta,
                                  genomicVs.fasta,  J_gene_CDR3_anchors.csv,
                                  V_gene_CDR3_anchors.csv

  -w, --set_wd TEXT               Sets the working directory to path
                                  [default: ./]

  -b, --set_batch TEXT            Sets batchname to identify run. If not set
                                  random name is generated

  -D, --set_database TEXT         Igor database created with database script.
  --help                          Show this message and exit.

Commands:
  db-attach         Attach tables to database.
  db-export         Export database model in igor formatted files
  db-ls             List tables in database by groups and show number of...
  db-rm             Delete tables in database by groups.
  db-test           Manipulations of models
  igor-align        IGoR's call to aligns
  igor-evaluate     IGoR's call to evaluate input sequences
  igor-generate     IGoR's call to generate sequences
  igor-infer        IGoR's call to infer model from input sequences and...
  igor-pgen         IGoR's call to calculate pgen of input sequences
  igor-read-seqs    IGoR's call to read_seqs
  igor-scenarios    IGoR's call to get best scenarios.
  imgt-get-genomes  Get genomes from imgt website of specifing species and...
  model-create      Make a new default model VJ or VDJ with uniform...
  model-export      Export IGoR's models from txt (model_parms.txt,...
  model-plot        Plot real marginals of the bayesian network events


## pygor subcommands

### pygor imgt

Usage: pygor imgt-get-genomes [OPTIONS]

  Get genomes from imgt website of specifing species and chain in imgt
  format.

Options:
  --info                          List species and chain avialable in imgt
                                  website.

  -t, --recombination_type [VJ|VDJ]
                                  Igor recombination type.
  --imgt-species TEXT             IMGT species name for name specifications
                                  run imgt --info.

  --imgt-chain TEXT               IMGT chain name e.g. TRA, TRB.
  --help                          Show this message and exit.


### pygor model-create
Usage: pygor model-create [OPTIONS]

  Make a new default model VJ or VDJ with uniform probability distribution

Options:
  -t, --recombination_type [VJ|VDJ]
                                  Igor recombination type.
  --help                          Show this message and exit.


### pygor model-plot
Usage: pygor model-plot [OPTIONS]

  Plot real marginals of the bayesian network events

Options:
  -o, --output-prefix TEXT  Prefix to pdf files with model plots.
  --help                    Show this message and exit.


### pygor model-export
Usage: pygor model-export [OPTIONS]

  Export IGoR's models from txt (model_parms.txt, model_marginals.txt) files
  to db viceversa

Options:
  --from-txt TEXT...  Export Igor's model from txt files model_parms.txt and
                      model_marginals.txt.

  --from-db TEXT      Export Igor's model from database file.
  --to-txt TEXT...    Output filename of Igor recombination model to
                      <model_parms.txt> <model_marginals.txt>.

  --to-db TEXT        Output filename of Igor recombination model to
                      <model.db>.

  --help              Show this message and exit.



### pygor igor-evaluate
Usage: pygor igor-evaluate [OPTIONS]

  IGoR's call to evaluate input sequences

Options:
  -i, --input-sequences TEXT  Input sequences in FASTA, TXT or CSV formats.
  -o, --output-db TEXT        Output database file.
  --help                      Show this message and exit.




### pygor db-attach
Usage: pygor db-attach [OPTIONS]

  Attach tables to database.

Options:
  --from-db TEXT                  Database copy source filename.
  --from-batch TEXT               Database copy source filename.
  --from-genome-dir TEXT          Database copy source filename.
  --from-model-path TEXT          Database copy source filename.
  --igor-model-dir <model.db>     IGoR model database file.
  --igor-model-parms <model_parms.txt>
                                  IGoR model parms (or params) file.
  --igor-model-marginals <model_marginals.txt>
                                  IGoR model marginals file.
  --scenarios TEXT                If --from-db no need to add filename
  --pgen TEXT                     If --from-db no need to add filename
  --genomes TEXT                  Copy V, (D) and J genetic data to database
  --genomesV TEXT                 Copy just V genomes tables to database.
  --genomesD TEXT                 Copy just D genomes tables to database.
  --genomesJ TEXT                 Copy just J genomes tables to database.
  --genomesCDR3 TEXT              Copy just CDR3 anchors tables to database.
  --alignments TEXT               Copy all available alignments tables to
                                  database.

  --alignmentsV TEXT              Copy V alignments tables to database.
  --alignmentsD TEXT              Copy D alignments tables to database.
  --alignmentsJ TEXT              Copy J alignments tables to database.
  --alignmentsCDR3 TEXT           Copy indexed cdr3 table to database.
  --help                          Show this message and exit.



