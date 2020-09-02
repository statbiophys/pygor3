#!/usr/bin/env python3
# TODO: Script to load database given a collection of IGoR files
# 1. genomes
# 2. models
# 3. batch
# 3.1. indexed_sequences
# 3.2. alignments
# 3.3. best scenarios
# 3.4. pgen
# 3.5. coverage


import pygor3 as p3
# from optparse import OptionParser
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser()

    igor_models = parser.add_argument_group('IGoR default models')

    # Use IGoR default model
    igor_models.add_argument("-s", "--species", dest="species", help='Igor species')
    igor_models.add_argument("-c", "--chain", dest="chain", help='Igor chain')  # , type=str, choices=['TRB', 'TRA'])

    # Or specify directory in IGoR structure.
    parser.add_argument("-M", "--model_path", dest="model_path",
                        help='IGoR model directory path, this path include ref_genomes and model_parms')

    # parser.add_argument("-m", "--set_custom_model", dest="model", nargs=2, help='IGoR model_params.txt')
    parser.add_argument("-g", "--path_ref_genome", dest="path_ref_genome", help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv') #, default='./ref_genome')
    # parser.add_argument("--set_genomicVs", dest="genomicVs", help='V gene templates fasta file path')
    # parser.add_argument("--set_genomicDs", dest="genomicDs", help='D gene templates fasta file path')
    # parser.add_argument("--set_genomicJs", dest="genomicJs", help='J gene templates fasta file path')

    parser.add_argument("-w", "--WD", dest="working_directory", help="Path where files gonna be created.", default='./')

    # To load all data files use batchname
    parser.add_argument("-b", "--batch", dest="batch",
                        help='Batchname to identify run. If not set random name is generated', required=True)
    # parser.add_option(dest="path_ref_genome")
    # parser.add_argument("-o", "--output", dest="output", help='filename of csv file to export data')

    args = parser.parse_args()
    print(args)

    # load a task object to get all data in one structure
    task = p3.IgorTask()

    # Now create database to save genome references and also the data
    task.igor_batchname = args.batch  # "sample"  # options.batchname
    task.igor_wd = args.working_directory  # "./test"
    task.update_batch_filenames()
    print(task.igor_fln_db)
    task.create_db()

    task.load_db_from_indexed_sequences()
    task.load_db_from_indexed_cdr3()

    if (args.species is not None) and (args.chain is not None):
        print("species : ", args.species, " chain: ", args.chain )
        try:
            task.run_datadir()
            task.igor_species = args.species
            task.igor_chain = args.chain
            igor_model_path = task.igor_models_root_path + task.igor_species + "/" + p3.igor_option_path_dict[task.igor_chain] + "/"
            # igor_model_path = task.igor_models_root_path + task.igor_species + "/" \
            #                 + p3.igor_option_path_dict[task.igor_chain] + "/"
            args.model_path = igor_model_path
            print(args.model_path)
        except Exception as e:
            print("Default model not found!")
            print(e)

    # IF MODEL_PATH PROVIDED THEN
    # task.igor_model_dir_path = args.model_path
    if (args.model_path is not None):
        try:
            task.update_model_filenames(args.model_path)
            # task.igor_model_parms_file = task.igor_model_dir_path+"/models/model_parms.txt"
            # task.igor_model_marginals_file = task.igor_model_dir_path+"/models/model_marginals.txt"
            # task.igor_path_ref_genome = task.igor_model_dir_path+"/ref_genome/"
            # task.load_IgorModel()
            args.path_ref_genome = task.igor_path_ref_genome
            task.load_db_from_models()
        except Exception as e:
            print("Couldn't load models to database")
            print(e)

    # IF GENOME PATH PROVIDED
    if args.path_ref_genome is not None:
        try:
            task.igor_path_ref_genome = args.path_ref_genome  # "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/demo/thi/genomics_repseqio_F"
            task.load_IgorRefGenome()
            # print(task.igor_fln_indexed_sequences)
            # print(task.genomes.df_genomicVs) # is where all this data is collected
            # print(task.genomes.df_genomicDs)
            # print(task.genomes.df_genomicJs)

            task.load_db_from_genomes()
        except Exception as e:
            print("Couldn't load genome templates to database")
            print("ERROR: ", e)

    task.load_db_from_alignments()

    task.load_db_from_bestscenarios()

    task.load_db_from_pgen()

if __name__ == "__main__":
    main()
