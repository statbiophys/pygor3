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
    parser.add_argument("-g", "--path_ref_genome", dest="path_ref_genome", help='Directory where genome references are stored: genomicDs.fasta,  genomicJs.fasta,  genomicVs.fasta,  J_gene_CDR3_anchors.csv,  V_gene_CDR3_anchors.csv', default='./ref_genome')
    # parser.add_argument("--set_genomicVs", dest="genomicVs", help='V gene templates fasta file path')
    # parser.add_argument("--set_genomicDs", dest="genomicDs", help='D gene templates fasta file path')
    # parser.add_argument("--set_genomicJs", dest="genomicJs", help='J gene templates fasta file path')

    parser.add_argument("-w", "--WD", dest="working_directory", help="Path where files gonna be created.", default='./')

    # To load all data files use batchname
    parser.add_argument("-b", "--batch", dest="batch",
                        help='Batchname to identify run. If not set random name is generated')
    # parser.add_option(dest="path_ref_genome")
    # parser.add_argument("-o", "--output", dest="output", help='filename of csv file to export data')

    args = parser.parse_args()

    # load a task object to get all data in one structure
    task = p3.IgorTask()
    task.igor_batchname = args.batch #"sample"  # options.batchname
    task.igor_wd = args.working_directory #"./test"

    # TODO: IF MODEL DIRECTORY PATH IS SET THEN COMPLETE igor_model_{parms,marginals}_file path
    # IF MODEL_PATH PROVIDED THEN
    task.igor_model_dir_path = args.model_path
    task.igor_model_parms_file = task.igor_model_dir_path+"/models/model_parms.txt"
    task.igor_model_marginals_file = task.igor_model_dir_path+"/models/model_marginals.txt"
    task.igor_path_ref_genome = task.igor_model_dir_path+"/ref_genome/"

    # IF GENOME PATH PROVIDED
    task.igor_path_ref_genome = args.path_ref_genome # "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/demo/thi/genomics_repseqio_F"
    
    task.update_batch_filenames()
    task.load_IgorRefGenome()
    print(task.igor_fln_indexed_sequences)
    print(task.genomes.df_genomicVs) # is where all this data is collected
    print(task.genomes.df_genomicDs)
    print(task.genomes.df_genomicJs)

    # Now create database to save genome references and also the data
    print(task.igor_fln_db)
    task.create_db()
    task.load_db_from_indexed_sequences()
    task.load_db_from_indexed_cdr3()
    task.load_db_from_genomes()
    task.load_db_from_alignments()
    task.load_db_from_models()

    task.load_db_from_bestscenarios()

    # seq_index = 0
    # indexed_sequence = task.igor_db.get_IgorIndexedSeq_By_seq_index(seq_index)
    # indexed_sequence.offset = 0
    # print(indexed_sequence)

    # genomes = p3.IgorRefGenome()
    # genomes.path_ref_genome =
    # genomes.update_fln_names()
    # df_V, df_J = genomes.load_dataframes()
    #
    # task.igor_wd = "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/demo/thi/test"
    # task.igor_batchname = "sample"


    # db = p3.IgorSqliteDB()
    # flnIgorDB = task.igor_batchname + ".db"
    # db.flnIgorSQL = "/home/alfaceor/Dropbox/PosDoc/IGoR/dev/MyGithub/pygor3/pygor3/IgorDB.sql"
    # db.createSqliteDB(flnIgorDB)
    # # db.load_VDJ_Database(
    # #     task.igor_fln_indexed_sequences,
    # #     genomes.fln_genomicVs,
    # #     genomes.fln_genomicDs,
    # #     genomes.fln_genomicJs,
    # #     task.igor_fln_align_V_alignments,
    # #     task.igor_fln_align_D_alignments,
    # #     task.igor_fln_align_J_alignments)
    #
    # seq_index = 0
    # indexed_sequence = db.get_IgorIndexedSeq_By_seq_index(seq_index)
    # indexed_sequence.offset = 0
    #
    # best_v_align_data = db.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
    # best_d_align_data = db.get_best_IgorAlignment_data_By_seq_index('D', indexed_sequence.seq_index)
    # best_j_align_data = db.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)
    #
    # vdj_naive_alignment = {'V': best_v_align_data,
    #                        'D': best_d_align_data,
    #                        'J': best_j_align_data}
    #
    # ofile = open(task.igor_batchname + '__' + str(indexed_sequence.seq_index) + '_na.fasta', 'w')
    # str_fasta = generate_str_fasta(indexed_sequence, vdj_naive_alignment)
    # ofile.write(str_fasta)
    # ofile.close()
    #
    # ofile = open('sequences2align.fasta', 'w')
    # str_fasta = generate_str_fasta_simple(indexed_sequence, vdj_naive_alignment)
    # ofile.write(str_fasta)
    # ofile.close()

    # import xarray as xr
    # import numpy as np
    # list_nt_lbl = ['A', 'C', 'G', 'T', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N']
    # da_heavy_pen_nuc44_vect = xr.DataArray(np.array(p3.utils.heavy_pen_nuc44_vect).reshape(15, 15), \
    #                                dims=('x', 'y'))
    # print(len(list_nt_lbl))
    # strDim = 'x'
    # da_heavy_pen_nuc44_vect[strDim] = range(len(list_nt_lbl))
    # strCoord = 'lbl__' + strDim
    # da_heavy_pen_nuc44_vect[strCoord] = (strDim, list_nt_lbl)
    #
    # strDim = 'y'
    # da_heavy_pen_nuc44_vect[strDim] = range(len(list_nt_lbl))
    # strCoord = 'lbl__' + strDim
    # da_heavy_pen_nuc44_vect[strCoord] = (strDim, list_nt_lbl)
    # print(da_heavy_pen_nuc44_vect)
    # import pygor3.utils as utils
    # print(utils.da_heavy_pen_nuc44_vect)

    #
    #  =

    # labels = self.parms.Event_dict[key]['value'].values
    #

    # strDim = 'y'
    # self.xdata[key][strDim] = range(len(self.xdata[key][strDim]))
    # strCoord = 'lbl__' + strDim
    # self.xdata[key][strCoord] = (strDim, labels)
    # xr.DataArray(aaa, dims=)
    #

    # 1. Print only sequence with corrected offset
    # 2. identify mistmatches
    # 3. show only aligment sequences

    """
    # print(seq_str)
    input_seq_align = p3.IgorAlignment_data()
    input_seq_align.seq_index = input.seq_index
    input_seq_align.strGene_seq = input.sequence
    input_seq_align.offset = 0
    input_seq_align.strGene_name = str( input.seq_index )
    input_seq_align.strGene_class = "InputSequence"
    input_seq_align.length = len(input_seq_align.strGene_seq)

    print("=" * 70)
    print(input_seq_align.strGene_seq)
    print("-" * 70)
    # print(best_v_align_data.strGene_seq)
    # print(best_d_align_data.strGene_seq)

    str_gap = '-'
    str_prefix = str_gap * (best_j_align_data.offset)
    str_tmp = str_prefix+best_j_align_data.strGene_seq
    print(str_tmp)
    print(best_j_align_data)
    # list(str_tmp)
    # from offset to offset_5_p
    str_tmp_2 = str_tmp[best_j_align_data.offset_5_p+1:best_j_align_data.offset_3_p+1]
    print(len(str_tmp_2))
    print(str_tmp_2)


    print("-" * 70)
    fasta_data = [input_seq_align, best_v_align_data, best_d_align_data, best_j_align_data]
    min_offset = min(fasta_data, key=lambda x: x.offset).offset
    print(min_offset)

    seq_list = list(input_seq_align.strGene_seq)
    vdj_mismatches = best_v_align_data.mismatches+best_d_align_data.mismatches+best_j_align_data.mismatches
    for ii in vdj_mismatches:
        seq_list[ii] = seq_list[ii].lower()
    input_seq_align.strGene_seq = ''.join(seq_list)
    fasta_data[0] = input_seq_align

    for item in fasta_data:
        str_prefix = 'X'*(item.offset - min_offset)
        seq_tmp = str_prefix + item.strGene_seq
        print(seq_tmp)
        aaa = seq_tmp[(item.offset_5_p - min_offset):(item.offset_3_p - min_offset)+1]
        print(item.strGene_class, item.offset_5_p, item.offset_3_p)
        print(aaa)

    # load Alignments from file


    # load Anchors with anchors extract CDR3 and naive alignment

    """


if __name__ == "__main__":
    main()
