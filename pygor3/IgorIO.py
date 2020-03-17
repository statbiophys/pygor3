#      Author: Quentin Marcou & Carlos Olivares
#
#  This source code is distributed as part of the IGoR software.
#  IGoR (Inference and Generation of Repertoires) is a versatile software to
#  analyze and model immune receptors generation, selection, mutation and all
#  other processes.
#   Copyright (C) 2019  Quentin Marcou & Carlos Olivares
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import pandas as pd
import xarray as xr
import networkx as nx

from .IgorDictionaries import *
from .IgorDefaults import *
from .IgorSqliteDB import *
from .IgorSqliteDBBestScenarios import *
### load IGoR sequences database
from pygor3 import rcParams
import subprocess

from .IgorSQL import *


### GENERIC FUNCTIONS
# Generation of label for a simple identification of genomic template sequence.
def genLabel(strName):
    aaa = strName.split("|")
    if len(aaa) > 1 :
        return aaa[1]
    else:
        return strName

def command_from_dict_options(dicto:dict):
    """ Return igor options from dictionary"""
    cmd = ''
    print()
    for key in dicto.keys():
        if dicto[key]['active']:
            cmd = cmd + " " + key + " " + dicto[key]['value']
            if dicto[key]['dict_options'] is not None:
                #print(key, dicto[key]['dict_options'])
                cmd = cmd + " " + command_from_dict_options(dicto[key]['dict_options'])
    return cmd

def run_command(cmd):
    """from http://blog.kagesenshi.org/2008/02/teeing-python-subprocesspopen-output.html
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = []
    while True:
        line = p.stdout.readline()
        line = line.decode("utf-8")
        stdout.append(line)
        print (line, end='')
        if line == '' and p.poll() != None:
            break
    return ''.join(stdout)

def run_command_no_output(cmd):
    """from http://blog.kagesenshi.org/2008/02/teeing-python-subprocesspopen-output.html
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return p

# FIXME: IT IS BETTER TO USE DECORATORS FOR VARIABLES LIKE igor_batchname and update the dependencies on that automatically?
class IgorTask:
    """
    This class should encapsulate all
    the input parameters and output files when IGoR run.
    """
    def __init__(self):
        self.igor_exec_path = ""
        self.igor_datadir = ""
        self.igor_modelspath = ""

        self.igor_modeldirpath = ""

        self.igor_wd = ""
        self.igor_batchname = ""
        self.igor_specie = ""
        self.igor_chain = ""
        self.igor_model_parms_file = ""
        self.igor_model_marginals_file = ""
        self.igor_read_seqs = ""
        self.igor_threads = ""

        # read
        self.igor_fln_indexed_sequences = ""
        self.igor_fln_indexed_CDR3 = ""
        # inference
        self.igor_fln_infer_final_marginals = ""
        self.igor_fln_infer_final_parms = ""
        # evaluate
        self.igor_fln_evaluate_final_marginals = ""
        self.igor_fln_evaluate_infer_final_parms = ""
        # output
        self.igor_fln_output_pgen = ""
        self.igor_fln_output_scenarios = ""
        self.igor_fln_output_coverage = ""

        self.igor_fln_db = ""

        # TODO: experimental dictionary to check status of igor batch associated files
        # almost each of these files correspond to a sql table
        self.batch_data = igor_batch_dict


        self.igor_db = None #IgorSqliteDB()
        self.igor_db_bs = None

        self.b_read_seqs = False
        self.b_align = False
        self.b_infer = False
        self.b_evaluate = False
        self.b_generate = False

        self.mdl = IgorModel()

        # FIXME: THIS OPTIONS SHOULD BE AT RC PARMS
        tmp_dict_options = {
                         '---thresh': {'active': False, 'value': '15', 'dict_options': {}},
                         '---matrix': {'active': False, 'value': 'path/to/file', 'dict_options': {}},
                         '---gap_penalty': {'active': False, 'value': 'X', 'dict_options': {}},
                         '---best_align_only': {'active': False, 'value': '', 'dict_options': {}}
                     }
        self.igor_align_dict_options = igor_align_dict_options

        self.igor_output_dict_options = igor_output_dict_options

        try:
            p = subprocess.Popen("head /dev/urandom | tr -dc A-Za-z0-9 | head -c10", shell=True, stdout=subprocess.PIPE)
            line = p.stdout.readline()
            self.igor_batchname = "dataIGoR"+line.decode("utf-8").replace('\n', '')
        except Exception as e:
            print(e)
            raise e

        try:
            p = subprocess.Popen("pwd", shell=True, stdout=subprocess.PIPE)
            line = p.stdout.readline()
            self.igor_wd = line.decode("utf-8").replace('\n', '')
        except Exception as e:
            print(e)
            raise e

        try:
            p = subprocess.Popen("which igor", shell=True, stdout=subprocess.PIPE)
            line = p.stdout.readline()
            self.igor_exec_path = line.decode("utf-8").replace('\n', '')
        except Exception as e:
            print(e)
            raise e

        try:
            self.run_datadir()
        except Exception as e:
            print(e)
            raise e

    def load_IgorModel(self):
        if (self.igor_specie == "" or self.igor_chain == ""):
            self.mdl = IgorModel(model_parms_file = self.igor_model_parms_file, model_marginals_file=self.igor_model_marginals_file)
        else :
            self.mdl = IgorModel.load_default(self.igor_specie, igor_option_path_dict[self.igor_chain])

    @classmethod
    def default_model(cls, specie, chain, model_parms_file=None, model_marginals_file=None):
        """Return an IgorTask object"""
        cls.igor_specie = specie
        cls.igor_chain = chain
        cls.igor_modeldirpath =  model_parms_file
        cls.run_datadir()

        if model_parms_file is None:
            cls.igor_model_parms_file = cls.igor_modelspath +"/" + cls.igor_specie + "/" + igor_option_path_dict[cls.igor_chain]

    def update_batch_filenames(self):
        self.igor_fln_indexed_sequences = self.igor_wd + "/aligns/" + self.igor_batchname + "_indexed_sequences.csv"
        self.igor_fln_indexed_CDR3 = self.igor_wd + "/aligns/" + self.igor_batchname + "_indexed_CDR3.csv"

        # inference
        tmpstr = self.igor_wd + "/" + self.igor_batchname + "_inference/"
        self.igor_fln_infer_final_parms = tmpstr + "final_parms.txt"
        self.igor_fln_infer_final_marginals = tmpstr + "final_marginals.txt"

        # evaluate
        tmpstr = self.igor_wd + "/" + self.igor_batchname + "_evaluate/"
        self.igor_fln_evaluate_final_parms = tmpstr + "final_parms.txt"
        self.igor_fln_evaluate_final_marginals = tmpstr + "final_marginals.txt"

        # output
        tmpstr = self.igor_wd + "/" + self.igor_batchname + "_output/"
        self.igor_fln_output_pgen = tmpstr + "Pgen_counts.csv"
        self.igor_fln_output_scenarios = tmpstr + "best_scenarios_counts.csv"
        self.igor_fln_output_coverage = tmpstr + "coverage.csv"

        # Set all files as not existing by default
        import os.path
        for file_id in  igor_file_id_list:
            self.batch_data[file_id]['status'] = os.path.isfile(self.batch_data[file_id]['filename'])
        # database
        self.igor_fln_db = self.igor_wd + "/" + self.igor_batchname+".db"

        tmp_prefix_aligns = self.igor_wd + "/aligns/" + self.igor_batchname
        self.batch_data['indexed_sequences']['filename'] = tmp_prefix_aligns  + "_indexed_sequences.csv"
        self.batch_data['indexed_CDR3']['filename'] = tmp_prefix_aligns + "_indexed_CDR3.csv"
        self.batch_data['aligns_V_alignments']['filename'] = tmp_prefix_aligns + "_V_alignments.csv"
        self.batch_data['aligns_D_alignments']['filename'] = tmp_prefix_aligns + "_D_alignments.csv"
        self.batch_data['aligns_J_alignments']['filename'] = tmp_prefix_aligns + "_J_alignments.csv"

        tmp_prefix = self.igor_wd + "/" + self.igor_batchname
        self.batch_data['infer_final_parms']['filename'] = tmp_prefix + "_inference/" + "final_parms.txt"
        self.batch_data['infer_final_marginals']['filename'] = tmp_prefix + "_inference/" + "final_marginals.txt"
        self.batch_data['evaluate_final_parms']['filename'] = tmp_prefix + "_evaluate/" + "final_parms.txt"
        self.batch_data['evaluate_final_marginals']['filename'] = tmp_prefix + "_evaluate/" + "final_marginals.txt"
        self.batch_data['output_pgen']['filename'] = tmp_prefix + "_output/" + "Pgen_counts.csv"
        self.batch_data['output_scenarios']['filename'] = tmp_prefix + "_output/" + "best_scenarios_counts.csv"
        self.batch_data['output_coverage']['filename'] = tmp_prefix + "_output/" + "coverage.csv"

    # def update_igor_filenames_by_modeldirpath(self, modeldirpath=None):
    #     if modeldirpath is None:
    #         modeldirpath = self.igor_modelspath + "/" + self.igor_specie + "/" + igor_option_path_dict[self.igor_chain]
    #
    #     self.batch_data['genomicVs']['filename'] = tmp_prefix
    #     self.batch_data['genomicDs']['filename'] = tmp_prefix
    #     self.batch_data['genomicJs']['filename'] = tmp_prefix
    #
    #     self.batch_data['V_gene_CDR3_anchors']['filename'] = tmp_prefix
    #     self.batch_data['J_gene_CDR3_anchors']['filename'] = tmp_prefix
    #
    #     self.batch_data['model_parms']['filename'] = tmp_prefix
    #     self.batch_data['model_marginals']['filename'] = tmp_prefix

    def update_batchname(self, batchname):
        self.igor_batchname = batchname
        self.update_batch_filenames()

    @classmethod
    def load_from_batchname(cls, batchname, wd=None):
        cls = IgorTask()
        if wd is None:
            cls.igor_wd = "."
        else:
            cls.igor_wd = wd
        cls.update_batchname(batchname)

        try:
            cls.run_datadir()
        except Exception as e:
            print(e)
            raise e

        return cls

    def run_demo(self):
        cmd = self.igor_exec_path+ " -run_demo"
        run_command(cmd)

    def run_datadir(self):
        cmd = self.igor_exec_path+ " -datadir"
        self.igor_datadir = run_command(cmd).replace('\n','')
        self.igor_modelspath = self.igor_datadir + "/models/"

    def run_read_seqs(self):
        "igor -set_wd $WDPATH -batch foo -read_seqs ../demo/murugan_naive1_noncoding_demo_seqs.txt"
        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        cmd = cmd + " -read_seqs " + self.igor_read_seqs
        # TODO: if self.igor_read_seqs extension fastq then convert to csv and copy and create the file in aligns. Overwrite if necesserasy
        print(cmd)
        run_command(cmd)
        self.b_read_seqs = True # FIXME: If run_command success then True

    def run_align(self):
        #"igor -set_wd ${tmp_dir} -batch ${randomBatch} -species
        # ${species} -chain ${chain} -align --all"
        if self.b_read_seqs is False:
            self.run_read_seqs()
        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        # TODO: USE COSTUM MODEL OR USE SPECIFIED SPECIES?
        # I think that the safests is to use the
        cmd = cmd + " -species " + self.igor_specie
        cmd = cmd + " -chain " + self.igor_chain
        cmd = cmd + " -align " + command_from_dict_options(self.igor_align_dict_options)
        #return cmd
        print(cmd)
        run_command(cmd)
        #run_command_no_output(cmd)
        self.b_align = True # FIXME: If run_command success then True

    def run_evaluate(self):
        #"igor -set_wd $WDPATH -batch foo -species human -chain beta
        # -evaluate -output --scenarios 10"
        if self.b_align is False:
            self.run_align()

        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        # TODO: USE COSTUM MODEL OR USE SPECIFIED SPECIES?
        # I think that the safests is to use the
        cmd = cmd + " -species " + self.igor_specie
        cmd = cmd + " -chain " + self.igor_chain
        # here the evaluation
        cmd = cmd + " -evaluate -output " + command_from_dict_options(self.igor_output_dict_options)
        #return cmd
        print(cmd)
        # FIXME: REALLY BIG FLAW USE DICTIONARY FOR THE SPECIE AND CHAIN
        self.mdl = IgorModel.load_default(self.igor_specie, igor_option_path_dict[self.igor_chain], modelpath=self.igor_modelspath)
        run_command(cmd)
        #run_command_no_output(cmd)
        #self.b_evaluate = True # FIXME: If run_command success then Truerun_infer

    def run_infer(self):
        #"igor -set_wd $WDPATH -batch foo -species human -chain beta
        # -evaluate -output --scenarios 10"
        if self.b_align is False:
            self.run_align()

        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        # TODO: USE COSTUM MODEL OR USE SPECIFIED SPECIES?
        # I think that the safests is to use the
        cmd = cmd + " -species " + self.igor_specie
        cmd = cmd + " -chain " + self.igor_chain
        # here the evaluation
        cmd = cmd + " -infer " + command_from_dict_options(self.igor_output_dict_options)
        #return cmd
        print(cmd)
        # FIXME: REALLY BIG FLAW USE DICTIONARY FOR THE SPECIE AND CHAIN
        self.mdl = IgorModel.load_default(self.igor_specie, igor_option_path_dict[self.igor_chain], modelpath=self.igor_modelspath)
        run_command(cmd)
        #run_command_no_output(cmd)
        self.b_infer = True # FIXME: If run_command success then True
        

    def run_clean_batch(self):
        cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_evaluate"
        run_command_no_output(cmd)
        cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_output"
        run_command_no_output(cmd)
        cmd = "rm " + self.igor_wd + "/aligns/" + self.igor_batchname + "*.csv"
        run_command_no_output(cmd)

    def load_VDJ_database(self, flnIgorSQL):
        self.flnIgorSQL = flnIgorSQL
        self.igor_db = IgorSqliteDB(flnIgorSQL)
        # FIXME :EVERYTHING
        flnIgorIndexedSeq = self.igor_wd+"/aligns/"+self.igor_batchname+"_indexed_sequences.csv"
        # FIXME PATH AND OPTIONS NEED TO BE CONSISTENT
        IgorModelPath = self.igor_modelspath + self.igor_specie + "/" \
                        + igor_option_path_dict[self.igor_chain] + "/"
        IgorRefGenomePath = IgorModelPath + "ref_genome/"

        flnVGeneTemplate = IgorRefGenomePath + "genomicVs.fasta"
        flnDGeneTemplate = IgorRefGenomePath + "genomicDs.fasta"
        flnJGeneTemplate = IgorRefGenomePath + "genomicJs.fasta"

        flnVGeneCDR3Anchors = IgorRefGenomePath + "V_gene_CDR3_anchors.csv"
        flnJGeneCDR3Anchors = IgorRefGenomePath + "J_gene_CDR3_anchors.csv"

        ### IGoR Alignments files
        flnVAlignments = self.igor_wd + "/aligns/" + self.igor_batchname + "_V_alignments.csv"
        flnDAlignments = self.igor_wd + "/aligns/" + self.igor_batchname + "_D_alignments.csv"
        flnJAlignments = self.igor_wd + "/aligns/" + self.igor_batchname + "_J_alignments.csv"

        ### IGoR ouptut files
        flnModelParms = IgorModelPath + "models/model_parms.txt"
        flnModelMargs = IgorModelPath + "models/model_marginals.txt"
        flnIgorBestScenarios = self.igor_wd + self.igor_batchname + "_output/best_scenarios_counts.csv"

        flnIgorDB = self.igor_batchname+".db"
        self.igor_db.createSqliteDB(flnIgorDB)
        self.igor_db.load_VDJ_Database(flnIgorIndexedSeq, \
                             flnVGeneTemplate, flnDGeneTemplate, flnJGeneTemplate, \
                             flnVAlignments, flnDAlignments, flnJAlignments)

        # ### load IGoR model parms and marginals.
        # # FIXME: THIS IS REDUNDANT IN SOME PLACE check it out.
        # self.mdl = IgorModel(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
        # mdlParms = IgorModel.Model_Parms(flnModelParms)  # mdl.parms
        # mdlMargs = IgorModel.Model_Marginals(flnModelMargs)  # mdl.marginals
        #
        # # load IGoR best scenarios file.
        # db_bs = IgorSqliteDBBestScenarios.IgorSqliteDBBestScenariosVDJ()
        # db_bs.createSqliteDB("chicagoMouse_bs.db")
        # db_bs.load_IgorBestScenariosVDJ_FromCSV(flnIgorBestScenarios)

    def load_VDJ_BS_database(self, flnIgorBSSQL):
        flnIgorBestScenarios = self.igor_wd+"/"+self.igor_batchname+"_output/best_scenarios_counts.csv"
        self.igor_db_bs = IgorSqliteDBBestScenariosVDJ(flnIgorBSSQL) #IgorDBBestScenariosVDJ.sql
        self.igor_db_bs.createSqliteDB(self.igor_batchname+"_bs.db")
        self.igor_db_bs.load_IgorBestScenariosVDJ_FromCSV(flnIgorBestScenarios)

    def get_pgen_pd(self):
        #load pgen file
        import pandas as pd
        df = pd.read_csv(self.igor_fln_output_pgen, sep=';')
        df = df.set_index('seq_index')
        df = df.sort_index()
        df_seq = pd.read_csv(self.igor_fln_indexed_sequences, sep=';')
        df_seq = df_seq.set_index('seq_index').sort_index()
        df_cdr3 = pd.read_csv(self.igor_fln_indexed_CDR3, sep=';')
        df_cdr3 = df_cdr3.set_index('seq_index').sort_index()
        df = df.merge(df_seq, left_index=True, right_index=True)
        df = df.merge(df_cdr3, left_index=True, right_index=True)
        return df





    ### IGOR INPUT SEQUENCES  ####


class IgorIndexedSequence:
    """
    Return a IgorIndexedSequence instance
    """
    def __init__(self, seq_index=-1, sequence=''):
        self.seq_index     = seq_index
        self.sequence      = sequence

    def __str__(self):
        return str(self.to_dict())

    def to_dict(self):
        """
        Return a IgorIndexedSequence instance as a python dictionary.
        """
        dictIndexedSequence       =  {
            "seq_index"     : self.seq_index  , \
            "sequence"      : self.sequence
            }
        
        return dictIndexedSequence

    @classmethod
    def load(cls, seq_index, sequence):
        cls = IgorIndexedSequence()
        try:
            cls.seq_index = seq_index
            cls.sequence  = sequence
        except Exception as e:
            print(e)
            raise e
        return cls


    @classmethod
    def load_FromCSVline(cls, csvline, delimiter=";"):
        """
        Return a IgorIndexedSequence instance from a line of IGoR indexed_sequences.csv file.
        :param csvline: String line of a csv file.
        :param delimiter: Character to delimitate csv file.
        :return: IgorIndexedSequence object
        """
        cls = IgorIndexedSequence()
        csvsplit = csvline.replace("\n", "").split(";")
        try:
            cls.seq_index    = int  (csvsplit[0])
            cls.sequence     = int  (csvsplit[1])
        except Exception as e:
            print(e)
            raise e
        return cls
    
    @classmethod
    def load_FromSQLRecord(cls, sqlRecord):
        """
        Return a IgorIndexedSequence instance from a database record accordingly.
        with the database specification.
        :param sqlRecord: sqlite record of one entry.
        :return: IgorIndexedSequence object.
        """
        cls = IgorIndexedSequence()
        try:
            cls.seq_index    = int  (sqlRecord[0])
            cls.sequence     = int  (sqlRecord[1])
        except Exception as e:
            print(e)
            raise e
        return cls

### IGOR ALIGNMENTS  ####
class IgorAlignment_data:
    def __init__(self):
        self.seq_index     = -1
        self.gene_id       = -1
        self.score         = -1
        self.offset        = 0
        self.insertions    = list()
        self.deletions     = list() 
        self.mismatches    = list()
        self.length        = 0
        self.offset_5_p    = 0
        self.offset_3_p    = 0
        
        self.strGene_name  = ""
        self.strGene_class = ""
        self.strGene_seq   = ""

    def to_dict(self):
        dictAlignment_data       =  {
            "seq_index"     : self.seq_index  , \
            "gene_id"       : self.gene_id    , \
            "score"         : self.score      , \
            "offset"        : self.offset     , \
            "insertions"    : self.insertions , \
            "deletions"     : self.deletions  , \
            "mismatches"    : self.mismatches , \
            "length"        : self.length     , \
            "offset_5_p"    : self.offset_5_p , \
            "offset_3_p"    : self.offset_3_p , \
            "strGene_name"  : self.strGene_name    , \
            "strGene_class" : self.strGene_class   , \
            "strGene_seq"   : self.strGene_seq
            }
        
        return dictAlignment_data
    
    @classmethod
    def load_FromCSVLine(cls, csvline, strGene_name="", delimiter=";"):
        #seq_index;gene_name;score;offset;insertions;deletions;mismatches;length;5_p_align_offset;3_p_align_offset
        cls = IgorAlignment_data()
        csvsplit = csvline.replace("\n", "").split(";")
        try:            
            cls.seq_index    = int  (csvsplit[0])
            cls.strGene_name = str  (csvsplit[1])
            cls.score        = float(csvsplit[2])
            cls.offset       = int  (csvsplit[3])
            cls.insertions   = eval (csvsplit[4])
            cls.deletions    = eval (csvsplit[5])
            cls.mismatches   = eval (csvsplit[6])
            cls.length       = int  (csvsplit[7])
            cls.offset_5_p   = int  (csvsplit[8])
            cls.offset_3_p   = int  (csvsplit[9])
        except Exception as e:
            print(e)
            raise e
        return cls

    @classmethod
    def load_FromSQLRecord(cls, sqlRecordAlign, strGene_name=""):
        """
        Return a IgorAlignment_data instance from a IgorSqlRecord.
        :param sqlRecordAlign: record of a sql database table.
        :param strGene_name: gene_name associated to the record.
        :return: IgorAlignment_data instance
        """
        cls = IgorAlignment_data()
        try:
            cls.seq_index    = int  (sqlRecordAlign[0])
            cls.gene_id      = int  (sqlRecordAlign[1])
            cls.score        = float(sqlRecordAlign[2])
            cls.offset       = int  (sqlRecordAlign[3])
            cls.insertions   = eval (sqlRecordAlign[4])
            cls.deletions    = eval (sqlRecordAlign[5])
            cls.mismatches   = eval (sqlRecordAlign[6])
            cls.length       = int  (sqlRecordAlign[7])
            cls.offset_5_p   = int  (sqlRecordAlign[8])
            cls.offset_3_p   = int  (sqlRecordAlign[9])
            # TODO: Bestway to retrieve the name of the gene_name
            if strGene_name == None:
                cls.strGene_name = str(cls.gene_id)
            else:
                cls.strGene_name = strGene_name
            return cls
        except Exception as e:
            print(e)
            raise e

### IGOR MODEL ####
class IgorModel:
    def __init__(self, model_parms_file=None, model_marginals_file=None):
        self.parms = IgorModel_Parms()
        self.marginals = IgorModel_Marginals()
        self.xdata = dict()
        self.metadata = dict()
        self.specie = ""
        self.chain = ""


        # check input files
        flag_parms = (model_parms_file is not None)
        flag_marginals = (model_marginals_file is not None)
        flag_xdata = (flag_parms and flag_marginals)

        if flag_parms:
            self.parms.read_model_parms(model_parms_file)
        if flag_marginals:
            self.marginals.read_model_marginals(model_marginals_file)
        if flag_xdata:
            self.generate_xdata()

    def __str__(self):
        return ".xdata" + str(self.get_events_nicknames_list())

    # TODO: finish this method to load model with default installed igor.
    @classmethod
    def load_default(cls, IgorSpecie, IgorChain, modelpath=rcParams['paths.igor_models']):
        """
        :return IgorModel loaded with the default location for specie and chain
        """        
        # IGoR run parameters
        #IgorSpecie    = specie #"mouse"
        #IgorChain     = chain #"tcr_beta"
        IgorModelPath = modelpath+"/"+IgorSpecie+"/"+IgorChain+"/"

        # FIXME: FIND A WAY TO GENERALIZE THIS WITH SOMEKIND OF STANDARD NAME
        flnModelParms = IgorModelPath + "models/model_parms.txt"
        flnModelMargs = IgorModelPath + "models/model_marginals.txt"
        print("Parms filename: ", flnModelParms)
        print("Margs filename: ", flnModelMargs)
        print("-"*50)

#        IgorRefGenomePath = IgorModelPath+"ref_genome/"
#        flnVGeneTemplate = IgorRefGenomePath+"genomicVs.fasta"
#        flnDGeneTemplate = IgorRefGenomePath+"genomicDs.fasta"
#        flnJGeneTemplate = IgorRefGenomePath+"genomicJs.fasta"
#        
#        flnVGeneCDR3Anchors = IgorRefGenomePath+"V_gene_CDR3_anchors.csv"
#        flnJGeneCDR3Anchors = IgorRefGenomePath+"J_gene_CDR3_anchors.csv"
        cls = IgorModel(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
        cls.specie = IgorSpecie
        cls.chain = IgorChain

        return cls

    @classmethod
    def load_from_networkx(cls, IgorSpecie, IgorChain):
        """
        :return IgorModel loaded with the default location for specie and chain
        """
        cls = IgorModel(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
        return cls

    def generate_xdata(self):
        Event_Genechoice_List = ['v_choice', 'j_choice', 'd_gene']
        Event_Dinucl_List = ['vd_dinucl', 'dj_dinucl', 'vj_dinucl']
        Event_Insertion_List = ['vd_ins', 'dj_ins', 'vj_ins']
        Event_Deletion_List = ['v_3_del', 'j_5_del', 'd_3_del', 'd_5_del']
        
        for key in self.marginals.marginals_dict:
            self.xdata[key] = xr.DataArray(self.marginals.marginals_dict[key], \
                          dims=tuple(self.marginals.network_dict[key]))
            #print "key: ", key, self.xdata[key].dims
            strCoord = "priority"
            self.xdata[key][strCoord] = self.parms.get_Event(key).priority
            for strDim in self.xdata[key].dims:
                self.xdata[key][strDim] = range(len(self.xdata[key][strDim]))
                if strDim in Event_Genechoice_List:
                    #print strDim
                    labels = self.parms.Event_dict[strDim]['name'].map(genLabel).values
                    #print type(labels)
                    strCoord = 'lbl__'+strDim
                    self.xdata[key][strCoord] = (strDim, labels) # range(len(self.xdata[key][coord]))
                elif not (strDim in Event_Dinucl_List):
                    labels = self.parms.Event_dict[strDim]['value'].values
                    #print strDim
                    #print labels
                    strCoord = 'lbl__'+strDim
                    self.xdata[key][strCoord] = (strDim, labels) # range(len(self.xdata[key][coord]))

    def get_Event_Marginal(self, event_nickname: str):
        """Returns an xarray with the marginal probability of the event given the nickname"""
        # FIXME: add new way to make the recursion.
        if event_nickname in self.parms.get_EventsNickname_list():
            da_event = self.xdata[event_nickname]
            dependencies = self.parms.Edges_dict[event_nickname]
            #1. Sort the dependencies by priority, then by dependencie
            if queue is not empty:
                self.get_Event_Marginal(nicki)
            return da_event
        else:
            print("Event nickname : " + event_nickname + " is not an event in this IGoR model.")
            return list()

    def plot_Event_Marginal(self, event_nickname:str):
        event = self.parms.get_Event(event_nickname, by_nickname=True)
        lblEvent = event_nickname.replace("_", " ")
        xEtiqueta = lblEvent
        yEtiqueta = "P"

        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.set_xlabel(xEtiqueta)
        ax.set_ylabel(yEtiqueta)

        if event.event_type == 'GeneChoice' :
            # TODO: plot using bar plot, and add the names on that.
            XX = self.xdata[event_nickname][event_nickname].values
            YY = self.xdata[event_nickname].values
            ax.bar(XX, YY)
            #self.xdata[event_nickname].plot(ax=ax)
            #ax.legend()
            # mdl.xdata[strEvent]['lbl__'+strEvent]
            ax.set_xticks(self.xdata[event_nickname][event_nickname].values)
            ax.set_xticklabels(self.xdata[event_nickname]['lbl__' + event_nickname].values, rotation=90)
            fig.tight_layout()

        elif event.event_type == 'Insertion':
            # Use labels as a coordinate.
            # Insertions are in principle independent,
            # FIXME: but if not what to do.
            YY = self.xdata[event_nickname].values
            XX = self.xdata[event_nickname]['lbl__' + event_nickname].values
            ax.plot(XX, YY)
        elif event.event_type == 'Deletion':
            #YY = self.xdata[event_nickname].values
            #XX = self.xdata[event_nickname]['lbl__' + event_nickname].values
            #ax.plot(XX, YY)
            self.xdata[event_nickname].plot(ax=ax)

        else:
            print(event)
            if len (self.xdata[event_nickname].shape ) == 2:
                self.xdata[event_nickname].plot(ax=ax)
                ax.legend()
                ax.set_xlabel(xEtiqueta)
                ax.set_ylabel(yEtiqueta)
                # mdl.xdata[strEvent]['lbl__'+strEvent]
                ax.set_xticks(self.xdata[event_nickname][event_nickname].values)
                ax.set_xticklabels(self.xdata[event_nickname]['lbl__' + event_nickname].values, rotation=90)
                ax.set_xticks(self.xdata[event_nickname][event_nickname].values)
                ax.set_xticklabels(self.xdata[event_nickname]['lbl__' + event_nickname].values, rotation=90)
                print("Accepted Events nicknames are : "+str(self.get_events_nicknames_list()))
        #return self.get_Event_Marginal(nickname)

    def get_events_types_list(self):
        "Return list of event types in current model"
        # The event list should be extracted from the Event_list
        events_set = set()
        for event in self.parms.Event_list:
            events_set.add(event.event_type)
        return list(events_set)

    def get_events_nicknames_list(self):
        "Return list of event nicknames in current model"
        # The event list should be extracted from the Event_list
        events_set = set()
        for event in self.parms.Event_list:
            events_set.add(event.nickname)
        return list(events_set)

    # def infer(self, batchname=None, iterations=5):
    #     import subprocess
    #     igor_exec = rcParams['paths.igor_exec']
    #     wd = "."
    #     cmd = igor_exec +" -set_wd " + wd + " -set_custom_model " + self.parms.model_parms_file + " -infer --N_iter "+str(iterations)
    #     print(cmd)

    def export_event_to_csv(self, strEvent, *args, **kargs):
        # if path_or_buf is None:
        #    path_or_buf = 'event__'+strEvent+".csv"
        # strEvent = 'd_3_del'
        df = self.xdata[strEvent].to_dataframe(name="Prob").drop('priority', 1)
        df.to_csv(*args, **kargs)

    def plot_dumm_report(self, strEvent):
        # strEvent = 'd_gene'
        import matplotlib.pyplot as plt
        fig, ax =  plt.subplots()
        dependencias = list(self.xdata[strEvent].dims)
        dependencias.remove(strEvent)
        dependencias_dim = [self.xdata[strEvent][dep].shape[0] for dep in dependencias]
        # eventos = eventos.remove(strEvent)
        lista = list()
        import numpy as np
        # np.ndindex()
        for index in np.ndindex(*dependencias_dim):
            dictionary = dict(zip(dependencias, index))
            # TODO: PLOT EACH DAMM FIGURE
            self.xdata[strEvent][dictionary].plot()
            aaa = [str(key) + "__" + str(dictionary[key]) for key in dictionary.keys()]
            lbl_file = "___".join(aaa)
            df = self.xdata[strEvent][dictionary].to_dataframe("P").drop('priority', 1)
            df.plot.bar(x="lbl__"+strEvent, y='P', ax=ax)
            print("*"*10)
            print(lbl_file)
            print(df)
            #df.to_csv(lbl_file+".csv")
            #fig.savefig(lbl_file+".png")
            #ax.clear()
        return fig

class IgorModel_Parms:
    """
    Class to get a list of Events directly from the *_parms.txt
    :param model_parms_file: Igor parms file path.
    """
    def __init__(self, model_parms_file=None):
        ## Parms file representation
        self.Event_list = list() # list of Rec_event
        self.Edges      = list()
        self.ErrorRate  = list()

        ## pygor definitions
        self.Event_dict = dict()
        self.Edges_dict = dict()
        self.dictNameNickname = dict()
        self.dictNicknameName = dict()
        self.G = nx.DiGraph()
        self.preMarginalDF = pd.DataFrame()
        self.model_parms_file = ""

        if model_parms_file is not None:
            print(model_parms_file)
            self.read_model_parms(model_parms_file)
            #self.get_EventDict_DataFrame()

    def __str__(self):
        tmpstr = "{ 'len Event_list': "+str(len(self.Event_list)) \
                +", 'len Egdes': "+str(len(self.Edges)) \
                +", 'len ErrorRate': "+str(len(self.ErrorRate))+" }"
        return tmpstr
        #return "{ Event_list, Egdes, ErrorRate}"

#    def __eq__(self, other):
#        if isinstance(self, other.__class__):
#            return (self.Event_list == other.Event_list) and (self.Edges == other.Edges) and (self.ErrorRate == other.ErrorRate)
#        else:
#            return NotImplemented
#
#    def __hash__(self):
#        return 3;


    @classmethod
    def from_network_dict(cls, network_dict:dict):
        # outfile << event_type<< ";" <<
        # SingleErrorRate
        cls = IgorModel_Parms()
        # 1. Create Event_list
        cls.Event_list = list()
        for nickname in network_dict.keys():
            # FIXME: Use default values
            # Create a default event by nickname
            dict_IgorRec_Event = IgorRec_Event_default_dict[nickname]
            event = IgorRec_Event.from_dict(dict_IgorRec_Event)
            try:
                cls.Event_list.append(event)
                print("New event has been added: ", dict_IgorRec_Event)
            except Exception as e:
                raise e
        # 2. Fill Events with realizations
        # 2.1. if event.event_type == 'GeneChoice':
        #       load file
        #mdl0.parms.Event_list[0].realizations[0])

        # load events from default dictionary.
        #
        return cls

    def load_events_from_dict(self, dicto):
        print(dicto)

    # TODO: Check how the imgt functions return data
    def load_GeneChoice_realizations_by_nickname(self, event_nickname:str, flnGenomic):
        event = self.get_Event(event_nickname)
        from Bio import SeqIO

        if event.event_type == 'GeneChoice':
            for index, record in enumerate(SeqIO.parse(flnGenomic, "fasta")):
                event_realization = IgorEvent_realization()
                event_realization.index = index
                event_realization.value = record.seq
                event_realization.name = record.description
                event.add_realization(IgorEvent_realization)
            print(event_nickname, " from file : ", flnGenomic)

    def load_Deletion_realizations_by_nickname(self, event_nickname: str, limits=(-4, 20)):
        event = self.get_Event(event_nickname)
        if event.event_type == 'Deletion':
            start, end = limits
            for index, ndels in enumerate(range(start, end)):
                event_realization = IgorEvent_realization()
                event_realization.index = index
                event_realization.value = ndels
            print(event_nickname, " limits : ", limits)

    def load_Insertion_realizations_by_nickname(self, event_nickname: str, limits=(0, 24)):
        event = self.get_Event(event_nickname)
        if event.event_type == 'Insertion':
            start, end = limits
            # FIXME: VALIDATE FOR POSITIVE VALUES
            for index, nins in enumerate(range(start, end)):
                event_realization = IgorEvent_realization()
                event_realization.index = index
                event_realization.value = nins
            print(event_nickname, " limits : ", limits)

    def load_DinucMarkov_realizations_by_nickname(self, event_nickname: str):
        event = self.get_Event(event_nickname)
        if event.event_type == 'DinucMarkov':
            for index, nt_char in enumerate(['A', 'C', 'G', 'T']):
                event_realization = IgorEvent_realization()
                event_realization.index = index
                event_realization.value = nt_char


    # FIXME: FINISH THIS METHOD
    def write_model_parms(self, filename="tmp_mdl_parms.txt"):
        """Writes a model graph structure from a model params object.
        Note that for now this method does not read the error rate information.
        """
        
        # Sort events in list with the your specific preference.
        # FIXME: FIND ANOTHER WAY TO WRITE IN A CORRECT ORDER
        #igor_nickname_list = ["v_choice", "j_choice", "d_gene", "v_3_del"]
        #self.get_Event(nicknameList)
        #self.Event_list
        strSepChar=";"
        with open(filename, "w") as ofile:
            #1. Write events
            ofile.write("@Event_list\n")
            #for event in self.Event_list:
            for nickname in Igor_nickname_list: # Igor_nicknameList is in IgorDefaults.py
                try:
                    event = self.get_Event(nickname)
                    strLine = "#" + \
                    str(event.event_type) + strSepChar + \
                    str(event.seq_type) + strSepChar + \
                    str(event.seq_side) + strSepChar + \
                    str(event.priority) + strSepChar + \
                    str(event.nickname) + "\n"
                    ofile.write(strLine)
                    # WRITE THE LIST OF REALIZATIONS adding character '%'
                    df = event.get_realization_DataFrame()
                    str_df = df.to_csv(sep=strSepChar, header=False)
                    str_realization_list=""
                    for strLine in str_df.split("\n"):
                        str_realization_list = str_realization_list + "%"+strLine + "\n"
                    ofile.write(str_realization_list)

                # FIXME: NOT GOOD TO MAKE A BLIND PASS BUT A TEMPORAL SOLUTION UNTIL I FIND A WAY TO ESTABLISH PRIORITIES ON THE EVENTS
                except Exception as e:
                    pass

            #2. Write Edges
            ofile.write("@Edges\n")
            
            #3. Write ErrorRate
            ofile.write("@ErrorRate\n")
    
    def read_model_parms(self, filename):
        """Reads a model graph structure from a model params file.
        Note that for now this method does not read the error rate information.
        """
        with open(filename, "r") as ofile:
            # dictionary containing recarrays?
            line = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)

            if strip_line == "@Event_list":
                self.read_Event_list(ofile)

            line = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
            if strip_line == "@Edges":
                self.read_Edges(ofile)

            # FIXME: ErrorRate added
            line = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
            if strip_line == "@ErrorRate" :
                self.read_ErrorRate(ofile)
        self.model_parms_file = filename
        self.get_EventDict_DataFrame()

    # save in Event_list
    def read_Event_list(self, ofile):
        lastPos    = ofile.tell()
        line = ofile.readline()
        strip_line = line.rstrip('\n')  # Remove end of line character
        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        #event = Rec_Event()

        while strip_line[0] == '#':
            # get the metadata of the event list
            event_metadata = strip_line[1:].split(";") #GeneChoice;V_gene;Undefined_side;7;v_choice
            event_metadata[3] = int(event_metadata[3]) # change priority to integer
            event = IgorRec_Event(*event_metadata)
            #self.G.add_node(event.nickname)
            # Now read the realizations (or possibilities)
            lastPos    = ofile.tell()
            line       = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)

            while strip_line[0] == '%':
                realization = IgorEvent_realization()
                realizData  = strip_line[1:].split(";")
                if event.event_type == "GeneChoice":
                    realization.name  = realizData[0]
                    realization.value = realizData[1]
                    realization.index = int(realizData[2])
                elif event.event_type == "DinucMarkov":
                    realization.value = realizData[0]
                    realization.index = int(realizData[1])
                else:
                    realization.value = int(realizData[0])
                    realization.index = int(realizData[1])

                event.add_realization(realization)
                # next line
                lastPos  = ofile.tell()
                line     = ofile.readline()
                strip_line = line.rstrip('\n')  # Remove end of line character
                strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)

            self.Event_list.append(event)
        ofile.seek(lastPos)

    def read_Edges(self, ofile):
        #print "read_Edges"
        lastPos  = ofile.tell()
        line = ofile.readline()
        strip_line = line.rstrip('\n')  # Remove end of line character
        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        while strip_line[0] == '%':
            edge = strip_line[1:].split(';')
            self.Edges.append(edge)
            #read nextline
            lastPos  = ofile.tell()
            line = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        ofile.seek(lastPos)

    def read_ErrorRate(self, ofile):
        lastPos  = ofile.tell()
        line = ofile.readline()
        strip_line = line.rstrip('\n')  # Remove end of line character
        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        while strip_line[0] == '#':
            if 'SingleErrorRate' == strip_line[1:] :
                lastPos  = ofile.tell()
                line = ofile.readline()
                strip_line = line.rstrip('\n').rstrip()
                error = strip_line
                self.ErrorRate = {"SingleErrorRate" : error }
        ofile.seek(lastPos)

    def get_EventsNickname_list(self):
        return [event.nickname for event in self.Event_list]

    def get_EventsName_list(self):
        return [event.name for event in self.Event_list]

    def get_Event(self, event_nickname, by_nickname=True):
        """Returns the RecEvent with corresponding name or nickname."""
        if by_nickname:
            for ev in self.Event_list:
                if ev.nickname == event_nickname:
                    return ev
            raise Exception(
                'RecEvent with nickname \"' + event_nickname + "\" not found.")
        else:
            for ev in self.Event_list:
                if ev.name == event_nickname:
                    return ev
            raise Exception(
                'RecEvent with name \"' + event_nickname + "\" not found.")

    def get_EventDict_DataFrame(self):
        self.Event_dict = dict()
        self.dictNameNickname = dict()
        #dictio = dict()
        for event in self.Event_list:
            #dictio[event.nickname] = event.get_realization_DataFrame()
            self.Event_dict[event.nickname] = event.get_realization_DataFrame()
            self.dictNameNickname[event.name] = event.nickname
        #return dictio
        self.dictNicknameName = {v: k for k, v in self.dictNameNickname.items()}
        self.getBayesGraph()
    
    def getBayesGraph(self):
        self.G = nx.DiGraph()
        for rec_event in self.Event_list:
            self.G.add_node(rec_event.nickname)
            self.Edges_dict[rec_event.nickname] = list()
            self.dictNameNickname[rec_event.name] = rec_event.nickname
        
        for edge in self.Edges:
            # Graph to get the dependecies
            self.G.add_edge(self.dictNameNickname[edge[0]], self.dictNameNickname[edge[1]])
            self.Edges_dict[self.dictNameNickname[edge[1]]].append(self.dictNameNickname[edge[0]])
        #self.G = self.G.reverse()
        
    def genPreMarginalDF(self):
        data=[]
        for event in self.Event_list:    
            #print (parms.dictNameNickname[event.name])
            #parms.Edges
            #tmpDict = dict()
            lista = []
            for edge in self.Edges:
                if edge[1] == event.name:
                    #print(parms.dictNameNickname[edge[0]])
                    #print(edge[0])
                    lista.append(self.dictNameNickname[edge[0]])
            tmpDict = {'event': event.nickname, 'priority': event.priority, 'Edges': lista }
            data.append(tmpDict)
        
        self.preMarginalDF = pd.DataFrame(data) #.set_index('event')
        self.preMarginalDF['nEdges'] = self.preMarginalDF['Edges'].map(len)
        self.preMarginalDF.sort_values(['priority', 'nEdges'], ascending=[False,True])

    def genMarginalFile(self, model_marginals_file=None):
        self.genPreMarginalDF()
        #self.preMarginalDF
        if model_marginals_file == None:
            model_marginals_file = "model_marginals.txt"
        ofile = open(model_marginals_file, "w")
        for index, row in self.preMarginalDF.iterrows():
            nickname = row['event']
            ofile.write("@"+nickname+"\n")
            #DimEvent = len(parms.Event_dict[event.nickname])
            #DimEdges = len(parms.Edges_dict[event.nickname])
            
            DimEvent = len(self.Event_dict[nickname])
            strDimLine = "$Dim["
            DimList = []
            if row['nEdges'] == 0:
                strDimLine = strDimLine +str(DimEvent)
                strDimLine = strDimLine +"]"
            else:
                for evNick in row['Edges']: #parms.Edges_dict[event.nickname]:
                    Dim = len(self.Event_dict[evNick])
                    strDimLine = strDimLine +str(Dim)+","
                    DimList.append(Dim)
                strDimLine = strDimLine + str(DimEvent)
                strDimLine = strDimLine +"]"
            ofile.write(strDimLine+"\n")
            
            
            lista = row['Edges'] # self.Event_dict[nickname]
            for indices in np.ndindex(tuple(DimList)):
                #print indices
                strTmp = "#"
                for ii in range(len(lista)):
                    strTmp = strTmp+"["+lista[ii]+","+str(indices[ii])+"]"
                    if not (ii == len(lista)-1):
                        strTmp = strTmp + ","
                ofile.write(strTmp+"\n")
                ofile.write("%")
                unifProb = (1./DimEvent)
                for jj in range(DimEvent):
                    ofile.write(str(unifProb))
                    if not (jj == DimEvent-1):
                        ofile.write(",")
                ofile.write("\n")
                               
        ofile.close()

    def plot_Graph(self): # FIXME: ALLOW the possibility to pass an ax like ax=None):
        """Return a plot of the bayesian network """
        #if ax is None:

        pos = nx.spring_layout(self.G)
        # priorities up
        prio_dict = dict()
        for event in self.Event_list:
            if not (event.priority in prio_dict):
                prio_dict[event.priority] = list()
            prio_dict[event.priority].append(event)
        #print(str(prio_dict))
        xwidth = 240
        yfactor = 40
        for key in prio_dict:
            lenKey = len(prio_dict[key])
            if lenKey == 1:
                pos[prio_dict[key][0].nickname] = np.array([float(xwidth) / 2.0, float(key) * yfactor])
            else:
                xx = np.linspace(0, xwidth, lenKey)
                for ii, ev in enumerate(prio_dict[key]):
                    xpos = xx[ii]  # float(xwidth)*float(ii)/float(lenKey)
                    pos[ev.nickname] = np.array([xpos, float(key) * yfactor])

        #### import matplotlib.pyplot as plt
        #### fig, ax = plt.subplots()
        #### ax.set_aspect('equal')
        #### nx.draw(self.G, pos=pos, ax=ax, with_labels=True, arrows=True, arrowsize=20,
        ####         node_size=800, font_size=10, font_weight='bold')  # FIXME: make a better plot: cutting edges.

        #### plt.show()

        import hvplot.networkx as hvnx
        graph = hvnx.draw(self.G, with_labels=True, FontSize=10, pos=pos, alpha=0.5,
                        arrowstyle='fancy', arrowsize=2000, node_size=1000, width=400, height=400)
                        ##, arrows=True, arrowsize=20, node_size=800, font_size=10, font_weight='bold')
        return graph

class IgorRec_Event:
    """Recombination event class containing event's name, type, realizations,
    etc... Similar to IGoR's C++ RecEvent class.
    """
    def __init__(self, event_type, seq_type, seq_side, priority,
                 nickname):
        self.event_type = event_type
        self.seq_type = seq_type
        self.seq_side = seq_side
        self.priority = priority
        self.realizations = list()
        self.name = ""
        self.nickname = nickname
#        if nickname is not None:
#            self.nickname = nickname
        self.update_name()

    def to_dict(self):
        dictIgorRec_Event = {
            "event_type": self.event_type, \
            "seq_type": self.seq_type, \
            "seq_side": self.seq_side, \
            "priority": self.priority, \
            "realizations": self.realizations, \
            "name": self.name, \
            "nickname": self.nickname
        }

        return dictIgorRec_Event


    @classmethod
    def from_dict(cls, dict_IgorRec_Event:dict):
        """Returns a IgorRec_Event based on dictionary
        """
        cls = IgorRec_Event(dict_IgorRec_Event["event_type"], dict_IgorRec_Event["seq_type"],
                            dict_IgorRec_Event["seq_side"], dict_IgorRec_Event["priority"], dict_IgorRec_Event["nickname"])
        # 'event_type', 'seq_type', 'seq_side', 'priority', and 'nickname'
        # FIXME: Is better to make this class as an extension of a dictionary container?
        # Given the nickname complete the events with
        #cls.nickname = dict_IgorRec_Event["nickname"]
        #cls.event_type = dict_IgorRec_Event["event_type"]
        #cls.seq_type = dict_IgorRec_Event["seq_type"]
        #cls.seq_side = dict_IgorRec_Event["seq_side"]
        #cls.priority = dict_IgorRec_Event["priority"]
        cls.realizations = dict_IgorRec_Event["realizations"] # TODO: CREATE FUNCTION TO GENERATE realizations vector
        #cls.name = dict_IgorRec_Event["name"]

        return cls

    @classmethod
    def from_default_nickname(cls, nickname:str):
        cls = IgorRec_Event.to_dict(IgorRec_Event_default_dict[nickname])
        return cls

    def __str__(self):
        return str(self.to_dict())

    def __lt__(self, other):
        return self.priority < other.priority 
#        if ( self.priority < other.priority ):
#            return True
#        elif ( self.priority == other.priority  ):
#            # FIXME: less dependencies should be on top

    def add_realization(self, realization):
        """Add a realization to the RecEvent realizations list."""
        self.realizations.append(realization)
        self.update_name()

    def update_name(self):
        """Updates the name of the event (will have no effect if the RecEvent
        has not been modified since the last call).

        """
        if self.event_type == "DinucMarkov":
            self.name = self.event_type + "_" + self.seq_type + "_" + \
                        self.seq_side + "_prio" + \
                        str(self.priority) + "_size" + \
                        str(len(self.realizations) ** 2)
        else:
            self.name = self.event_type + "_" + self.seq_type + "_" + \
                        self.seq_side + "_prio" + \
                        str(self.priority) + "_size" + \
                        str(len(self.realizations))

    # TODO: Create a realization vector from a fasta file
    def set_realization_vector(self):
        if self.event_type == 'GeneChoice':
            print('GeneChoice')

    def set_realization_vector_GeneChoice(self, flnGenomic:str):
        """
        Sets a realization vector from a filename
        :param flnGenomic: fasta file with the genomic template IMGT or other template.
        """
        #FIXME: FINISH IT
        # TODO: Add realizations from fasta file.
        from Bio import SeqIO
        #for record in list(SeqIO.parse(flnGenomic, "fasta")):

    def get_realization_vector(self):
        """This methods returns the event realizations sorted by the
        realization index as a list.
        """
        if self.event_type == 'GeneChoice':
            tmp = [""] * len(self.realizations)  # empty(, dtype = str)
        else:
            tmp = np.empty(len(self.realizations),
                              dtype=type(self.realizations[0].value))
        # print("Unfinished method get realization vector")
        processed_real_indices = []
        for real in self.realizations:
            if processed_real_indices.count(real.index) == 0:
                if real.name != "":
                    tmp[real.index] = real.name
                else:
                    tmp[real.index] = real.value
                processed_real_indices.append(real.index)
            else:
                print("REALIZATION INDICES ARE DEGENERATE")

        return tmp


    def get_realization_DataFrame(self):
        """ return an Event realizations as a pandas DataFrame to manipulate it.

        """
        return pd.DataFrame.from_records([realiz.to_dict() for realiz in self.realizations], index='index').sort_index()

class IgorEvent_realization:
    """A small class storing for each RecEvent realization its name, value and
    corresponding index.

    """
    def __init__(self):
        self.name  = "" #name
        self.value = "" #value
        self.index = "" #index

    def __lt__(self, other):
        return self.index < other.index

    def __gt__(self, other):
        return self.index > other.index

    def __str__(self):
        if self.name == "":
            return self.value+";"+self.index
        else:
            return self.name+";"+self.value+";"+str(self.index)

    def __repr__(self):
        return "Event_realization(" + str(self.index) + ")"

    def to_dict(self):
        return {
            'index' : self.index,
            'value' : self.value,
            'name'  : self.name
            }

    @classmethod
    def from_dict(self, event_dict:dict):
        cls = IgorEvent_realization()
        cls.index = event_dict['index']
        cls.value = event_dict['value']
        cls.name = event_dict['name']
        return cls
#    def __eq__(self, other):
#        if isinstance(self, other.__class__):
#            return (    (self.event_type   == other.event_type) \
#                    and (self.seq_type     == other.seq_type) \
#                    and (self.seq_side     == other.seq_side) \
#                    and (self.priority     == other.priority) \
#                    and (self.realizations == other.realizations) \
#                    and (self.name         == other.name) \
#                    and (self.nickname     == other.nickname) \
#                    )
#        else:
#            return NotImplemented
#
#    def __hash__(self):
#        return hash((self.event_type, \
#                     self.seq_type, \
#                     self.seq_side, \
#                     self.priority, \
#                     self.realizations, \
#                     self.name, \
#                     self.nickname))

#    def __str__(self):
#        return self.event_type+";"+self.seq_type+";"+self.seq_side+";"+self.priority+";"+self.nickname
#
#    def __repr__(self):
#        return "Rec_event(" + self.nickname + ")"



### FIXME:
# @recombinationEvent
# $Dim
# #Indices of the realizations of the parent events
# %1d probability array.

class IgorModel_Marginals:
    """
    Class to get a list of Events directly from the *_parms.txt
    :param model_marginals_file: Igor marginals file.
    """
    def __init__(self, model_marginals_file=None):
#        self.Event_list = list() # list of Rec_event
#        self.Edges      = list()
#        self.error_rate = list()
        self.marginals_dict = {}
        self.network_dict = {}
        #self.G          = nx.DiGraph()
        self.model_marginals_file = ""
        if model_marginals_file is not None:
            self.read_model_marginals(model_marginals_file)

    #  @d_3_del
    #  $Dim[3,21,21]
    #  #[d_gene,0],[d_5_del,0]
    #  %0,0,0,1.6468e-08,0.00482319,1.08101e-09,0.0195311,0.0210679,0.0359338,0.0328678,2.25686e-05,4.97463e-07,0,9.31048e-08,1.01642e-05,0.000536761,0.0260845,0.0391021,0.319224,0.289631,0.211165
    #  #[d_gene,0],[d_5_del,1]
    #  %0,0,6.86291e-08,2.00464e-09,0.00163832,2.02919e-06,0.0306066,0.0126832,0.000872623,0.016518,0.00495292,0.000776747,4.45576e-05,0.000667902,0.00274004,0.00435049,0.300943,0.182499,0.13817,0.302534,0

    def read_model_marginals(self, filename, dim_names=False):
        """Reads a model marginals file. Returns a tuple containing a dict
        containing the individual events probabilities indexed by the events
        nicknames and a dict containing the list of dimension names/ordering for
        each event.
        """
        with open(filename, "r") as ofile:
            # Model parameters are stored inside a dictionnary of ndarrays
#            marginals_dict = {}
#            network_dict = {}
            element_name = ""
            first = True
            first_dim_line = False
            element_marginal_array = []
            indices_array = []

            for line in ofile:
                strip_line = line.rstrip("\n")  # Remove end of line character
                if strip_line[0] == "@":
                    first_dim_line = True
                    if not first:
                        # Add the previous to the dictionnary
                        self.marginals_dict[element_name] = element_marginal_array
                    else:
                        first = False

                    element_name = strip_line[1:]
                # print element_name

                if strip_line[0] == "$":
                    # define array dimensions
                    coma_index = strip_line.find(",")
                    dimensions = []

                    # Get rid of $Dim[
                    previous_coma_index = 4
                    while coma_index != -1:
                        dimensions.append(
                            int(strip_line[previous_coma_index + 1:coma_index]))
                        previous_coma_index = coma_index
                        coma_index = strip_line.find(",", coma_index + 1)

                    # Add last dimension and get rid of the closing bracket
                    dimensions.append(int(strip_line[previous_coma_index + 1:-1]))

                    element_marginal_array = np.ndarray(shape=dimensions)

                if strip_line[0] == "#":
                    if first_dim_line:
                        dimensions_names = []
                        if len(dimensions) > 1:
                            comma_index = strip_line.find(",")
                            opening_bracket_index = strip_line.find("[")
                            while opening_bracket_index != -1:
                                dimensions_names.append(
                                    strip_line[
                                        opening_bracket_index + 1:comma_index])
                                opening_bracket_index = strip_line.find(
                                    "[", comma_index)
                                comma_index = strip_line.find(
                                    ",", opening_bracket_index)
                        first_dim_line = False
                        dimensions_names.append(element_name)
                        self.network_dict[element_name] = dimensions_names

                    # update indices
                    indices_array = []
                    if len(dimensions) > 1:
                        comma_index = strip_line.find(",")
                        closing_brack_index = strip_line.find("]")
                        while closing_brack_index != -1:
                            indices_array.append(int(
                                strip_line[comma_index + 1:closing_brack_index]))
                            opening_bracket_index = strip_line.find(
                                "[", closing_brack_index)
                            comma_index = strip_line.find(
                                ",", opening_bracket_index)
                            closing_brack_index = strip_line.find(
                                "]", closing_brack_index + 1)

                if strip_line[0] == "%":
                    # read doubles
                    coma_index = strip_line.find(",")
                    marginals_values = []

                    # Get rid of the %
                    previous_coma_index = 0
                    while coma_index != -1:
                        marginals_values.append(
                            float(strip_line[previous_coma_index + 1:coma_index]))
                        previous_coma_index = coma_index
                        coma_index = strip_line.find(",", coma_index + 1)

                    # Add last dimension and get rid of the closing bracket
                    marginals_values.append(
                        float(strip_line[previous_coma_index + 1:]))
                    if len(marginals_values) != dimensions[-1]:
                        print("problem")
                    element_marginal_array[tuple(indices_array)] = marginals_values
            self.marginals_dict[element_name] = element_marginal_array

        self.model_marginals_file = filename

        #return marginals_dict, network_dict



    def write_model_marginals_from_parms(self, filename, model_parms_file=None, model_marginals_file=None):
        self.marginals_dict = {}
        self.network_dict = {}

### IGOR BEST SCENARIOS VDJ ###
class IgorBestScenariosVDJ:
    def __init__(self):
        self.seq_index = -1
        self.scenario_rank = -1
        self.scenario_proba_cond_seq = -1
        self.id_v_choice = -1
        self.id_j_choice = -1
        self.id_d_gene = -1
        self.id_v_3_del = -1
        self.id_d_5_del = -1
        self.id_d_3_del = -1
        self.id_j_5_del = -1
        self.id_vd_ins = -1
        self.vd_dinucl = list()
        self.id_dj_ins = -1
        self.dj_dinucl = list()
        self.mismatches = list()
        self.mismatcheslen = -1

        # Case of use of the IgorModel.Model_Parms class.
        self.flnModelParms = ""
        self.mdlParms = ""
        self.mdl = ""
        # self.IgorModelParms
        self.strSeq_index = ""

    #        # To fetch data from datase connection to a database
    #        self.IgorDB = ""

    def __str__(self):
        return str(self.to_dict())

    def setModel_Parms(self, flnModelParms):
        self.flnModelParms = flnModelParms
        self.mdlParms = IgorModel_Parms(model_parms_file=self.flnModelParms)

    def to_dict(self):
        dictBestScenario = {
            "seq_index": self.seq_index, \
            "scenario_rank": self.scenario_rank, \
            "scenario_proba_cond_seq": self.scenario_proba_cond_seq, \
            "v_choice": self.id_v_choice, \
            "j_choice": self.id_j_choice, \
            "d_gene": self.id_d_gene, \
            "v_3_del": self.id_v_3_del, \
            "d_5_del": self.id_d_5_del, \
            "d_3_del": self.id_d_3_del, \
            "j_5_del": self.id_j_5_del, \
            "vd_ins": self.id_vd_ins, \
            "vd_dinucl": self.vd_dinucl, \
            "dj_ins": self.id_dj_ins, \
            "dj_dinucl": self.dj_dinucl, \
            "mismatches": self.mismatches, \
            "mismatcheslen": self.mismatcheslen
        }

        return dictBestScenario

    def to_dict_names(self):
        dictBestScenario = {
            "seq_index": self.seq_index, \
            "scenario_rank": self.scenario_rank, \
            "scenario_proba_cond_seq": self.scenario_proba_cond_seq, \
            "v_choice": self.getV_gene_name(), \
            "j_choice": self.getJ_gene_name(), \
            "d_gene": self.getD_gene_name(), \
            "v_3_del": self.getV_3_dels(), \
            "d_5_del": self.getD_5_dels(), \
            "d_3_del": self.getD_3_dels(), \
            "j_5_del": self.getJ_5_dels(), \
            "vd_ins": self.getVD_ins(), \
            "vd_dinucl": self.vd_dinucl, \
            "dj_ins": self.getDJ_ins(), \
            "dj_dinucl": self.dj_dinucl, \
            "mismatches": self.mismatches, \
            "mismatcheslen": self.mismatcheslen
        }

        return dictBestScenario

    def to_dict_ntsequences(self):
        dictBestScenario = {
            "seq_index": self.seq_index, \
            "scenario_rank": self.scenario_rank, \
            "scenario_proba_cond_seq": self.scenario_proba_cond_seq, \
            "v_choice": self.getV_ntsequence(), \
            "j_choice": self.getJ_ntsequence(), \
            "d_gene": self.getD_ntsequence(), \
            "v_3_del": self.getV_3_dels(), \
            "d_5_del": self.getD_5_dels(), \
            "d_3_del": self.getD_3_dels(), \
            "j_5_del": self.getJ_5_dels(), \
            "vd_ins": self.getVD_ins(), \
            "vd_dinucl": self.getVD_Region(), \
            "dj_ins": self.getDJ_ins(), \
            "dj_dinucl": self.getDJ_Region(), \
            "mismatches": self.mismatches, \
            "mismatcheslen": self.mismatcheslen
        }

        return dictBestScenario

    @classmethod
    def load_FromLineBestScenario(cls, line, delimiter=";"):
        # seq_index;scenario_rank;scenario_proba_cond_seq;GeneChoice_V_gene_Undefined_side_prio7_size35;GeneChoice_J_gene_Undefined_side_prio7_size14;GeneChoice_D_gene_Undefined_side_prio6_size2;Deletion_V_gene_Three_prime_prio5_size21;Deletion_D_gene_Five_prime_prio5_size21;Deletion_D_gene_Three_prime_prio5_size21;Deletion_J_gene_Five_prime_prio5_size23;Insertion_VD_genes_Undefined_side_prio4_size31;DinucMarkov_VD_genes_Undefined_side_prio3_size16;Insertion_DJ_gene_Undefined_side_prio2_size31;DinucMarkov_DJ_gene_Undefined_side_prio1_size16;Mismatches
        cls = IgorBestScenariosVDJ()
        linesplit = line.split(delimiter)
        linesplit = line.split(";")
        for ii in range(len(linesplit)):
            # TODO: find a better way to do this, if is a list keep it as list
            if (ii in [11, 13, 14]):
                linesplit[ii] = linesplit[ii]
            else:
                linesplit[ii] = linesplit[ii].replace("(", "").replace(")", "")

        try:
            # 1;1;0.123596;(14);(9);(1);(4);(8);(7);(11);(0);();(9);(0,2,0,1,2,3,2,0,0);(122,123,124)
            cls.seq_index = int(linesplit[0])
            cls.scenario_rank = int(linesplit[1])
            cls.scenario_proba_cond_seq = float(linesplit[2])
            print(linesplit[3], type(linesplit[3]), len(linesplit[3]))
            cls.id_v_choice = int(linesplit[3])
            cls.id_j_choice = int(linesplit[4])
            cls.id_d_gene = int(linesplit[5])
            cls.id_v_3_del = int(linesplit[6])
            cls.id_d_5_del = int(linesplit[7])
            cls.id_d_3_del = int(linesplit[8])
            cls.id_j_5_del = int(linesplit[9])
            cls.id_vd_ins = int(linesplit[10])
            cls.vd_dinucl = eval(linesplit[11])
            cls.id_dj_ins = int(linesplit[12])
            cls.dj_dinucl = eval(linesplit[13])
            cls.mismatches = eval(linesplit[14])
            cls.mismatcheslen = int(len(cls.mismatches))

            return cls

        except Exception as e:
            print(e)
            raise e

    @classmethod
    def load_FromDict(cls, dictBestScenarios):
        """
        Return a IgorBestScenariosVDJ instance from a IgorSqlRecord.
        :param sqlRecordAlign: record of a sql database table.
        :param strGene_name: gene_name associated to the record.
        :return: IgorAlignment_data instance
        """
        cls = IgorBestScenariosVDJ()
        try:
            #            cls.seq_index      = int(sqlRecordBestScenarios[ 0])
            #            cls.scenario_rank  = int(sqlRecordBestScenarios[ 1])
            #            cls.scenario_proba_cond_seq = float(sqlRecordBestScenarios[2])
            cls.id_v_choice = int(dictBestScenarios["v_choice"])
            cls.id_j_choice = int(dictBestScenarios["j_choice"])
            cls.id_d_gene = int(dictBestScenarios["d_gene"])
            cls.id_v_3_del = int(dictBestScenarios["v_3_del"])
            cls.id_d_5_del = int(dictBestScenarios["d_5_del"])
            cls.id_d_3_del = int(dictBestScenarios["d_3_del"])
            cls.id_j_5_del = int(dictBestScenarios["j_5_del"])
            cls.id_vd_ins = int(dictBestScenarios["vd_ins"])
            cls.vd_dinucl = eval(dictBestScenarios["vd_dinucl"])
            cls.id_dj_ins = int(dictBestScenarios["dj_ins"])
            cls.dj_dinucl = eval(dictBestScenarios["dj_dinucl"])
            #            cls.mismatches     = eval(dictBestScenarios["mismatches"])
            #            cls.mismatcheslen  = int(len(cls.mismatches) )

            return cls
        except Exception as e:
            print(e)
            raise e

    @classmethod
    def load_FromSQLRecord(cls, sqlRecordBestScenarios):
        """
        Return a IgorBestScenariosVDJ instance from a IgorSqlRecord.
        :param sqlRecordAlign: record of a sql database table.
        :param strGene_name: gene_name associated to the record.
        :return: IgorAlignment_data instance
        """
        cls = IgorBestScenariosVDJ()
        try:
            cls.seq_index = int(sqlRecordBestScenarios[0])
            cls.scenario_rank = int(sqlRecordBestScenarios[1])
            cls.scenario_proba_cond_seq = float(sqlRecordBestScenarios[2])
            cls.id_v_choice = int(sqlRecordBestScenarios[3])
            cls.id_j_choice = int(sqlRecordBestScenarios[4])
            cls.id_d_gene = int(sqlRecordBestScenarios[5])
            cls.id_v_3_del = int(sqlRecordBestScenarios[6])
            cls.id_d_5_del = int(sqlRecordBestScenarios[7])
            cls.id_d_3_del = int(sqlRecordBestScenarios[8])
            cls.id_j_5_del = int(sqlRecordBestScenarios[9])
            cls.id_vd_ins = int(sqlRecordBestScenarios[10])
            cls.vd_dinucl = eval(sqlRecordBestScenarios[11])
            cls.id_dj_ins = int(sqlRecordBestScenarios[12])
            cls.dj_dinucl = eval(sqlRecordBestScenarios[13])
            cls.mismatches = eval(sqlRecordBestScenarios[14])
            cls.mismatcheslen = int(sqlRecordBestScenarios[15])

            return cls
        except Exception as e:
            print(e)
            raise e

    # TODO: finish this class
    @classmethod
    def load_FromEventNameValues(cls, mdl, seq_index, strSeq_index, scenario_dict):
        # v_3_del = 2
        # d_5_del = 6
        # d_3_del = 1
        # vd_ins  = 1 and should be a "C"
        # dj_ins  = 3 and should be a "TCT"
        # FIXME: CORRECT ERROR MESSAGES and ADD documentation and VALIDATIONS OF EVENTS
        """
        Return a IgorBestScenariosVDJ instance from a dict of names or values.
        :param strGene_name: gene_name associated to the record.
        :return: IgorAlignment_data instance
        """
        ##### The event I think is the best one
        cls = IgorBestScenariosVDJ()  # .load_FromSQLRecord(record_bs[ii])
        # cls.setModel_Parms(flnModelParms)
        cls.mdl = mdl
        cls.seq_index = seq_index  # 59
        cls.strSeq_index = strSeq_index  # db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]
        Event_GeneChoice = ['v_choice', 'j_choice', 'd_gene']
        Event_Deletions = ['v_3_del', 'd_5_del', 'd_3_del', 'j_5_del']
        Event_Insertions = ['vd_ins', 'dj_ins']
        Event_Dinucl = ['vd_dinucl', 'dj_dinucl']
        for event_nickname in scenario_dict.keys():
            if event_nickname in Event_GeneChoice:
                pd_event = cls.mdl.parms.Event_dict[event_nickname]
                gene_name = scenario_dict[event_nickname]  # 'TRBV17*01'
                gene_id = pd_event.loc[pd_event['name'] == gene_name].index.values[0]
                if event_nickname == 'v_choice':
                    cls.id_v_choice = gene_id
                elif event_nickname == 'j_choice':
                    cls.id_j_choice = gene_id
                elif event_nickname == 'd_gene':
                    cls.id_d_gene = gene_id
                else:
                    print("Something vey bad happen with " + str(scenario_dict[event_nickname]))

            elif event_nickname in Event_Deletions:
                pd_event = cls.mdl.parms.Event_dict[event_nickname]
                realiz_name = scenario_dict[event_nickname]
                realiz_id = pd_event.loc[pd_event['value'] == realiz_name].index.values[0]
                if event_nickname == 'v_3_del':
                    cls.id_v_3_del = realiz_id
                elif event_nickname == 'd_5_del':
                    cls.id_d_5_del = realiz_id
                elif event_nickname == 'd_3_del':
                    cls.id_d_3_del = realiz_id
                elif event_nickname == 'j_5_del':
                    cls.id_j_5_del = realiz_id
                else:
                    print("Something vey bad happen with " + str(scenario_dict[event_nickname]))

            elif event_nickname in Event_Insertions:
                pd_event = cls.mdl.parms.Event_dict[event_nickname]
                realiz_name = scenario_dict[event_nickname]
                realiz_id = pd_event.loc[pd_event['value'] == realiz_name].index.values[0]
                if event_nickname == 'v_3_del':
                    cls.id_vd_ins = realiz_id
                elif event_nickname == 'd_5_del':
                    cls.id_dj_ins = realiz_id
                else:
                    print("Something vey bad happen with " + str(scenario_dict[event_nickname]))

            elif event_nickname in Event_Dinucl:
                if event_nickname == 'vd_dinucl':
                    pd_event = cls.mdl.parms.Event_dict[event_nickname]
                    str_sequence = scenario_dict[event_nickname]
                    list_id_seq = list()
                    for str_nt in str_sequence:
                        realiz_name = str_nt
                        realiz_id = pd_event.loc[pd_event['value'] == realiz_name].index.values[0]
                        list_id_seq.append(realiz_id)
                    cls.vd_dinucl = list_id_seq

                elif event_nickname == 'dj_dinucl':
                    pd_event = cls.mdl.parms.Event_dict[event_nickname]
                    str_sequence = scenario_dict[event_nickname]
                    list_id_seq = list()
                    for str_nt in str_sequence:
                        realiz_name = str_nt
                        realiz_id = pd_event.loc[pd_event['value'] == realiz_name].index.values[0]
                        list_id_seq.append(realiz_id)
                    cls.dj_dinucl = list_id_seq

                else:
                    print("Something wrong with " + str(event_nickname))

            elif event_nickname in ['mismatches']:
                if isinstance(scenario_dict[event_nickname], list):
                    cls.mismatches = scenario_dict[event_nickname]
                    cls.mismatcheslen = len(cls.mismatches)
                else:
                    print("mismatches aren't list")
            else:
                print("Something bad happen!")

        return cls

    def save_scenario_fasta(self, outfilename):
        ofileScen = open(outfilename, "w")
        ofileScen.write(self.str_scenario_fasta())
        #        ofileScen.write("> "+str(self.seq_index)+", rank: "+str(self.scenario_rank)+ ", prob: "+str(self.scenario_proba_cond_seq)+"\n")
        #        ofileScen.write( self.strSeq_index  + "\n" )
        #        ofileScen.write( self.getV_fasta()  + "\n" )
        #        ofileScen.write( self.getVD_fasta() + "\n" )
        #        ofileScen.write( self.getD_fasta()  + "\n" )
        #        ofileScen.write( self.getDJ_fasta() + "\n" )
        #        ofileScen.write( self.getJ_fasta()  + "\n" )
        ofileScen.close()

    def str_scenario_fasta(self):
        strScenarioFasta = ""
        strScenarioFasta = strScenarioFasta + ">" + str(self.seq_index) + ", rank: " + str(
            self.scenario_rank) + ", prob: " + str(self.scenario_proba_cond_seq) + "\n"
        strScenarioFasta = strScenarioFasta + self.strSeq_index + "\n"
        strScenarioFasta = strScenarioFasta + self.getV_fasta() + "\n"
        strScenarioFasta = strScenarioFasta + self.getVD_fasta() + "\n"
        strScenarioFasta = strScenarioFasta + self.getD_fasta() + "\n"
        strScenarioFasta = strScenarioFasta + self.getDJ_fasta() + "\n"
        strScenarioFasta = strScenarioFasta + self.getJ_fasta() + "\n"

        # FIXME: TEMPORAL
        strScenarioFasta = strScenarioFasta + "> v_choice\n"
        strScenarioFasta = strScenarioFasta + self.getV_ntsequence() + "\n"

        strScenarioFasta = strScenarioFasta + "> d_gene\n"
        strScenarioFasta = strScenarioFasta + self.getD_ntsequence() + "\n"

        strScenarioFasta = strScenarioFasta + "> j_choice\n"
        strScenarioFasta = strScenarioFasta + self.getJ_ntsequence() + "\n"

        return strScenarioFasta

    #### V region methods
    def getV_fasta(self):
        strV_fasta = ""
        strV_fasta = strV_fasta + ">" + str(self.id_v_choice) + ": " + self.getV_gene_name() + ", dels 3' = " + str(
            self.getV_3_dels()) + "\n"
        strV_fasta = strV_fasta + self.getV_Region() + "\n"
        return strV_fasta

    def getV_gene_name(self):
        strEv = 'v_choice'
        name_V = self.mdlParms.Event_dict[strEv].loc[self.id_v_choice]['name']
        return name_V

    def getV_ntsequence(self):
        strEv = 'v_choice'
        seq_V = self.mdlParms.Event_dict[strEv].loc[self.id_v_choice]['value']
        return seq_V

    def getV_3_dels(self):
        strEv = 'v_3_del'
        n_v_3_del = self.mdlParms.Event_dict[strEv].loc[self.id_v_3_del]['value']
        return n_v_3_del

    def getV_Region(self):
        # seq_id=59
        strEv = 'v_choice'
        seq_V = self.mdlParms.Event_dict[strEv].loc[self.id_v_choice]['value']
        n_v_3_del = self.getV_3_dels()
        if n_v_3_del == 0:
            return (seq_V)
        # FIXME: ADD palindromic insertions
        elif n_v_3_del < 0:
            return (seq_V + 'X' * n_v_3_del)
        else:
            return (seq_V[:-n_v_3_del])

    #### J region methods
    def getJ_fasta(self):
        strJ_fasta = ""
        strJ_fasta = strJ_fasta + ">" + str(self.id_j_choice) + ": " + self.getJ_gene_name() + ", dels 5' = " + str(
            self.getJ_5_dels()) + "\n"
        strJ_fasta = strJ_fasta + self.getJ_Region() + "\n"
        return strJ_fasta

    def getJ_gene_name(self):
        strEv = 'j_choice'
        name_J = self.mdlParms.Event_dict[strEv].loc[self.id_j_choice]['name']
        return name_J

    def getJ_ntsequence(self):
        strEv = 'j_choice'
        seq_J = self.mdlParms.Event_dict[strEv].loc[self.id_j_choice]['value']
        return seq_J

    def getJ_5_dels(self):
        strEv = 'j_5_del'
        n_j_5_del = self.mdlParms.Event_dict[strEv].loc[self.id_j_5_del]['value']
        return n_j_5_del

    def getJ_Region(self):
        # seq_id=59
        strEv = 'j_choice'
        seq_J = self.mdlParms.Event_dict[strEv].loc[self.id_j_choice]['value']  # .values
        n_j_5_del = self.getJ_5_dels()
        if n_j_5_del == 0:
            return (seq_J)
        # FIXME: ADD palindromic insertions
        elif n_j_5_del < 0:
            return (seq_J + 'X' * n_j_5_del)
        else:
            return (seq_J[n_j_5_del:])

    #### D region methods
    def getD_fasta(self):
        strD_fasta = ""
        strD_fasta = strD_fasta + ">" + str(self.id_d_gene) + ": " + self.getD_gene_name() \
                     + ", " + str(self.id_d_5_del) + " dels 5' = " + str(self.getD_5_dels()) \
                     + ", " + str(self.id_d_3_del) + " dels 3' = " + str(self.getD_3_dels()) + "\n"
        strD_fasta = strD_fasta + self.getD_Region() + "\n"
        return strD_fasta

    def getD_gene_name(self):
        strEv = 'd_gene'
        name_D = self.mdlParms.Event_dict[strEv].loc[self.id_d_gene]['name']
        return name_D

    def getD_ntsequence(self):
        strEv = 'd_gene'
        seq_D = self.mdlParms.Event_dict[strEv].loc[self.id_d_gene]['value']
        return seq_D

    def getD_5_dels(self):
        strEv = 'd_5_del'
        n_d_5_del = self.mdlParms.Event_dict[strEv].loc[self.id_d_5_del]['value']
        return n_d_5_del

    def getD_3_dels(self):
        strEv = 'd_3_del'
        n_d_3_del = self.mdlParms.Event_dict[strEv].loc[self.id_d_3_del]['value']
        return n_d_3_del

    def getD_Region(self):
        strEv = 'd_gene'
        seq_D = self.mdlParms.Event_dict[strEv].loc[self.id_d_gene]['value']  # .values
        n_d_5_del = self.getD_5_dels()
        n_d_3_del = self.getD_3_dels()
        if n_d_3_del == 0:
            return (seq_D[n_d_5_del:])
        # FIXME: ADD palindromic insertions
        elif n_d_3_del < 0:
            return (seq_D[n_d_5_del:] + 'X' * n_d_3_del)
        else:
            return (seq_D[n_d_5_del:-n_d_3_del])

    #### VD region methods
    def getVD_fasta(self):
        strVD_fasta = ""
        strVD_fasta = strVD_fasta + ">" + str(self.id_vd_ins) + ", VD insertions = " + str(self.getVD_ins()) + "\n"
        strVD_fasta = strVD_fasta + self.getVD_Region() + "\n"
        return strVD_fasta

    def getVD_ins(self):
        strEv = 'vd_ins'
        n_vd_ins = self.mdlParms.Event_dict[strEv].loc[self.id_vd_ins]['value']
        return n_vd_ins

    def getVD_Region(self):
        strEv = 'vd_dinucl'
        seq_VD_dinucl = self.mdlParms.Event_dict[strEv].loc[self.vd_dinucl]['value'].values
        return (''.join(seq_VD_dinucl.tolist()))

        #### DJ region methods

    def getDJ_fasta(self):
        strDJ_fasta = ""
        strDJ_fasta = strDJ_fasta + ">" + str(self.id_dj_ins) + ", DJ insertions = " + str(self.getDJ_ins()) + "\n"
        strDJ_fasta = strDJ_fasta + self.getDJ_Region() + "\n"
        return strDJ_fasta

    def getDJ_ins(self):
        strEv = 'dj_ins'
        n_dj_ins = self.mdlParms.Event_dict[strEv].loc[self.id_dj_ins]['value']
        return n_dj_ins

    def getDJ_Region(self):
        strEv = 'dj_dinucl'
        seq_DJ_dinucl = self.mdlParms.Event_dict[strEv].loc[self.dj_dinucl]['value'].values
        return (''.join(seq_DJ_dinucl.tolist()))

        # FIXME: CHANGE NAME TO get_ScenarioProb.

    def get_EventProb(self):
        ###### v_choice
        strEvent = 'v_choice'
        da_event = self.mdl.xdata[strEvent]
        p_V = da_event[{strEvent: self.id_v_choice}]
        # print("p_V = ", p_V)
        # p_V = da_event.where(da_event['lbl__'+strEvent] == bs.getV_gene_name() ,drop=True)

        ###### j_choice
        strEvent = 'j_choice'
        da_event = self.mdl.xdata[strEvent]
        p_J = da_event[{strEvent: self.id_j_choice}]
        # print("p_J = ", p_J)
        # p_J = da_event.where(da_event['lbl__'+strEvent] == bs.getJ_gene_name() ,drop=True)

        ###### d_gene
        strEvent = 'd_gene'
        da_event = self.mdl.xdata[strEvent]
        p_DgJ = da_event[{'d_gene': self.id_d_gene, 'j_choice': self.id_j_choice}]
        # print("p_DgJ = ", p_DgJ)
        # p_DgJ

        ###### v_3_del
        strEvent = 'v_3_del'
        da_event = self.mdl.xdata[strEvent]
        p_V_3_del = da_event[{'v_3_del': self.id_v_3_del, 'v_choice': self.id_v_choice}]
        # print("p_V_3_del = ", p_V_3_del)
        # p_V_3_del

        ###### j_5_del
        strEvent = 'j_5_del'
        da_event = self.mdl.xdata[strEvent]
        p_J_5_del = da_event[{'j_5_del': self.id_j_5_del, 'j_choice': self.id_j_choice}]
        # print("p_J_5_del = ", p_J_5_del)
        # p_J_5_del

        ###### d_5_del
        strEvent = 'd_5_del'
        da_event = self.mdl.xdata[strEvent]
        p_D_5_del = da_event[{'d_5_del': self.id_d_5_del, 'd_gene': self.id_d_gene}]
        # print("p_D_5_del = ", p_D_5_del)
        # p_D_5_del

        ###### d_3_del
        strEvent = 'd_3_del'
        da_event = self.mdl.xdata[strEvent]
        p_D_3_del = da_event[{'d_3_del': self.id_d_3_del, 'd_5_del': self.id_d_5_del, 'd_gene': self.id_d_gene}]
        # print("p_D_3_del = ", p_D_3_del)
        # p_D_3_del

        ###### vd_ins
        strEvent = 'vd_ins'
        da_event = self.mdl.xdata[strEvent]
        p_VD_ins = da_event[{'vd_ins': self.id_vd_ins}]
        # print("p_VD_ins = ", p_VD_ins)

        ###### vd_dinucl
        strEvent = 'vd_dinucl'
        da_event = self.mdl.xdata[strEvent]
        # Get the last nucleotide of V region (after deletions)

        str_prev_nt = self.getV_Region()[-1]
        pd_tmp = self.mdl.parms.Event_dict[strEvent]
        prev_nt = pd_tmp.loc[pd_tmp['value'] == str_prev_nt].index.values[0]

        # for each nucleotide on inserted list
        Num_nt = 4  # 4 nucleotides A, C, G, T
        p_VD_dinucl = 1
        for curr_nt in self.vd_dinucl:
            id_dinucl = prev_nt * Num_nt + curr_nt
            prob_tmp = da_event[{'vd_dinucl': id_dinucl}]
            p_VD_dinucl = p_VD_dinucl * prob_tmp
            # print(prev_nt, curr_nt, id_dinucl, prob_tmp, p_VD_dinucl)
            prev_nt = curr_nt

        #
        # print("p_VD_dinucl = ", p_VD_dinucl)

        ###### dj_ins
        strEvent = 'dj_ins'
        da_event = self.mdl.xdata[strEvent]
        p_DJ_ins = da_event[{'dj_ins': self.id_dj_ins}]
        # print("p_DJ_ins = ", p_DJ_ins)

        ###### dj_dinucl
        strEvent = 'dj_dinucl'
        da_event = self.mdl.xdata[strEvent]
        # Get the last nucleotide of V region (after deletions)

        #        str_prev_nt = (self.getV_Region() + self.getVD_Region() + self.getD_Region() )[-1]
        str_prev_nt = (self.getJ_Region())[0]
        pd_tmp = self.mdl.parms.Event_dict[strEvent]
        prev_nt = pd_tmp.loc[pd_tmp['value'] == str_prev_nt].index.values[0]
        # print("prev_nt : ", prev_nt)

        # for each nucleotide on inserted list
        Num_nt = 4  # 4 nucleotides A, C, G, T
        p_DJ_dinucl = 1
        # self.dj_dinucl = self.dj_dinucl[::-1]
        for curr_nt in self.dj_dinucl:
            id_dinucl = prev_nt * Num_nt + curr_nt
            prob_tmp = da_event[{'dj_dinucl': id_dinucl}]
            p_DJ_dinucl = p_DJ_dinucl * prob_tmp
            # print(prev_nt, curr_nt, id_dinucl, prob_tmp, p_DJ_dinucl)
            prev_nt = curr_nt

        # print("p_DJ_dinucl = ", p_DJ_dinucl)

        p_vecE = p_V * p_J * p_DgJ * p_V_3_del * p_J_5_del * p_D_5_del * p_D_3_del * \
                 p_VD_ins * p_VD_dinucl * p_DJ_ins * p_DJ_dinucl

        return p_vecE.values

    def get_DictNicknameProbs(self):
        dictNicknameProbs = dict()
        {
            "v_choice": self.id_v_choice, \
            "j_choice": self.id_j_choice, \
            "d_gene": self.id_d_gene, \
            "v_3_del": self.id_v_3_del, \
            "d_5_del": self.id_d_5_del, \
            "d_3_del": self.id_d_3_del, \
            "j_5_del": self.id_j_5_del, \
            "vd_ins": self.id_vd_ins, \
            "vd_dinucl": self.vd_dinucl, \
            "dj_ins": self.id_dj_ins, \
            "dj_dinucl": self.dj_dinucl, \
            "mismatches": self.mismatches, \
            "mismatcheslen": self.mismatcheslen
        }
        ###### v_choice
        strEvent = 'v_choice'
        da_event = self.mdl.xdata[strEvent]
        p_V = da_event[{strEvent: self.id_v_choice}]
        dictNicknameProbs[strEvent] = p_V
        # print("p_V = ", p_V)
        # p_V = da_event.where(da_event['lbl__'+strEvent] == bs.getV_gene_name() ,drop=True)

        ###### j_choice
        strEvent = 'j_choice'
        da_event = self.mdl.xdata[strEvent]
        p_J = da_event[{strEvent: self.id_j_choice}]
        dictNicknameProbs[strEvent] = p_J
        # print("p_J = ", p_J)
        # p_J = da_event.where(da_event['lbl__'+strEvent] == bs.getJ_gene_name() ,drop=True)

        ###### d_gene
        strEvent = 'd_gene'
        da_event = self.mdl.xdata[strEvent]
        p_DgJ = da_event[{'d_gene': self.id_d_gene, 'j_choice': self.id_j_choice}]
        dictNicknameProbs[strEvent] = p_DgJ
        # print("p_DgJ = ", p_DgJ)
        # p_DgJ

        ###### v_3_del
        strEvent = 'v_3_del'
        da_event = self.mdl.xdata[strEvent]
        p_V_3_del = da_event[{'v_3_del': self.id_v_3_del, 'v_choice': self.id_v_choice}]
        dictNicknameProbs[strEvent] = p_V_3_del
        # print("p_V_3_del = ", p_V_3_del)
        # p_V_3_del

        ###### j_5_del
        strEvent = 'j_5_del'
        da_event = self.mdl.xdata[strEvent]
        p_J_5_del = da_event[{'j_5_del': self.id_j_5_del, 'j_choice': self.id_j_choice}]
        dictNicknameProbs[strEvent] = p_J_5_del

        ###### d_5_del
        strEvent = 'd_5_del'
        da_event = self.mdl.xdata[strEvent]
        p_D_5_del = da_event[{'d_5_del': self.id_d_5_del, 'd_gene': self.id_d_gene}]
        dictNicknameProbs[strEvent] = p_D_5_del

        ###### d_3_del
        strEvent = 'd_3_del'
        da_event = self.mdl.xdata[strEvent]
        p_D_3_del = da_event[{'d_3_del': self.id_d_3_del, 'd_5_del': self.id_d_5_del, 'd_gene': self.id_d_gene}]
        dictNicknameProbs[strEvent] = p_D_3_del

        ###### vd_ins
        strEvent = 'vd_ins'
        da_event = self.mdl.xdata[strEvent]
        p_VD_ins = da_event[{'vd_ins': self.id_vd_ins}]
        dictNicknameProbs[strEvent] = p_VD_ins

        ###### vd_dinucl
        strEvent = 'vd_dinucl'
        da_event = self.mdl.xdata[strEvent]
        # Get the last nucleotide of V region (after deletions)

        str_prev_nt = self.getV_Region()[-1]
        pd_tmp = self.mdl.parms.Event_dict[strEvent]
        prev_nt = pd_tmp.loc[pd_tmp['value'] == str_prev_nt].index.values[0]

        # for each nucleotide on inserted list
        Num_nt = 4  # 4 nucleotides A, C, G, T
        p_VD_dinucl = 1
        for curr_nt in self.vd_dinucl:
            id_dinucl = prev_nt * Num_nt + curr_nt
            prob_tmp = da_event[{'vd_dinucl': id_dinucl}]
            p_VD_dinucl = p_VD_dinucl * prob_tmp
            # print(prev_nt, curr_nt, id_dinucl, prob_tmp, p_VD_dinucl)
            prev_nt = curr_nt

        dictNicknameProbs[strEvent] = p_VD_dinucl

        ###### dj_ins
        strEvent = 'dj_ins'
        da_event = self.mdl.xdata[strEvent]
        p_DJ_ins = da_event[{'dj_ins': self.id_dj_ins}]
        dictNicknameProbs[strEvent] = p_DJ_ins

        ###### dj_dinucl
        strEvent = 'dj_dinucl'
        da_event = self.mdl.xdata[strEvent]
        # Get the last nucleotide of V region (after deletions)

        str_prev_nt = (self.getV_Region() + self.getVD_Region() + self.getD_Region())[-1]
        pd_tmp = self.mdl.parms.Event_dict[strEvent]
        prev_nt = pd_tmp.loc[pd_tmp['value'] == str_prev_nt].index.values[0]
        # print("prev_nt : ", prev_nt)

        # for each nucleotide on inserted list
        Num_nt = 4  # 4 nucleotides A, C, G, T
        p_DJ_dinucl = 1
        for curr_nt in self.dj_dinucl:
            id_dinucl = prev_nt * Num_nt + curr_nt
            prob_tmp = da_event[{'dj_dinucl': id_dinucl}]
            p_DJ_dinucl = p_DJ_dinucl * prob_tmp
            # print(prev_nt, curr_nt, id_dinucl, prob_tmp, p_DJ_dinucl)
            prev_nt = curr_nt

        dictNicknameProbs[strEvent] = p_DJ_dinucl

        return dictNicknameProbs

    def get_ErrorProb(self):
        r = float(self.mdl.parms.ErrorRate['SingleErrorRate'])
        L = len(self.strSeq_index)
        print("error rate: ", r, "n mismatches", self.mismatcheslen)
        # return r**(self.mismatcheslen)
        return (r / 3) ** (self.mismatcheslen)  # * (1-r)**( L - self.mismatcheslen)

