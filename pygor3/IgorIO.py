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

from .utils import *
from .IgorSQL import *
import collections

### GENERIC FUNCTIONS
# Generation of label for a simple identification of genomic template sequence.
def genLabel(strName):
    aaa = strName.split("|")
    if len(aaa) > 1 :
        return aaa[1]
    else:
        return strName
v_genLabel = np.vectorize(genLabel)

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
    # print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout = []
    while True:
        line = p.stdout.readline()
        line = line.decode("utf-8")
        stdout.append(line)
        #print (line, end='')
        if line == '' and p.poll() != None:
            break
    return ''.join(stdout)


def execute_command_generator(cmd):
    popen = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

def run_command_print(cmd):
    std_output_str = ""
    for path in execute_command_generator(cmd):
        print(path, end="")
        std_output_str = std_output_str + '\n'

    return std_output_str

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
    def __init__(self, igor_exec_path=None, igor_datadir=None,
                 igor_models_root_path=None, igor_species=None, igor_chain=None,
                 igor_model_dir_path=None,
                 igor_path_ref_genome=None, fln_genomicVs=None, fln_genomicDs=None, fln_genomicJs=None, fln_V_gene_CDR3_anchors=None, fln_J_gene_CDR3_anchors=None,
                 igor_wd=None, igor_batchname=None,
                 igor_model_parms_file=None, igor_model_marginals_file=None,
                 igor_read_seqs=None,
                 igor_threads=None,
                 igor_fln_indexed_sequences=None,
                 igor_fln_indexed_CDR3=None,
                 igor_fln_align_V_alignments=None,
                 igor_fln_align_D_alignments=None,
                 igor_fln_align_J_alignments=None,
                 igor_fln_infer_final_marginals=None,
                 igor_fln_infer_final_parms=None,
                 igor_fln_evaluate_final_marginals=None,
                 igor_fln_evaluate_final_parms=None,
                 igor_fln_output_pgen=None,
                 igor_fln_output_scenarios=None,
                 igor_fln_output_coverage=None,
                 igor_fln_generated_realizations_werr=None,
                 igor_fln_generated_seqs_werr=None,
                 igor_fln_generation_info=None,
                 igor_fln_db=None
                 ):
        # To execute IGoR externally
        self.igor_exec_path = igor_exec_path
        self.igor_datadir = igor_datadir

        # To load default models and genomic templates
        self.igor_models_root_path = igor_models_root_path # igor models paths where all species and chains are stored.
        self.igor_species = igor_species
        self.igor_chain = igor_chain

        self.igor_model_dir_path = igor_model_dir_path

        # genome references alignments
        self.igor_path_ref_genome = igor_path_ref_genome
        self.fln_genomicVs = fln_genomicVs
        self.fln_genomicDs = fln_genomicDs
        self.fln_genomicJs = fln_genomicJs
        self.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
        self.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

        self.igor_wd = igor_wd
        self.igor_batchname = igor_batchname

        self.igor_model_parms_file = igor_model_parms_file
        self.igor_model_marginals_file = igor_model_marginals_file

        self.igor_read_seqs = igor_read_seqs
        self.igor_threads = igor_threads

        # read
        self.igor_fln_indexed_sequences = igor_fln_indexed_sequences
        # aligns
        self.igor_fln_indexed_CDR3 = igor_fln_indexed_CDR3
        self.igor_fln_align_V_alignments = igor_fln_align_V_alignments
        self.igor_fln_align_J_alignments = igor_fln_align_J_alignments
        self.igor_fln_align_D_alignments = igor_fln_align_D_alignments
        # inference
        self.igor_fln_infer_final_marginals = igor_fln_infer_final_marginals
        self.igor_fln_infer_final_parms = igor_fln_infer_final_parms
        # evaluate
        self.igor_fln_evaluate_final_marginals = igor_fln_evaluate_final_marginals
        self.igor_fln_evaluate_final_parms = igor_fln_evaluate_final_parms
        # output
        self.igor_fln_output_pgen = igor_fln_output_pgen
        self.igor_fln_output_scenarios = igor_fln_output_scenarios
        self.igor_fln_output_coverage = igor_fln_output_coverage

        # TODO: NO DATABASE FIELDS FOR THIS GUYS
        self.igor_fln_generated_realizations_werr = igor_fln_generated_realizations_werr
        self.igor_fln_generated_seqs_werr = igor_fln_generated_seqs_werr
        self.igor_fln_generation_info = igor_fln_generation_info

        self.igor_fln_db = igor_fln_db

        # TODO: experimental dictionary to check status of igor batch associated files
        # almost each of these files correspond to a sql table
        self.batch_data = igor_batch_dict


        self.igor_db = IgorSqliteDB()
        self.igor_db_bs = None

        self.b_read_seqs = False
        self.b_align = False
        self.b_infer = False
        self.b_evaluate = False
        self.b_generate = False

        self.mdl = IgorModel()
        self.genomes = IgorRefGenome() #{ 'V' : IgorRefGenome(), 'D' : IgorRefGenome(), 'J' : IgorRefGenome() }

        self.igor_align_dict_options = igor_align_dict_options

        self.igor_evaluate_dict_options = igor_evaluate_dict_options

        self.igor_output_dict_options = igor_output_dict_options


        try:
            if self.igor_batchname is None:
                self.gen_random_batchname()
        except Exception as e:
            print(e)
            raise e

        try:
            if self.igor_wd is None:
                self.gen_igor_wd()
        except Exception as e:
            print(e)
            raise e

        try:
            if self.igor_exec_path is None:
                p = subprocess.Popen("which igor", shell=True, stdout=subprocess.PIPE)
                line = p.stdout.readline()
                self.igor_exec_path = line.decode("utf-8").replace('\n', '')
        except Exception as e:
            print(e)
            raise e

        try:
            if self.igor_datadir is None:
                self.run_datadir()
        except Exception as e:
            print(e)
            raise e

    def __repr__(self):
        str_repr = ""
        for key, value in self.to_dict().items():
            str_repr = str_repr + str(key) +" = "+ str(value) + "\n"
        return str_repr

    def to_dict(self):
        dicto = {
            "igor_species": self.igor_species,
            "igor_chain": self.igor_chain,
            "igor_model_dir_path": self.igor_model_dir_path,
            "igor_path_ref_genome": self.igor_path_ref_genome,
            "fln_genomicVs": self.fln_genomicVs,
            "fln_genomicDs": self.fln_genomicDs,
            "fln_genomicJs": self.fln_genomicJs,
            "fln_V_gene_CDR3_anchors": self.fln_V_gene_CDR3_anchors,
            "fln_J_gene_CDR3_anchors": self.fln_J_gene_CDR3_anchors,
            "igor_wd": self.igor_wd,
            "igor_batchname": self.igor_batchname,
            "igor_model_parms_file": self.igor_model_parms_file,
            "igor_model_marginals_file": self.igor_model_marginals_file,
            "igor_read_seqs": self.igor_read_seqs,
            "igor_threads": self.igor_threads,
            "igor_fln_indexed_sequences": self.igor_fln_indexed_sequences,
            "igor_fln_indexed_CDR3": self.igor_fln_indexed_CDR3,
            "igor_fln_align_V_alignments": self.igor_fln_align_V_alignments,
            "igor_fln_align_J_alignments": self.igor_fln_align_J_alignments,
            "igor_fln_align_D_alignments": self.igor_fln_align_D_alignments,
            "igor_fln_infer_final_marginals": self.igor_fln_infer_final_marginals,
            "igor_fln_infer_final_parms": self.igor_fln_infer_final_parms,
            "igor_fln_evaluate_final_marginals": self.igor_fln_evaluate_final_marginals,
            "igor_fln_evaluate_final_parms": self.igor_fln_evaluate_final_parms,
            "igor_fln_output_pgen": self.igor_fln_output_pgen,
            "igor_fln_output_scenarios": self.igor_fln_output_scenarios,
            "igor_fln_output_coverage": self.igor_fln_output_coverage,

            "igor_fln_generated_realizations_werr": self.igor_fln_generated_realizations_werr,
            "igor_fln_generated_seqs_werr": self.igor_fln_generated_seqs_werr,
            "igor_fln_generation_info": self.igor_fln_generation_info,

            "igor_fln_db": self.igor_fln_db,
            "b_read_seqs": self.b_read_seqs,
            "b_align" : self.b_align,
            "b_infer": self.b_infer,
            "b_evaluate": self.b_evaluate,
            "b_generate": self.b_generate
        }
        return dicto

    def load_IgorRefGenome(self, igor_path_ref_genome=None):
        # FIXME: THERE ARE 2 OPTIONS HERE:
        # 1. From template directory self.igor_path_ref_genome
        if igor_path_ref_genome is not None:
            self.igor_path_ref_genome = igor_path_ref_genome
        self.genomes = IgorRefGenome.load_from_path(self.igor_path_ref_genome)
        # TODO: FIND A BETTER WAY TO SYNCHRONIZE NAMES (FORWARD AND BACKWARD)
        self.fln_genomicVs = self.genomes.fln_genomicVs
        self.fln_genomicDs = self.genomes.fln_genomicDs
        self.fln_genomicJs = self.genomes.fln_genomicJs

        self.fln_V_gene_CDR3_anchors = self.genomes.fln_V_gene_CDR3_anchors
        self.fln_J_gene_CDR3_anchors = self.genomes.fln_J_gene_CDR3_anchors

        # # 2. Or directly from files.
        # self.genomes = IgorRefGenome()
        # self.genomes.fln_genomicVs = self.fln_genomicVs
        # self.genomes.fln_genomicDs = self.fln_genomicDs
        # self.genomes.fln_genomicJs = self.fln_genomicJs

    def make_model_default_VJ_from_genomes(self, igor_path_ref_genome=None):
        try:
            self.load_IgorRefGenome(igor_path_ref_genome=igor_path_ref_genome)
            mdl_parms = IgorModel_Parms.make_default_VJ(self.genomes.df_genomicVs, self.genomes.df_genomicJs)
            mdl_marginals = IgorModel_Marginals.make_uniform_from_parms(mdl_parms)
            self.mdl = IgorModel.load_from_parms_marginals_object(mdl_parms, mdl_marginals)
        except Exception as e:
            print("ERROR: ", e)

    def make_model_default_VDJ_from_genomes(self, igor_path_ref_genome=None):
        try:
            self.load_IgorRefGenome(igor_path_ref_genome=igor_path_ref_genome)
            mdl_parms = IgorModel_Parms.make_default_VDJ(self.genomes.df_genomicVs, self.genomes.df_genomicDs, self.genomes.df_genomicJs)
            mdl_marginals = IgorModel_Marginals.make_uniform_from_parms(mdl_parms)
            self.mdl = IgorModel.load_from_parms_marginals_object(mdl_parms, mdl_marginals)
        except Exception as e:
            print("ERROR: ", e)

    def load_IgorModel(self):
        if ( (self.igor_species is None ) or (self.igor_chain is None)):
            self.mdl = IgorModel(model_parms_file = self.igor_model_parms_file, model_marginals_file=self.igor_model_marginals_file)
        else :
            self.mdl = IgorModel.load_default(self.igor_species, igor_option_path_dict[self.igor_chain])

    def load_IgorModel_from_infer_files(self):
        try:
            self.mdl = IgorModel(model_parms_file = self.igor_fln_infer_final_parms, model_marginals_file=self.igor_fln_infer_final_marginals)
        except Exception as e:
            print("ERROR: IgorTask.load_IgorModel_inferred:")
            print(e)

    @classmethod
    def default_model(cls, specie, chain, model_parms_file=None, model_marginals_file=None):
        """Return an IgorTask object"""
        cls = IgorTask()
        cls.igor_species = specie
        cls.igor_chain = chain
        #cls.igor_modeldirpath =  model_parms_file
        cls.run_datadir()
        cls.igor_model_dir_path = cls.igor_models_root_path +"/" + cls.igor_species + "/" + igor_option_path_dict[cls.igor_chain]

        if model_parms_file is None:
            cls.igor_model_parms_file = cls.igor_model_dir_path+ "/models/model_parms.txt"
            cls.igor_model_marginals_file = cls.igor_model_dir_path + "/models/model_marginals.txt"
            cls.mdl = IgorModel(model_parms_file=cls.igor_model_parms_file, model_marginals_file=cls.igor_model_marginals_file)
        return cls

    def gen_igor_wd(self):
        p = subprocess.Popen("pwd", shell=True, stdout=subprocess.PIPE)
        line = p.stdout.readline()
        self.igor_wd = line.decode("utf-8").replace('\n', '')

    def gen_random_batchname(self):
        p = subprocess.Popen("head /dev/urandom | tr -dc A-Za-z0-9 | head -c10", shell=True, stdout=subprocess.PIPE)
        line = p.stdout.readline()
        self.igor_batchname = "dataIGoR" + line.decode("utf-8").replace('\n', '')

    def update_model_filenames(self, model_path=None):

        try:
            # if model_path is None use the self.igor_model_dir_path
            if model_path is None:
                # use previously defined igor_model_dir_path
                if self.igor_model_dir_path is None:
                    # if wasn't defined use the current directory
                    model_path = "."
                    if (not (self.igor_species is None) ) and (not (self.igor_chain is None)):
                        self.run_datadir()
                        self.igor_model_dir_path = self.igor_models_root_path + "/" + self.igor_species + "/" + \
                                                       igor_option_path_dict[self.igor_chain]
            else:
                # if a model_path is provided then override it
                self.igor_model_dir_path = model_path

            self.igor_model_parms_file = self.igor_model_dir_path + "/models/model_parms.txt"
            self.igor_model_marginals_file = self.igor_model_dir_path + "/models/model_marginals.txt"
            self.igor_path_ref_genome = self.igor_model_dir_path + "/ref_genome/"
        except Exception as e:
            print("WARNING: IgorTask.update_model_filenames", e, self.igor_model_dir_path)


    def update_ref_genome(self, igor_path_ref_genome=None):
        if igor_path_ref_genome is not None:
            self.igor_path_ref_genome = igor_path_ref_genome
        self.genomes.update_fln_names(path_ref_genome=self.igor_path_ref_genome) #fln_genomicVs = self.igor_path_ref_genome + ""
        self.fln_genomicVs = self.genomes.fln_genomicVs
        self.fln_genomicJs = self.genomes.fln_genomicJs
        self.fln_genomicDs = self.genomes.fln_genomicDs
        self.fln_V_gene_CDR3_anchors = self.genomes.fln_V_gene_CDR3_anchors
        self.fln_J_gene_CDR3_anchors = self.genomes.fln_J_gene_CDR3_anchors


    def update_batch_filenames(self):
        # reads
        if self.igor_wd is None:
            self.gen_igor_wd()
        if self.igor_batchname is None:
            self.gen_random_batchname()

        self.igor_fln_indexed_sequences = self.igor_wd + "/aligns/" + self.igor_batchname + "_indexed_sequences.csv"

        # aligns
        self.igor_fln_indexed_CDR3 = self.igor_wd + "/aligns/" + self.igor_batchname + "_indexed_CDR3s.csv"

        self.igor_fln_align_V_alignments = self.igor_wd + "/aligns/" + self.igor_batchname + "_V_alignments.csv"
        self.igor_fln_align_J_alignments = self.igor_wd + "/aligns/" + self.igor_batchname + "_J_alignments.csv"
        self.igor_fln_align_D_alignments = self.igor_wd + "/aligns/" + self.igor_batchname + "_D_alignments.csv"

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
        return run_command(cmd)

    def run_datadir(self):
        cmd = self.igor_exec_path+ " -getdatadir"
        self.igor_datadir = run_command(cmd).replace('\n','')
        self.igor_models_root_path = self.igor_datadir + "/models/"

    def run_read_seqs(self, igor_read_seqs=None):
        if igor_read_seqs is not None:
            self.igor_read_seqs = igor_read_seqs

        "igor -set_wd $WDPATH -batch foo -read_seqs ../demo/murugan_naive1_noncoding_demo_seqs.txt"
        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        cmd = cmd + " -read_seqs " + self.igor_read_seqs
        # TODO: if self.igor_read_seqs extension fastq then convert to csv and copy and create the file in aligns. Overwrite if necesserasy
        print(cmd)
        cmd_stdout = run_command(cmd)
        self.b_read_seqs = True # FIXME: If run_command success then True
        return cmd_stdout

    def run_align(self, igor_read_seqs=None):
        #"igor -set_wd ${tmp_dir} -batch ${randomBatch} -species
        # ${species} -chain ${chain} -align --all"
        import os.path

        if igor_read_seqs is not None:
            self.igor_read_seqs = igor_read_seqs

        if self.b_read_seqs is False:
            self.run_read_seqs(igor_read_seqs=igor_read_seqs)
        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        # TODO: USE COSTUM MODEL OR USE SPECIFIED SPECIES?
        # I think that the safests is to use the
        # FIXME: CHANGE TO CUSTOM GENOMICS
        cmd = cmd + " -set_genomic "

        if os.path.isfile(self.genomes.fln_genomicVs):
            cmd = cmd + " --V " + self.genomes.fln_genomicVs
        if os.path.isfile(self.genomes.fln_genomicDs):
            cmd = cmd + " --D " + self.genomes.fln_genomicDs
        if os.path.isfile(self.genomes.fln_genomicJs):
            cmd = cmd + " --J " + self.genomes.fln_genomicJs

        cmd = cmd + " -set_CDR3_anchors "

        if os.path.isfile(self.genomes.fln_V_gene_CDR3_anchors):
            cmd = cmd + " --V " + self.genomes.fln_V_gene_CDR3_anchors
        if os.path.isfile(self.genomes.fln_J_gene_CDR3_anchors):
            cmd = cmd + " --J " + self.genomes.fln_J_gene_CDR3_anchors

        cmd = cmd + " -align " + command_from_dict_options(self.igor_align_dict_options)
        #return cmd
        print(cmd)
        cmd_stdout = run_command_print(cmd)
        #run_command_no_output(cmd)
        self.b_align = True # FIXME: If run_command success then True
        return cmd_stdout

    def run_evaluate(self, igor_read_seqs=None, N_scenarios=None):
        # "igor -set_wd $WDPATH -batch foo -species human -chain beta
        # -evaluate -output --scenarios 10"
        print(self.to_dict())
        import os.path
        if igor_read_seqs is not None:
            self.igor_read_seqs = igor_read_seqs

        if self.b_align is False:
            self.run_align(igor_read_seqs=self.igor_read_seqs)

        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        # TODO: USE COSTUM MODEL OR USE SPECIFIED SPECIES?
        # I think that the safests is to use the
        # cmd = cmd + " -species " + self.igor_species
        # cmd = cmd + " -chain " + self.igor_chain
        cmd = cmd + " -set_custom_model " + self.igor_model_parms_file + " " + self.igor_model_marginals_file

        # here the evaluation
        self.igor_output_dict_options["--scenarios"]['active'] = True
        if N_scenarios is not None:
            self.igor_output_dict_options["--scenarios"]['value'] = str(N_scenarios)
        self.igor_output_dict_options["--Pgen"]['active'] = True
        cmd = cmd + " -evaluate " + command_from_dict_options(self.igor_evaluate_dict_options)
        cmd = cmd + " -output " + command_from_dict_options(self.igor_output_dict_options)
        # return cmd
        print(cmd)
        # FIXME: REALLY BIG FLAW USE DICTIONARY FOR THE SPECIE AND CHAIN
        # self.mdl = IgorModel.load_default(self.igor_species, igor_option_path_dict[self.igor_chain], modelpath=self.igor_models_root_path)
        run_command(cmd)
        # run_command_no_output(cmd)
        # self.b_evaluate = True # FIXME: If run_command success then Truerun_infer

    def run_pgen(self, igor_read_seqs=None):
        # "igor -set_wd $WDPATH -batch foo -species human -chain beta
        # -evaluate -output --scenarios 10"
        print(self.to_dict())
        import os.path
        if igor_read_seqs is not None:
            self.igor_read_seqs = igor_read_seqs

        if self.b_align is False:
            self.run_align(igor_read_seqs=self.igor_read_seqs)

        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        # TODO: USE COSTUM MODEL OR USE SPECIFIED SPECIES?
        # I think that the safests is to use the
        cmd = cmd + " -set_custom_model " + self.igor_model_parms_file + " " + self.igor_model_marginals_file

        # here the evaluation
        self.igor_output_dict_options["--scenarios"]['active'] = False
        self.igor_output_dict_options["--Pgen"]['active'] = True
        cmd = cmd + " -evaluate " + command_from_dict_options(self.igor_evaluate_dict_options)
        cmd = cmd + " -output " + command_from_dict_options(self.igor_output_dict_options)
        # return cmd
        print(cmd)
        # FIXME: REALLY BIG FLAW USE DICTIONARY FOR THE SPECIE AND CHAIN
        # self.mdl = IgorModel.load_default(self.igor_species, igor_option_path_dict[self.igor_chain], modelpath=self.igor_models_root_path)
        run_command(cmd)
        # run_command_no_output(cmd)
        # self.b_evaluate = True # FIXME: If run_command success then Truerun_infer

    def run_scenarios(self, igor_read_seqs=None, N_scenarios=None):
        #"igor -set_wd $WDPATH -batch foo -species human -chain beta
        # -evaluate -output --scenarios 10"
        print(self.to_dict())
        import os.path
        if igor_read_seqs is not None:
            self.igor_read_seqs = igor_read_seqs

        if self.b_align is False:
            self.run_align(igor_read_seqs=self.igor_read_seqs)

        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        # TODO: USE COSTUM MODEL OR USE SPECIFIED SPECIES?
        # I think that the safests is to use the
        cmd = cmd + " -set_custom_model " + self.igor_model_parms_file + " " + self.igor_model_marginals_file

        # here the evaluation
        self.igor_output_dict_options["--scenarios"]['active'] = True
        if N_scenarios is not None:
            self.igor_output_dict_options["--scenarios"]['value'] = str(N_scenarios)
        self.igor_output_dict_options["--Pgen"]['active'] = False
        cmd = cmd + " -evaluate " + command_from_dict_options(self.igor_evaluate_dict_options)
        cmd = cmd + " -output " + command_from_dict_options(self.igor_output_dict_options)
        #return cmd
        print(cmd)
        # FIXME: REALLY BIG FLAW USE DICTIONARY FOR THE SPECIE AND CHAIN
        # self.mdl = IgorModel.load_default(self.igor_species, igor_option_path_dict[self.igor_chain], modelpath=self.igor_models_root_path)
        # run_command(cmd)
        run_command_print(cmd)

    def run_infer(self, igor_read_seqs=None):
        #"igor -set_wd $WDPATH -batch foo -species human -chain beta
        # -evaluate -output --scenarios 10"

        if igor_read_seqs is not None:
            self.igor_read_seqs = igor_read_seqs

        if self.b_align is False:
            self.run_align(igor_read_seqs=igor_read_seqs)
            print("Alignment finished!")

        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        # TODO: USE COSTUM MODEL OR USE SPECIFIED SPECIES?
        # I think that the safests is to use the
        # cmd = cmd + " -species " + self.igor_species
        # cmd = cmd + " -chain " + self.igor_chain
        cmd = cmd + " -set_custom_model " + self.igor_model_parms_file + " " + self.igor_model_marginals_file
        # here the evaluation
        cmd = cmd + " -infer " #+ command_from_dict_options(self.igor_output_dict_options)
        #return cmd
        print(cmd)
        # FIXME: REALLY BIG FLAW USE DICTIONARY FOR THE SPECIE AND CHAIN
        # self.mdl = IgorModel.load_default(self.igor_species, igor_option_path_dict[self.igor_chain], modelpath=self.igor_models_root_path)
        self.mdl = IgorModel(model_parms_file=self.igor_model_parms_file, model_marginals_file=self.igor_model_marginals_file)
        # output = run_command(cmd)
        output = run_command_print(cmd)
        #run_command_no_output(cmd)
        self.b_infer = True # FIXME: If run_command success then True
        return output

    def run_generate(self, N_seqs=None):
        cmd = self.igor_exec_path
        cmd = cmd + " -set_wd " + self.igor_wd
        cmd = cmd + " -batch " + self.igor_batchname
        cmd = cmd + " -set_custom_model " + self.igor_model_parms_file + " " + self.igor_model_marginals_file
        if N_seqs is not None:
            cmd = cmd + " -generate " + str(N_seqs)
        else:
            cmd = cmd + " -generate "
        print(cmd)

        # run_command(cmd)
        run_command_print(cmd)
        path_generated = self.igor_wd + "/" + self.igor_batchname + "_generated/"
        self.igor_fln_generated_realizations_werr = path_generated + "generated_realizations_werr.csv"
        self.igor_fln_generated_seqs_werr = path_generated + "generated_seqs_werr.csv"
        self.igor_fln_generation_info = path_generated + "generated_seqs_werr.out"
        self.b_generate = True

        # FIXME: LOAD TO DATABASE CREATE PROPER TABLES FOR THIS
        # import pandas as pd
        df = pd.read_csv(self.igor_fln_generated_seqs_werr, delimiter=';').set_index('seq_index')
        return df

    def run_generate_to_dataframe(self, N):
        self.run_generate(self, N)

        # FIXME: LOAD TO DATABASE CREATE PROPER TABLES FOR THIS
        # import pandas as pd
        df = pd.read_csv(self.igor_fln_generated_seqs_werr, delimiter=';').set_index('seq_index')
        return df

    def run_clean_batch(self):
        cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_evaluate"
        run_command_no_output(cmd)
        cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_output"
        run_command_no_output(cmd)
        cmd = "rm " + self.igor_wd + "/aligns/" + self.igor_batchname + "*.csv"
        run_command_no_output(cmd)
        cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_generated"
        run_command_no_output(cmd)

    def create_db(self, igor_fln_db=None):
        if igor_fln_db is not None:
            self.igor_fln_db = igor_fln_db
        self.igor_db = IgorSqliteDB.create_db(self.igor_fln_db)

    def load_db_from_indexed_sequences(self):
        self.igor_db.load_IgorIndexedSeq_FromCSV(self.igor_fln_indexed_sequences)

    # load genome templates from fasta and csv files.
    def load_db_from_genomes(self):
        print("Loading Gene templates ...")
        self.igor_db.load_IgorGeneTemplate_FromFASTA("V", self.genomes.fln_genomicVs)
        self.igor_db.load_IgorGeneTemplate_FromFASTA("J", self.genomes.fln_genomicJs)
        try:
            self.igor_db.load_IgorGeneTemplate_FromFASTA("D", self.genomes.fln_genomicDs)
        except Exception as e:
            print(e)
            print("No D gene template found in batch files structure")
            pass
        # load
        print("loading Anchors data ...")
        # try:
        self.igor_db.load_IgorGeneAnchors_FromCSV("V", self.genomes.fln_V_gene_CDR3_anchors)
        self.igor_db.load_IgorGeneAnchors_FromCSV("J", self.genomes.fln_J_gene_CDR3_anchors)
        # except Exception as e:
        #     print("ERROR : ", e)

    def load_db_from_alignments(self):
        print(self.igor_fln_align_V_alignments)
        self.igor_db.load_IgorAlignments_FromCSV("V", self.igor_fln_align_V_alignments)
        self.igor_db.load_IgorAlignments_FromCSV("J", self.igor_fln_align_J_alignments)
        try:
            self.igor_db.load_IgorAlignments_FromCSV("D", self.igor_fln_align_D_alignments)
        except Exception as e:
            print(e)
            print("Couldn't load D gene alignments!")
            pass
        print("Alignments loaded in database in "+str(self.igor_fln_db))

    def load_db_from_models(self, mdl=None):
        # self.load_IgorModel()
        try:
            if self.igor_db.Q_model_in_db():
                print("WARNING: Overwriting previous model in database ", self.igor_fln_db)
                self.igor_db.delete_IgorModel_Tables()
            if mdl is None:
                self.igor_db.load_IgorModel(self.mdl)
            else:
                self.igor_db.load_IgorModel(mdl)
        except Exception as e:
            print("Couldn't load model to database from IgorModel object")
            print("ERROR: ", e)

    def load_db_from_inferred_model(self):
        self.load_IgorModel_from_infer_files()
        try:
            self.igor_db.load_IgorModel(self.mdl)
        except Exception as e:
            print("Couldn't load model to database from IgorModel object")
            print("ERROR: ", e)

    def load_db_from_indexed_cdr3(self):
        print(self.igor_fln_indexed_CDR3)
        self.igor_db.load_IgorIndexedCDR3_FromCSV(self.igor_fln_indexed_CDR3)

    def load_db_from_bestscenarios(self):
        print(self.igor_fln_output_scenarios)
        self.igor_db.load_IgorBestScenarios_FromCSV(self.igor_fln_output_scenarios, self.mdl)

    def load_db_from_pgen(self):
        print(self.igor_fln_output_pgen)
        self.igor_db.load_IgorPgen_FromCSV(self.igor_fln_output_pgen)

    def load_mdl_from_db(self):
        try:
            self.mdl = self.igor_db.get_IgorModel()
        except Exception as e:
            print("WARNING: Igor Model was not found in ", self.igor_fln_db)
            pass
        # return self.mdl

    # def get_IgorModel_from_db(self):
    #     self.mdl = self.igor_db.get_IgorModel()
    #     return self.mdl

    # FIXME: this method should be deprecated!!!
    def load_VDJ_database(self, flnIgorSQL):
        self.flnIgorSQL = flnIgorSQL
        self.igor_db = IgorSqliteDB(flnIgorSQL)
        # FIXME :EVERYTHING
        flnIgorIndexedSeq = self.igor_wd+"/aligns/"+self.igor_batchname+"_indexed_sequences.csv"
        # FIXME PATH AND OPTIONS NEED TO BE CONSISTENT
        IgorModelPath = self.igor_models_root_path + self.igor_species + "/" \
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

    def from_db_get_naive_align_dict_by_seq_index(self, seq_index):
        indexed_sequence = self.igor_db.get_IgorIndexedSeq_By_seq_index(seq_index)
        indexed_sequence.offset = 0

        best_v_align_data = self.igor_db.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
        best_j_align_data = self.igor_db.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)

        try:
            best_d_align_data = self.igor_db.get_best_IgorAlignment_data_By_seq_index('D', indexed_sequence.seq_index)
            vdj_naive_alignment = {'V': best_v_align_data,
                                   'D': best_d_align_data,
                                   'J': best_j_align_data}
            v_align_data_list = self.igor_db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
            # print('V', len(v_align_data_list), [ii.score for ii in v_align_data_list])
            d_align_data_list = self.igor_db.get_IgorAlignment_data_list_By_seq_index('D', indexed_sequence.seq_index)
            # print('D', len(d_align_data_list), [ii.score for ii in d_align_data_list])
            j_align_data_list = self.igor_db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
            # print('J', len(j_align_data_list), [ii.score for ii in j_align_data_list])
            # 1. Choose the highest score then check if this one is the desire range.
            # if there is an overlap
            # calculate score without overlap. If overlap
            # if hightest score
            for i, d_align_data in enumerate(d_align_data_list):
                # Check if D is btwn V and J position
                if (best_v_align_data.offset_3_p <= d_align_data.offset_5_p) and (
                        d_align_data.offset_3_p <= best_j_align_data.offset_5_p):
                    # vdj_naive_alignment['D'+str(i)] = d_align_data
                    vdj_naive_alignment['D'] = d_align_data
                    break

        except Exception as e:
            print(e)
            print("No d gene alignments found!")
            vdj_naive_alignment = {'V': best_v_align_data,
                                   'J': best_j_align_data}
            v_align_data_list = self.igor_db.get_IgorAlignment_data_list_By_seq_index('V', indexed_sequence.seq_index)
            print('V', len(v_align_data_list), [ii.score for ii in v_align_data_list])
            j_align_data_list = self.igor_db.get_IgorAlignment_data_list_By_seq_index('J', indexed_sequence.seq_index)
            print('J', len(j_align_data_list), [ii.score for ii in j_align_data_list])
            pass

        return indexed_sequence, vdj_naive_alignment

    def from_db_str_fasta_naive_align_by_seq_index(self, seq_index):
        """ Given an Sequence index and the corresponding alignments vj/ vdj
            return a string with considering only offset"""

        fasta_list = list()
        indexed_sequence, vdj_alignments_dict = self.from_db_get_naive_align_dict_by_seq_index(seq_index)
        indexed_sequence.sequence = indexed_sequence.sequence.lower()
        # add mismatches in sequence.
        s = list(indexed_sequence.sequence)
        for key_align in vdj_alignments_dict.keys():
            for pos_mis in vdj_alignments_dict[key_align].mismatches:
                s[pos_mis] = s[pos_mis].upper()
        indexed_sequence.sequence = "".join(s)

        str_fasta = ""
        min_offset_key = min(vdj_alignments_dict.keys(), key=lambda x: vdj_alignments_dict[x].offset)  # .offset
        min_offset = vdj_alignments_dict[min_offset_key].offset
        min_offset = min(indexed_sequence.offset, min_offset)

        delta_offset = indexed_sequence.offset - min_offset
        str_prefix = '-' * (delta_offset)
        str_fasta_sequence = str_prefix + indexed_sequence.sequence
        # print(str_fasta_sequence)
        str_fasta = str_fasta + "> " + str(indexed_sequence.seq_index)
        str_fasta_description = str_fasta
        str_fasta = str_fasta + "\n"
        str_fasta = str_fasta + str_fasta_sequence + "\n"

        fasta_list.append([str_fasta_description, str_fasta_sequence])

        for key in vdj_alignments_dict.keys():
            vdj_alignments_dict[key].strGene_seq = vdj_alignments_dict[key].strGene_seq.lower()
            delta_offset = vdj_alignments_dict[key].offset - min_offset
            str_prefix = '-' * (delta_offset)
            str_fasta_sequence = str_prefix + vdj_alignments_dict[key].strGene_seq
            # print(str_fasta_sequence)
            str_fasta_description = "> " + key + ", " + vdj_alignments_dict[key].strGene_name
            str_fasta = str_fasta + str_fasta_description + "\n"
            str_fasta = str_fasta + str_fasta_sequence + "\n"

            fasta_list.append([str_fasta_description, str_fasta_sequence])

            offset_5_p = vdj_alignments_dict[key].offset_5_p - min_offset
            offset_3_p = vdj_alignments_dict[key].offset_3_p - min_offset
            # print("delta_offset : ", delta_offset)
            # print("offset_5_p : ", vdj_alignments_dict[key].offset_5_p, offset_5_p)
            # print("offset_3_p : ", vdj_alignments_dict[key].offset_3_p, offset_3_p)
            str_prefix_2 = '-' * (offset_5_p + 1)
            str_fasta_sequence2 = str_prefix_2 + str_fasta_sequence[offset_5_p + 1:offset_3_p + 1]
            str_fasta_description2 = "> " + vdj_alignments_dict[key].strGene_name + ", score : " + str(vdj_alignments_dict[key].score)
            str_fasta = str_fasta + str_fasta_description2 + "\n"
            str_fasta = str_fasta + str_fasta_sequence2 + "\n"

            fasta_list.append([str_fasta_description2, str_fasta_sequence2])

            # TODO ADD MISMATCHES
            # align = vdj_alignments_dict[key]
            # align mismatches are in indexed sequence reference I need to convert it to gene reference given the alignment
            # given the align.offset
            # pos_in_gene  = pos_in_seq - align.offset
            # pos_in_gene = cdr3 - align.offset

        # FIXME: make a list of tuples [(description_0, sequence_0), ..., (description_i, sequence_i), ..., (description_N, sequence_N)]

        sequence_len_list = list(map(lambda x: len(x[1]), fasta_list))
        max_seq_len = max(sequence_len_list)
        for fasta_rec in fasta_list:
            len_fasta_rec_seq = len(fasta_rec[1])
            if len_fasta_rec_seq < max_seq_len:
                #         print(fasta_rec)
                ngaps = max_seq_len - len_fasta_rec_seq
                str_ngaps = str(ngaps * '-')
                fasta_rec[1] = fasta_rec[1] + str_ngaps

        str_fasta = ""
        str_fasta = '\n'.join( [fasta_rec[0]+"\n"+fasta_rec[1] for fasta_rec in fasta_list] )
        return str_fasta #, fasta_list

    def from_db_plot_naive_align_by_seq_index(self, seq_index):
        import Bio.AlignIO
        import io
        aaa = self.from_db_str_fasta_naive_align_by_seq_index(seq_index)
        aln = Bio.AlignIO.read(io.StringIO(aaa), 'fasta')
        view_alignment(aln)

    def export_to_igorfiles(self):
        print("Export: ")
        #--- 1. Indexed Sequences
        if self.igor_db.Q_sequences_in_db() and not (self.igor_fln_indexed_sequences is None):
            try:
                self.igor_db.write_IgorIndexedSeq_to_CSV(self.igor_fln_indexed_sequences)
            except Exception as e:
                print("ERROR: write_IgorIndexedSeq_to_CSV", e)
        else:
            print("No IgorIndexedSeq Table not exported")

        #--- 2. Gene Templates
        if self.igor_db.Q_ref_genome_in_db_by_gene("V") and not (self.fln_genomicVs is None):
            try:
                self.igor_db.write_IgorGeneTemplate_to_fasta("V", self.fln_genomicVs)
            except Exception as e:
                print("ERROR: write_IgorGeneTemplate_to_fasta V", e)
        else:
            print("No IgorGeneTemplate V Table")

        if self.igor_db.Q_ref_genome_in_db_by_gene("J") and not (self.fln_genomicJs is None):
            try:
                self.igor_db.write_IgorGeneTemplate_to_fasta("J", self.fln_genomicJs)
            except Exception as e:
                print("ERROR: write_IgorGeneTemplate_to_fasta J", e)
        else:
            print("No IgorGeneTemplate J Table")

        if self.igor_db.Q_ref_genome_in_db_by_gene("D") and not (self.fln_genomicDs is None):
            try:
                self.igor_db.write_IgorGeneTemplate_to_fasta("D", self.fln_genomicDs)
            except Exception as e:
                print("ERROR: write_IgorGeneTemplate_to_fasta D", e)
        else:
            print("No IgorGeneTemplate D Table")

        if self.igor_db.Q_CDR3_Anchors_in_db("V") and not (self.fln_V_gene_CDR3_anchors is None):
            try:
                self.igor_db.write_IgorGeneAnchors_to_CSV("V", self.fln_V_gene_CDR3_anchors)
            except Exception as e:
                print("ERROR: write_IgorGeneAnchors_to_CSV V", e)
        else:
            print("No IgorGeneAnchors V Table")

        if self.igor_db.Q_CDR3_Anchors_in_db("J") and not (self.fln_J_gene_CDR3_anchors is None):
            try:
                self.igor_db.write_IgorGeneAnchors_to_CSV("J", self.fln_J_gene_CDR3_anchors)
            except Exception as e:
                print("ERROR: write_IgorGeneAnchors_to_CSV J", e)
        else:
            print("No IgorGeneAnchors J Table")



        #--- 3. Alignments
        if self.igor_db.Q_align_in_db():
            # b_igor_alignments
            if self.igor_db.Q_align_in_db_by_gene("V") and not (self.igor_fln_align_V_alignments is None):
                try:
                    self.igor_db.write_IgorAlignments_to_CSV("V", self.igor_fln_align_V_alignments)
                except Exception as e:
                    print("ERROR: write_IgorAlignments_to_CSV V", e)

            if self.igor_db.Q_align_in_db_by_gene("J") and not (self.igor_fln_align_J_alignments is None):
                try:
                    self.igor_db.write_IgorAlignments_to_CSV("J", self.igor_fln_align_J_alignments)
                except Exception as e:
                    print("ERROR: write_IgorAlignments_to_CSV J", e)

            if self.igor_db.Q_align_in_db_by_gene("D") and not (self.igor_fln_align_D_alignments is None):
                try:
                    self.igor_db.write_IgorAlignments_to_CSV("D", self.igor_fln_align_D_alignments)
                except Exception as e:
                    print("ERROR: write_IgorAlignments_to_CSV D", e)
            try:
                self.igor_db.write_IgorIndexedCDR3_to_CSV(self.igor_fln_indexed_CDR3)
            except Exception as e:
                print("WARNING: No indexed CDR3 files found", self.igor_fln_indexed_CDR3)
                print(e)
                pass

        # --- 4. Export Igor Model
        if self.igor_db.Q_model_in_db():
            if (not (self.igor_model_parms_file is None)) and (not (self.igor_model_marginals_file is None)):
                try:
                    self.igor_db.write_IgorModel_to_TXT(self.igor_model_parms_file, self.igor_model_marginals_file)
                except Exception as e:
                    print("ERROR: write_IgorModel_to_TXT ",e)
            else:
                print("ERROR: igor_model_parms_file or igor_model_marginals_file not specified.")
        else:
            print("No Models Tables")

        # --- 5. Export Igor Model
        if self.igor_db.Q_IgorPgen_in_db() and not (self.igor_fln_output_pgen is None):
            try:
                self.igor_db.write_IgorPgen_to_CSV(self.igor_fln_output_pgen)
            except Exception as e:
                print("ERROR: write_IgorPgen_to_CSV ", e)

        # --- 6. Export Igor Model
        # b_igor_scenarios
        if self.igor_db.Q_IgorBestScenarios_in_db() and not (self.igor_fln_output_scenarios is None):
            try:
                self.igor_db.write_IgorBestScenarios_to_CSV(self.igor_fln_output_scenarios)
            except Exception as e:
                print("ERROR: write_IgorBestScenarios_to_CSV ", e)

        # 1.1 if Alignments


    #### AIRR methods ###
    def parse_scenarios_to_airr(self, igor_fln_output_scenarios, airr_fln_output_scenarios):
        # 1. Read header of and make a list
        open(igor_fln_output_scenarios)
        # 2.
        pass


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
            cls.seq_index = int(sqlRecord[0])
            cls.sequence = str(sqlRecord[1]).replace('\n', '')
        except Exception as e:
            print(e)
            raise e
        return cls

### IGOR ALIGNMENTS  ####
class IgorRefGenome:
    def __init__(self):
        # FIXME: find a better way to add a default value for this and also the "/" separator
        self.path_ref_genome = "."

        self.fln_genomicVs = None # "genomicVs.fasta"
        self.fln_genomicDs = None # "genomicDs.fasta"
        self.fln_genomicJs = None # "genomicJs.fasta"

        self.fln_V_gene_CDR3_anchors = None # "V_gene_CDR3_anchors.csv"
        self.fln_J_gene_CDR3_anchors = None # "J_gene_CDR3_anchors.csv"

        self.df_genomicVs = None
        self.df_genomicDs = None
        self.df_genomicJs = None


        self.dict_genomicVs = None #(self.df_genomicVs.set_index('name').to_dict())['value']
        self.dict_genomicDs = None
        self.dict_genomicJs = None

        self.df_V_ref_genome = None
        self.df_J_ref_genome = None

    @classmethod
    def load_FromSQLRecord_list(cls, sqlrecords_genomicVs = None, sqlrecords_genomicDs = None, sqlrecords_genomicJs = None,
        sqlrecords_V_gene_CDR3_anchors = None, sqlrecords_J_gene_CDR3_anchors = None):
        cls = IgorRefGenome()
        # TODO: make query to database

        cls.df_genomicVs = pd.DataFrame.from_records(sqlrecords_genomicVs, columns=['id', 'name', 'value']).set_index('id')

        # Fasta to dataframe
        try:
            # df_V_anchors = pd.read_csv(self.fln_V_gene_CDR3_anchors, sep=';')
            df_V_anchors = pd.DataFrame.from_records(sqlrecords_V_gene_CDR3_anchors, columns=['id', 'gene', 'anchor_index']).set_index(('id'))

            cls.df_V_ref_genome = cls.df_genomicVs.set_index('name').join(df_V_anchors.set_index('gene')).reset_index()
            cls.dict_genomicVs = (cls.df_genomicVs.set_index('name').to_dict())['value']
        except Exception as e:
            print('No V genes were found.')
            print(e)
            pass

        # J genes
        cls.df_genomicJs = pd.DataFrame.from_records(sqlrecords_genomicJs, columns=['id', 'name', 'value']).set_index('id')
        try:
            df_J_anchors = pd.DataFrame.from_records(sqlrecords_J_gene_CDR3_anchors, columns=['id', 'gene', 'anchor_index']).set_index(('id'))
            cls.df_J_ref_genome = cls.df_genomicJs.set_index('name').join(df_J_anchors.set_index('gene')).reset_index()
            cls.dict_genomicJs = (cls.df_genomicJs.set_index('name').to_dict())['value']
        except Exception as e:
            print('No J genes were found.')
            print(e)
            pass

        # D genes
        try:
            cls.df_genomicDs = pd.DataFrame.from_records(sqlrecords_genomicDs,
                                                         columns=['id', 'name', 'value']).set_index('id')
            cls.dict_genomicDs = (cls.df_genomicDs.set_index('name').to_dict())['value']
            # TODO: SHOULD I BE REBUNDANT? or df_genomicDs is rebundant?
            # self.df_D_ref_genome
        except Exception as e:
            print('No D genes were found.')
            print(e)
            pass

        return cls

    @classmethod
    def load_from_path(cls, path_ref_genome):
        cls = IgorRefGenome()
        cls.path_ref_genome = path_ref_genome
        cls.update_fln_names()
        cls.load_dataframes()
        return cls


    def update_fln_names(self, path_ref_genome=None, fln_genomicVs=None, fln_genomicDs=None, fln_genomicJs=None, fln_V_gene_CDR3_anchors=None, fln_J_gene_CDR3_anchors=None):
        if path_ref_genome is not None:
            self.path_ref_genome = path_ref_genome

        if fln_genomicVs is None:
            self.fln_genomicVs = self.path_ref_genome + "/" + "genomicVs.fasta"
        else:
            self.fln_genomicVs = fln_genomicVs

        if fln_genomicDs is None:
            self.fln_genomicDs = self.path_ref_genome + "/" + "genomicDs.fasta"
        else:
            self.fln_genomicDs = fln_genomicDs

        if fln_genomicJs is None:
            self.fln_genomicJs = self.path_ref_genome + "/" + "genomicJs.fasta"
        else:
            self.fln_genomicJs = fln_genomicJs

        if fln_V_gene_CDR3_anchors is None:
            self.fln_V_gene_CDR3_anchors = self.path_ref_genome + "/" + "V_gene_CDR3_anchors.csv"
        else:
            self.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors

        if fln_J_gene_CDR3_anchors is None:
            self.fln_J_gene_CDR3_anchors = self.path_ref_genome + "/" + "J_gene_CDR3_anchors.csv"
        else:
            self.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

    # TODO: LOAD INSTANCE FROM DEFINED FILES, what is the difference btwn load_dataframes?
    def load_from_files(self, fln_genomicVs=None, fln_genomicDs=None, fln_genomicJs=None,
                        fln_V_gene_CDR3_anchors=None, fln_J_gene_CDR3_anchors=None):
        self.fln_genomicVs = fln_genomicVs
        self.fln_genomicDs = fln_genomicDs
        self.fln_genomicJs = fln_genomicJs
        self.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
        self.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

        self.load_dataframes()

    def load_dataframes(self):
        # Fasta to dataframe
        # V genes
        self.df_genomicVs = from_fasta_to_dataframe(self.fln_genomicVs)
        try:
            df_V_anchors = pd.read_csv(self.fln_V_gene_CDR3_anchors, sep=';')

            self.df_V_ref_genome = self.df_genomicVs.set_index('name').join(df_V_anchors.set_index('gene')).reset_index()
            self.dict_genomicVs = (self.df_genomicVs.set_index('name').to_dict())['value']
        except Exception as e:
            print('No V genes were found.')
            print(e)
            pass

        # J genes
        self.df_genomicJs = from_fasta_to_dataframe(self.fln_genomicJs)
        try:
            df_J_anchors = pd.read_csv(self.fln_J_gene_CDR3_anchors, sep=';')
            self.df_J_ref_genome = self.df_genomicJs.set_index('name').join(df_J_anchors.set_index('gene')).reset_index()
            self.dict_genomicJs = (self.df_genomicJs.set_index('name').to_dict())['value']
        except Exception as e:
            print('No J genes were found.')
            print(e)
            pass

        # D genes
        try:
            self.df_genomicDs = from_fasta_to_dataframe(self.fln_genomicDs)
            self.dict_genomicDs = (self.df_genomicDs.set_index('name').to_dict())['value']
            # TODO: SHOULD I BE REBUNDANT? or df_genomicDs is rebundant?
            # self.df_D_ref_genome
        except Exception as e:
            print('No D genes were found.')
            print(e)
            pass

        #return df_V_ref_genome, df_J_ref_genome

    def get_anchors_dict(self):
        dict_anchor_index = dict()
        dict_anchor_index['V'] = self.df_V_ref_genome.set_index('name')['anchor_index'].to_dict()
        dict_anchor_index['J'] = self.df_J_ref_genome.set_index('name')['anchor_index'].to_dict()
        return dict_anchor_index

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

        self.anchor_in_read = None

    def __str__(self):
        return str(self.to_dict())

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
            cls.insertions   = eval (csvsplit[4].replace("{","[").replace("}","]"))
            cls.deletions    = eval (csvsplit[5].replace("{","[").replace("}","]"))
            cls.mismatches   = eval (csvsplit[6].replace("{","[").replace("}","]"))
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

class IgorGeneTemplate:
    def __init__(self):
        self.flnGene = None
        self.dataframe = None

    def get_sequence(self, gene_name):
        # TODO: use dataframe return sequence
        sequence = ""
        return sequence

### IGOR MODEL ####
class IgorModel:
    def __init__(self, model_parms_file=None, model_marginals_file=None):
        self.parms = IgorModel_Parms()
        self.marginals = IgorModel_Marginals()
        self.genomic_dataframe_dict = dict()
        self.xdata = dict()
        self.factors = list()
        self.metadata = dict()
        self.specie = ""
        self.chain = ""

        self.Pmarginal = dict()

        # FIXME: But since DB is in refactor keep it for the moment
        self.BestScenariosHeaderList = list() # This is a ordered list store the nicknames of events in the header of the file
        # should be only necessary if no database present


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

        self.sequence_construction_event_list = list()

    def __getitem__(self, key):
        return self.xdata[key]

    def __str__(self):
        return ".xdata" + str(self.get_events_nicknames_list())

    # TODO: finish this method to load model with default installed igor.
    @classmethod
    def load_default(cls, IgorSpecie, IgorChain, modelpath=None): #rcParams['paths.igor_models']):
        """
        :return IgorModel loaded with the default location for specie and chain
        """        
        # IGoR run parameters
        #IgorSpecie    = specie #"mouse"
        #IgorChain     = chain #"tcr_beta"
        if modelpath is None:
            try:
                modelpath = run_igor_datadir() + "/models"
            except Exception as e:
                print("ERROR: getting default igor datadir.", e)




        IgorModelPath = modelpath+"/"+IgorSpecie+"/"+IgorChain+"/"
        print("Loading default IGoR model from path : ", IgorModelPath)
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
    def load_from_parms_marginals_object(cls, mdl_parms, mdl_marginals):
        cls = IgorModel()
        cls.parms = mdl_parms
        cls.marginals = mdl_marginals
        cls.generate_xdata()
        return cls

    # FIXME:
    @classmethod
    def load_from_networkx(cls, IgorSpecie, IgorChain):
        """
        :return IgorModel loaded with the default location for specie and chain
        """
        cls = IgorModel(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
        return cls

    def generate_xdata(self):
        # TODO: CHANGE TO
        Event_Genechoice_List = ['v_choice', 'j_choice', 'd_gene']
        Event_Dinucl_List = ['vd_dinucl', 'dj_dinucl', 'vj_dinucl']
        Event_Insertion_List = ['vd_ins', 'dj_ins', 'vj_ins']
        Event_Deletion_List = ['v_3_del', 'j_5_del', 'd_3_del', 'd_5_del']
        
        for key in self.marginals.marginals_dict:
            event = self.parms.get_Event(key)

            if event.event_type == 'DinucMarkov':
            #if key in Event_Dinucl_List:
                self.xdata[key] = xr.DataArray(self.marginals.marginals_dict[key].reshape(4,4), \
                          dims=('x', 'y'))
                labels = self.parms.Event_dict[key]['value'].values

                strDim = 'x'
                self.xdata[key][strDim] = range(len(self.xdata[key][strDim]))
                strCoord = 'lbl__' + strDim
                self.xdata[key][strCoord] = (strDim, labels)
                strDim = 'y'
                self.xdata[key][strDim] = range(len(self.xdata[key][strDim]))
                strCoord = 'lbl__' + strDim
                self.xdata[key][strCoord] = (strDim, labels)

            else:
                self.xdata[key] = xr.DataArray(self.marginals.marginals_dict[key], \
                                                   dims=tuple(self.marginals.network_dict[key]))
                #print "key: ", key, self.xdata[key].dims

                for strDim in self.xdata[key].dims:
                    self.xdata[key][strDim] = range(len(self.xdata[key][strDim]))
                    if strDim in Event_Genechoice_List:
                        #print strDim
                        #labels = self.parms.Event_dict[strDim]['name'].map(genLabel).values # FIXME: use the exact name defined in model_parms
                        labels = self.parms.Event_dict[strDim]['name'].values
                        strCoord = 'lbl__'+strDim
                        self.xdata[key][strCoord] = (strDim, labels) # range(len(self.xdata[key][coord]))

                        sequences = self.parms.Event_dict[strDim]['value'].values
                        strCoord = 'seq__' + strDim
                        self.xdata[key][strCoord] = (strDim, sequences)

                    elif not (strDim in Event_Dinucl_List):
                        labels = self.parms.Event_dict[strDim]['value'].values
                        strCoord = 'lbl__'+strDim
                        self.xdata[key][strCoord] = (strDim, labels) # range(len(self.xdata[key][coord]))
                        # event = self.parms.get_Event(key)
                        # print(event.event_type)
                        # self.xdata[key].attrs["event_type"] = event.event_type
                        # self.xdata[key].attrs["seq_type"] = event.seq_type
                        # self.xdata[key].attrs["seq_side"] = event.seq_side

            # Event attributes
            self.xdata[key].attrs["nickname"] = event.nickname
            self.xdata[key].attrs["event_type"] = event.event_type
            self.xdata[key].attrs["seq_type"] = event.seq_type
            self.xdata[key].attrs["seq_side"] = event.seq_side
            self.xdata[key].attrs["priority"] = event.priority

            self.xdata[key].attrs["parents"] = list(self.parms.G.predecessors(key))
            self.xdata[key].attrs["childs"] = list(self.parms.G.successors(key))


        self.generate_Pmarginals()

    def get_zero_xarray_from_list(self, strEvents_list:list):
        #strEvents_list = ['v_choice', 'j_choice']
        strEvents_tuple = tuple(strEvents_list)

        # Use model parms to create xarray with values
        da_shape_list = [len(self.parms.Event_dict[str_event_nickname]) for str_event_nickname in strEvents_list]
        da_shape_tuple = tuple(da_shape_list)
        da = xr.DataArray(np.zeros(da_shape_tuple), dims=strEvents_tuple)

        for event_nickname in strEvents_list:
            da[event_nickname] = self.parms.Event_dict[event_nickname].index.values
            labels = self.parms.Event_dict[event_nickname]['name'].values
            strCoord = 'lbl__' + event_nickname
            da[strCoord] = (event_nickname, labels)

        return da

    def VE_get_Pmarginals_initial_factors(self):
        factors = list()
        for da in self.xdata.values():
            if da.attrs["event_type"] == 'DinucMarkov':
                sarray = da.stack(z=('x', 'y'))
                sarray = sarray.rename({"z": da.attrs["nickname"]})
                # FIXME: I'm removing DinucMarkov in the factors because
                #  P(vd_dinucl) = P(y|x) and we don't have P(x)
                #  So we can't marginalize P(y) neither P(x,y)
            else:
                factors.append(da)

        return factors

    # Doesn't need to be a self method, but ...
    def VE_get_factors_by_sum_out_variable(self, var_to_eliminate, factors):
        #     var_to_eliminate = 'j_choice'
        factors_to_sum_out = list()
        _factors = list()

        # separate factors to sum-out
        for factor in factors:
            # if factor.attrs["event_type"] == 'DinucMarkov':
            #     lista = [factor.attrs["nickname"]] + factor.attrs["parents"]
            #     var_intersection = {var_to_eliminate}.intersection(set(lista))
            # else:
            var_intersection = {var_to_eliminate}.intersection(set(factor.dims))
            # print("var_intersection : ", var_intersection)
            if len(var_intersection) > 0:
                # print(var_intersection)
                factors_to_sum_out.append(factor)
            else:
                _factors.append(factor)

        # sum-out-variables
        da_sum_over_event = 1
        for factor in factors_to_sum_out:
            da_sum_over_event = da_sum_over_event * factor

        new_factor = da_sum_over_event.sum(dim=var_to_eliminate)
        factors = _factors + [new_factor]
        return factors

    def VE_get_Pmarginal_of_event(self, strEvent):
        # FIXME: use xdata instead of self.parms
        sorted_events = self.parms.get_Event_list_sorted()
        sorted_events_to_marginalize = [event for event in sorted_events if not event.event_type == "DinucMarkov"]
        sorted_events_to_marginalize_without_VE = [event for event in sorted_events_to_marginalize if
                                                   not event.nickname == strEvent]

        # Start eliminating events
        factors = self.VE_get_Pmarginals_initial_factors()
        for event_to_eliminate_VE in sorted_events_to_marginalize_without_VE:
            factors = self.VE_get_factors_by_sum_out_variable(event_to_eliminate_VE.nickname, factors)

        # Now multiply the remaining factors to get the marginal.
        Pmarginal = 1
        for factor in factors:
            Pmarginal = Pmarginal * factor
        return Pmarginal

    def generate_Pmarginals(self):
        # Apply Variable elimination method for this
        # 1. Get a list of event sorted by high priority and less number of parents
        self.Pmarginal = dict()
        # Get the marginal for each event
        for key, darray in self.xdata.items():
            # print(key, darray)
            if darray.attrs["event_type"] == "DinucMarkov":
                self.Pmarginal[key] = darray
            else:
                self.Pmarginal[key] = self.VE_get_Pmarginal_of_event(key)


    # # FIXME: MAKE IT GENERAL
    # def generate_Pmarginals(self):
    #     # FIXME: GENERALIZE FOR ANY NETWORK NOT ONLY FOR VDJ AND VJ
    #     strEvent = 'v_choice'
    #     self.Pmarginal[strEvent] = self.xdata[strEvent]
    #
    #     strEvent = 'j_choice'
    #     strEventParent01 = 'v_choice'
    #     self.Pmarginal[strEvent] = self.xdata[strEvent].dot(self.xdata[strEventParent01])
    #
    #     strEvent = 'v_3_del'
    #     strEventParent01 = 'v_choice'
    #     Pjoint_aux = self.Pmarginal[strEventParent01]
    #     self.Pmarginal[strEvent] = self.xdata[strEvent].dot(Pjoint_aux)
    #
    #     strEvent = 'j_5_del'
    #     strEventParent01 = 'j_choice'
    #     Pjoint_aux = self.Pmarginal[strEventParent01]
    #     self.Pmarginal[strEvent] = self.xdata[strEvent].dot(Pjoint_aux)
    #
    #     if 'd_gene' in self.xdata.keys():
    #         strEvent = 'd_gene'
    #         strEventParent01 = 'v_choice'
    #         strEventParent02 = 'j_choice'
    #         Pjoint_aux = self.xdata[strEventParent02] * self.Pmarginal[strEventParent01]
    #         self.Pmarginal[strEvent] = self.xdata[strEvent].dot(Pjoint_aux)
    #
    #         strEvent = 'd_gene'
    #         strEventParent01 = 'v_choice'
    #         strEventParent02 = 'j_choice'
    #         Pjoint_aux = self.xdata[strEventParent02] * self.Pmarginal[strEventParent01]
    #         self.Pmarginal[strEvent] = self.xdata[strEvent].dot(Pjoint_aux)
    #
    #         strEvent = 'd_5_del'
    #         strEventParent01 = 'd_gene'
    #         Pjoint_aux = self.Pmarginal[strEventParent01]
    #         self.Pmarginal[strEvent] = self.xdata[strEvent].dot(Pjoint_aux)
    #
    #         strEvent = 'd_3_del'
    #         strEventParent01 = 'd_gene'
    #         strEventParent02 = 'd_5_del'
    #         Pjoint_aux = self.xdata[strEventParent02] * self.Pmarginal[strEventParent01]
    #         self.Pmarginal[strEvent] = self.xdata[strEvent].dot(Pjoint_aux)
    #
    #         strEvent = 'vd_ins'
    #         self.Pmarginal[strEvent] = self.xdata[strEvent]
    #
    #         strEvent = 'vd_dinucl'
    #         self.Pmarginal[strEvent] = self.xdata[strEvent]
    #
    #         strEvent = 'dj_ins'
    #         self.Pmarginal[strEvent] = self.xdata[strEvent]
    #
    #         strEvent = 'dj_dinucl'
    #         self.Pmarginal[strEvent] = self.xdata[strEvent]
    #
    #     else: # VJ NETWORK
    #         strEvent = 'vj_ins'
    #         self.Pmarginal[strEvent] = self.xdata[strEvent]
    #
    #         strEvent = 'vj_dinucl'
    #         self.Pmarginal[strEvent] = self.xdata[strEvent]

    def export_csv(self, fln_prefix, sep=';'):
        # FIXME: TEMPORARY SOLUTION FOR VERY PARTICULAR CASES.
        #################################################################################
        strEvent = 'v_choice'
        da = self.xdata[strEvent]
        # print(list(self.parms.G.predecessors(strEvent)))
        evento = self.parms.get_Event(strEvent)
        df = pd.DataFrame(data=da.values, index=da['lbl__' + strEvent].values,
                          columns=["P"])  # da['lbl__' + strEvent].values
        lbl_file = fln_prefix + "P__" + strEvent + ".csv"
        df.to_csv(lbl_file, index_label=evento.seq_type, sep=sep)

        ### v_3_del
        strEvent = 'v_3_del'
        da = self.xdata[strEvent]
        # print(list(self.parms.G.predecessors(strEvent)))
        parents = list(self.parms.G.predecessors(strEvent))
        evento = self.parms.get_Event(strEvent)

        dependencias = list(self.xdata[strEvent].dims)
        # print("********", dependencias, strEvent)
        dependencias.remove(strEvent)
        dependencias_dim = [self.xdata[strEvent][dep].shape[0] for dep in dependencias]

        if len(parents) == 1:
            df = pd.DataFrame(data=da.values, index=da['lbl__' + dependencias[0]].values,
                              columns=da['lbl__' + strEvent].values)
            lbl_file = fln_prefix + "P__" + strEvent + "__G__" + dependencias[0] + ".csv"
            df.to_csv(lbl_file, sep=sep)  # , index_label=evento.seq_type)

        #################################################################################
        strEvent = 'j_choice'
        da = self.xdata[strEvent]
        parents = list(self.parms.G.predecessors(strEvent))
        evento = self.parms.get_Event(strEvent)

        dependencias = list(self.xdata[strEvent].dims)
        #print("********", dependencias, strEvent)
        dependencias.remove(strEvent)
        dependencias_dim = [self.xdata[strEvent][dep].shape[0] for dep in dependencias]

        if len(parents) == 0:
            df = pd.DataFrame(data=da.values, index=da['lbl__' + strEvent].values,
                              columns=["P"])  # da['lbl__' + strEvent].values
            lbl_file = fln_prefix + "P__" + strEvent + ".csv"
            df.to_csv(lbl_file, index_label=evento.seq_type, sep=sep)
        elif len(parents) == 1:
            df = pd.DataFrame(data=da.values, index=da['lbl__' + dependencias[0]].values,
                              columns=da['lbl__' + strEvent].values)
            lbl_file = fln_prefix + "P__" + strEvent + "__G__" + dependencias[0] + ".csv"
            df.to_csv(lbl_file, sep=sep) #, index_label=evento.seq_type)
        else:
            print("Recombination event "+strEvent+" has an export problem!")

        ### j_5_del
        strEvent = 'j_5_del'
        da = self.xdata[strEvent]
        # print(list(self.parms.G.predecessors(strEvent)))
        parents = list(self.parms.G.predecessors(strEvent))
        evento = self.parms.get_Event(strEvent)

        dependencias = list(self.xdata[strEvent].dims)
        # print("********", dependencias, strEvent)
        dependencias.remove(strEvent)
        dependencias_dim = [self.xdata[strEvent][dep].shape[0] for dep in dependencias]

        if len(parents) == 1:
            df = pd.DataFrame(data=da.values, index=da['lbl__' + dependencias[0]].values,
                              columns=da['lbl__' + strEvent].values)
            lbl_file = fln_prefix + "P__" + strEvent + "__G__" + dependencias[0] + ".csv"
            df.to_csv(lbl_file, sep=sep)  # , index_label=evento.seq_type)

        #################################################################################

        if 'd_gene' in self.xdata.keys():

            strEvent = 'd_gene'
            da = self.xdata[strEvent]
            parents = list(self.parms.G.predecessors(strEvent))
            print(parents)
            evento = self.parms.get_Event(strEvent)
            print(evento.event_type)
            print(evento.seq_type)
            dependencias = list(self.xdata[strEvent].dims)
            print("********", dependencias, strEvent)
            dependencias.remove(strEvent)
            dependencias_dim = [self.xdata[strEvent][dep].shape[0] for dep in dependencias]

            if len(parents) == 0:
                df = pd.DataFrame(data=da.values, index=da['lbl__' + strEvent].values,
                                  columns=["P"])  # da['lbl__' + strEvent].values
                lbl_file = fln_prefix + "P__" + strEvent + ".csv"
                df.to_csv(lbl_file, index_label=evento.seq_type, sep=sep)
            elif len(parents) == 1:
                df = pd.DataFrame(data=da.values, index=da['lbl__' + dependencias[0]].values,
                              columns=da['lbl__' + strEvent].values)
                lbl_file = fln_prefix + "P__" + strEvent + "__G__" + dependencias[0] + ".csv"
                df.to_csv(lbl_file, sep=sep)  # , index_label=evento.seq_type)
            elif len(parents) == 2:
                lbl_file = fln_prefix + "P__" + strEvent + "__G__" + dependencias[0] + "__" + dependencias[1] + ".csv"
                with open(lbl_file, 'w') as ofile:
                    for ii in da[strEvent].values:
                        title = "P(" + da["lbl__"+strEvent].values[ii] +"| "+dependencias[0] + "," + dependencias[1]+")"
                        ofile.write("\n"+title+"\n")
                        da_ii = da[{strEvent: ii}]
                        df = pd.DataFrame(data=da_ii.values, index=da['lbl__' + dependencias[0]].values,
                                          columns=da['lbl__' + dependencias[1]].values)

                        df.to_csv(ofile, mode='a', sep=sep) # , index_label=evento.seq_type)
            else:
                print("Recombination event " + strEvent + " has an export problem!")

            strEvent = 'd_gene'
            da = self.xdata[strEvent]
            parents = list(self.parms.G.predecessors(strEvent))
            evento = self.parms.get_Event(strEvent)
            dependencias = list(self.xdata[strEvent].dims)
            dependencias.remove(strEvent)
            dependencias_dim = [self.xdata[strEvent][dep].shape[0] for dep in dependencias]

            if len(parents) == 0:
                df = pd.DataFrame(data=da.values, index=da['lbl__' + strEvent].values,
                                  columns=["P"])  # da['lbl__' + strEvent].values
                lbl_file = fln_prefix + "P__" + strEvent + ".csv"
                df.to_csv(lbl_file, index_label=evento.seq_type, sep=sep)
            elif len(parents) == 1:
                df = pd.DataFrame(data=da.values, index=da['lbl__' + dependencias[0]].values,
                                  columns=da['lbl__' + strEvent].values)
                lbl_file = fln_prefix + "P__" + strEvent + "__G__" + dependencias[0] + ".csv"
                df.to_csv(lbl_file, sep=sep)  # , index_label=evento.seq_type)
            elif len(parents) == 2:
                lbl_file = fln_prefix + "P__" + strEvent + "__G__" + dependencias[0] + "__" + dependencias[1] + ".csv"
                with open(lbl_file, 'w') as ofile:
                    for ii in da[strEvent].values:
                        title = "P(" + da["lbl__" + strEvent].values[ii] + "| " + dependencias[0] + "," + dependencias[1] + ")"
                        ofile.write("\n" + title + "\n")
                        da_ii = da[{strEvent: ii}]
                        df = pd.DataFrame(data=da_ii.values, index=da['lbl__' + dependencias[0]].values,
                                          columns=da['lbl__' + dependencias[1]].values)

                        df.to_csv(ofile, mode='a', sep=sep)  # , index_label=evento.seq_type)
            else:
                print("Recombination event " + strEvent + " has an export problem!")

            #return df

            ## P(D3, D5 | D) = P( D3| D5,D) x P (D5,D)
            #### Deletions in D
            da = self.xdata['d_3_del']*self.xdata['d_5_del']

            ### DELETIONS
            strEvent = 'd_gene'
            da = self.xdata[strEvent]
            dependencias = list(da.dims)
            print("********", dependencias, strEvent)
            dependencias.remove(strEvent)
            dependencias_dim = [da[dep].shape[0] for dep in dependencias]

            lbl_file = fln_prefix + "P__" + strEvent + "__deletions" + ".csv"
            with open(lbl_file, 'w') as ofile:
                for ii in da[strEvent].values:
                    da_ii = da[{strEvent: ii}]
                    lbl_event_realization = da['lbl__' + strEvent].values[ii]
                    title = "_P(" + dependencias[0] + "," + dependencias[1] + "| " + strEvent + " = " + lbl_event_realization + ")"
                    ofile.write(title + "\n")
                    df = pd.DataFrame(data=da_ii.values, index=da['lbl__' + dependencias[0]].values,
                                      columns=da['lbl__' + dependencias[1]].values)

                    df.to_csv(ofile, mode='a', sep=sep)  # , index_label=evento.seq_type)
                    ofile.write("\n")

            # self.xdata['d_3_del'] # P( D3| D5,D)

            ### INSERTIONS
            strEvent = 'vd_ins'
            da = self.xdata[strEvent]
            df_vd = pd.DataFrame(data=da.values, index=da['lbl__' + strEvent].values,
                                 columns=["P("+strEvent+")"])  # da['lbl__' + strEvent].values
            strEvent = 'dj_ins'
            da = self.xdata[strEvent]
            df_dj = pd.DataFrame(data=da.values, index=da['lbl__' + strEvent].values,
                                 columns=["P("+strEvent+")"])  # da['lbl__' + strEvent].values

            df = df_vd.merge(df_dj, left_index=True, right_index=True)
            lbl_file = fln_prefix + "P__" + "insertions" + ".csv"
            df.to_csv(lbl_file, index_label="Insertions", sep=sep)

            ### DINUCL
            strEvent = 'vd_dinucl'
            da = self.xdata[strEvent]
            print(da)
            df = pd.DataFrame(data=da.values, index=da['lbl__x'].values,
                                 columns=da['lbl__y'].values)
            lbl_file = fln_prefix + "P__" + strEvent + ".csv"
            df.to_csv(lbl_file, index_label="From\To", sep=sep)

            strEvent = 'dj_dinucl'
            da = self.xdata[strEvent]
            print(da)
            df = pd.DataFrame(data=da.values, index=da['lbl__x'].values,
                              columns=da['lbl__y'].values)
            lbl_file = fln_prefix + "P__" + strEvent + ".csv"
            df.to_csv(lbl_file, index_label="From\To", sep=sep)

        else:
            ### INSERTIONS
            strEvent = 'vj_ins'
            da = self.xdata[strEvent]
            df = pd.DataFrame(data=da.values, index=da['lbl__' + strEvent].values,
                                 columns=["P(" + strEvent + ")"])  # da['lbl__' + strEvent].values


            lbl_file = fln_prefix + "P__" + "insertions" + ".csv"
            df.to_csv(lbl_file, index_label="Insertions", sep=sep)

            ### DINUCL
            strEvent = 'vj_dinucl'
            da = self.xdata[strEvent]
            print(da)
            df = pd.DataFrame(data=da.values, index=da['lbl__x'].values,
                              columns=da['lbl__y'].values)
            lbl_file = fln_prefix + "P__" + strEvent + ".csv"
            df.to_csv(lbl_file, index_label="From\To", sep=sep)

        #################################################################################

        # strEvent = 'j_choice'
        # strEventParent01 = 'v_choice'
        # self.Pmarginal[strEvent] = self.xdata[strEvent].dot(self.xdata[strEventParent01])
        #
        # strEvent = 'v_3_del'
        # strEventParent01 = 'v_choice'
        # Pjoint_aux = self.Pmarginal[strEventParent01]
        # self.Pmarginal[strEvent] = self.xdata[strEvent].dot(Pjoint_aux)

    # FIXME: THIS METHOD IS NOT FINISH!!
    def export_event_to_csv(self, event_nickname, fln_prefix, sep=';'):
        # if kwargs.get('sep') is None:
        #     kwargs['sep'] = ';'

        da = self.xdata[event_nickname]
        event = self.parms.get_Event(event_nickname)
        if da.event_type == 'GeneChoice':
            # print(list(self.parms.G.predecessors(strEvent)))
            df = pd.DataFrame(data=da.values, index=da['lbl__' + event_nickname].values,
                              columns=["P"])  # da['lbl__' + strEvent].values
            lbl_file = fln_prefix + "P__" + event_nickname + ".csv"
            df.to_csv(lbl_file, index_label=event.seq_type, sep=sep)

    def export_Pmarginal_to_csv(self, event_nickname:str, *args, **kwargs):

        if kwargs.get('sep') is None:
            kwargs['sep'] = ';'

        event = self.parms.get_Event(event_nickname, by_nickname=True)
        da = self.xdata[event_nickname]
        if event.event_type == 'GeneChoice' :
            df = da.to_dataframe(name="P") #.drop('priority', 1)
            df.to_csv(*args, **kwargs)
        elif event.event_type == 'Insertion':
            df = da.to_dataframe(name="P")  # .drop('priority', 1)
            df.to_csv(*args, **kwargs)
        elif event.event_type == 'Deletion':
            df = da.to_dataframe(name="P")  # .drop('priority', 1)
            df.to_csv(*args, **kwargs)

        elif event.event_type == 'DinucMarkov':
            ### FIXME:
            strEvent = 'dj_dinucl'
            df = pd.DataFrame(data=da.values, index=da['lbl__x'].values,
                              columns=da['lbl__y'].values)
            kwargs['index_label'] = "From\To"
            df.to_csv(*args, **kwargs) #, index_label=, sep=sep)

        else:
            print("Event nickname "+event_nickname+" is not present in this model.")
            print("Accepted Events nicknames are : "+str(self.get_events_nicknames_list()))

    # FIXME: CHANGE EVENT MARGINAL!!!
    def get_Event_Marginal(self, event_nickname: str):
        """Returns an xarray with the marginal probability of the event given the nickname"""
        # FIXME: add new way to make the recursion.
        # FIXME: FIRST without recursion, return
        if event_nickname in self.parms.get_EventsNickname_list():
            da_event = self.xdata[event_nickname]
            dependencies = self.parms.Edges_dict[event_nickname]
            event_parents = list( self.parms.G.predecessors(event_nickname) )
            #1. Sort the dependencies by priority, then by dependencie
            # if queue is not empty:
            #     self.get_Event_Marginal(nicki)
            if event_nickname == 'v_choice':
                da_marginal = da_event
                return da_marginal
            elif event_nickname == 'j_choice':
                print(event_parents)
                strEventParent = 'v_choice'
                da_event_parent = self.xdata[strEventParent]
                da_marginal = da_event_parent.dot(da_event)
                da_parents = 1
                for parent in event_parents:
                    da_parents = da_parents * self.xdata[parent]
                print('&'*20)
                print(da_event.dot(da_parents))
                return da_marginal
                # return da_event, da_event_parent, (da_event_parent.dot(da_event)), da_marginal.sum()
            elif event_nickname == 'd_gene':

                print(event_parents)
                strEventParent = 'v_choice'
                strEventParent2 = 'j_choice'
                da_event_parent = self.xdata[strEventParent]
                da_event_parent2 = self.xdata[strEventParent2]
                # da_marginal = da_event.dot()
                da_marginal = da_event_parent*da_event_parent2
                return da_marginal

        else:
            print("Event nickname : " + event_nickname + " is not an event in this IGoR model.")
            return list()

    def plot_event_GeneChoice(self, event_nickname:str, **kwargs):
        """ Return GeneChoice plot """
        # Default values in plot

        import numpy as np
        v_genLabel = np.vectorize(genLabel)
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        da = self.xdata[event_nickname]
        for parent_nickname in da.dims:  # attrs['parents']
            da["lbl__" + parent_nickname].values = v_genLabel(da["lbl__" + parent_nickname].values)

        parents_list = da.attrs['parents']
        if len(parents_list) == 0:
            # ONE DIMENSIONAL PLOT
            titulo = "$P($" + event_nickname + "$)$"
            fig, ax = plt.subplots(figsize=(18, 15))
            XX = da[event_nickname].values
            YY = da.values
            ax.bar(XX, YY, **kwargs)
            lbl_XX = da['lbl__' + event_nickname].values
            ax.set_xticks(XX)
            ax.set_xticklabels(v_genLabel(lbl_XX), rotation=90)
            ax.set_title(titulo)
            return fig, ax
        elif len(parents_list) == 1:
            if not 'cmap' in kwargs.keys():
                kwargs['cmap'] = 'gnuplot2_r'
            lbl_parents = ",".join(da.attrs['parents'])
            titulo = "$P($" + event_nickname + "$|$" + lbl_parents + "$)$"
            fig, ax = plt.subplots(figsize=(18, 15))
            XX = da[event_nickname].values
            YY = da.values
            da.plot(ax=ax, **kwargs)
            lbl_XX = da['lbl__' + event_nickname].values
            ax.set_title(titulo)
            ax.set_aspect('equal')
            return fig, ax
        elif len(parents_list) == 2:
            if not 'cmap' in kwargs.keys():
                kwargs['cmap'] = 'gnuplot2_r'
            # da = self.xdata[event_nickname]
            fig, ax = plt.subplots(*da[event_nickname].shape, figsize=(10, 20))
            for ii, ev_realiz in enumerate(da[event_nickname]):
                # print(ev_realiz.values)
                da[{event_nickname: ev_realiz.values}].plot(ax=ax[ii], cmap='gnuplot2_r')
                lbl_ev_realiz = str( ev_realiz["lbl__" + event_nickname].values )
                lbl_parents = str( ",".join(da.attrs['parents']) )
                titulo = "$P($" + event_nickname + "$ = $ " + lbl_ev_realiz + " $|$" + lbl_parents + "$)$"
                ax[ii].set_title(titulo)
            return fig, ax
        else:
            fig, ax = plt.subplots()
            ax.set_title("Dimensionality not supportted for event : ", event_nickname)
            return fig, ax

    def plot_event_Insertion(self, event_nickname:str, **kwargs):
        """ Return Insertion plot """

        import numpy as np
        v_genLabel = np.vectorize(genLabel)
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        da = self.xdata[event_nickname]

        parents_list = da.attrs['parents']
        if len(parents_list) == 0:
            # ONE DIMENSIONAL PLOT
            titulo = "$P($" + event_nickname + "$)$"
            fig, ax = plt.subplots(figsize=(18, 15))
            XX = da[event_nickname].values
            YY = da.values
            ax.bar(XX, YY, **kwargs)
            lbl_XX = da['lbl__' + event_nickname].values
            ax.set_xticks(XX)
            ax.set_xticklabels(lbl_XX, rotation=90)
            ax.set_title(titulo)
            return fig, ax
        elif len(parents_list) == 1:
            titulo = "$P($" + event_nickname + "$|$" + ",".join(da.attrs['parents']) + "$)$"
            fig, ax = plt.subplots(figsize=(18, 15))
            XX = da[event_nickname].values
            YY = da.values
            da.plot(ax=ax, **kwargs)
            lbl_XX = da['lbl__' + event_nickname].values
            ax.set_title(titulo)
            ax.set_aspect('equal')
            return fig, ax
        elif len(parents_list) == 2:
            # da = self.xdata[event_nickname]
            fig, ax = plt.subplots(*da[event_nickname].shape, figsize=(10, 20))
            for ii, ev_realiz in enumerate(da[event_nickname]):
                # print(ev_realiz.values)
                da[{event_nickname: ev_realiz.values}].plot(ax=ax[ii], cmap='gnuplot2_r')
                lbl_ev_realiz = str(ev_realiz["lbl__" + event_nickname].values)
                lbl_parents = str(",".join(da.attrs['parents']))
                print(lbl_ev_realiz, lbl_parents)
                titulo = "$P($" + event_nickname + "$ = $ " + lbl_ev_realiz + " $|$" + lbl_parents + "$)$"
                ax[ii].set_title(titulo)
            return fig, ax
        else:
            fig, ax = plt.subplots()
            ax.set_title("Dimensionality not supportted for event : ", event_nickname)
            return fig, ax

    def plot_event_Deletion(self, event_nickname:str, **kwargs):
        """ Return GeneChoice plot """

        import numpy as np
        v_genLabel = np.vectorize(genLabel)
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        da = self.xdata[event_nickname]

        parents_list = da.attrs['parents']
        if len(parents_list) == 0:
            # ONE DIMENSIONAL PLOT
            titulo = "$P($" + event_nickname + "$)$"
            fig, ax = plt.subplots(figsize=(18, 15))
            XX = da[event_nickname].values
            YY = da.values
            ax.bar(XX, YY, **kwargs)
            lbl_XX = da['lbl__' + event_nickname].values
            ax.set_xticks(XX)
            ax.set_xticklabels(lbl_XX, rotation=90)
            ax.set_title(titulo)
            return fig, ax
        elif len(parents_list) == 1:
            if not 'cmap' in kwargs.keys():
                kwargs['cmap'] = 'gnuplot2_r'
            titulo = "$P($" + event_nickname + "$|$" + ",".join(da.attrs['parents']) + "$)$"
            fig, ax = plt.subplots(figsize=(18, 15))
            XX = da[event_nickname].values
            YY = da.values
            da.plot(ax=ax, **kwargs)
            lbl_XX = da['lbl__' + event_nickname].values
            ax.set_title(titulo)
            ax.set_aspect('equal')
            return fig, ax
        elif len(parents_list) == 2:
            if not 'cmap' in kwargs.keys():
                kwargs['cmap'] = 'gnuplot2_r'
            # da = self.xdata[event_nickname]
            fig, ax = plt.subplots(*da[event_nickname].shape, figsize=(10, 50))
            for ii, ev_realiz in enumerate(da[event_nickname]):
                # print(ev_realiz.values)
                lbl_ev_realiz = str(ev_realiz["lbl__" + event_nickname].values)
                lbl_parents = str(",".join(da.attrs['parents']))
                da[{event_nickname: ev_realiz.values}].plot(ax=ax[ii], cmap='gnuplot2_r')
                titulo = "$P($" + event_nickname + "$ = $ " + lbl_ev_realiz + " $|$" + lbl_parents + "$)$"
                ax[ii].set_title(titulo)
            return fig, ax
        else:
            fig, ax = plt.subplots()
            ax.set_title("Dimensionality not supportted for event : ", event_nickname)
            return fig, ax

    def plot_event_DinucMarkov(self, event_nickname:str, **kwargs):
        """ Return GeneChoice plot """
        # Default values in plot
        if not 'cmap' in kwargs.keys():
            kwargs['cmap'] = 'gnuplot2_r'

        import numpy as np
        import matplotlib.pyplot as plt
        da = self.xdata[event_nickname]
        lblEvent = event_nickname.replace("_", " ")
        xEtiqueta = lblEvent
        yEtiqueta = "P"
        fig, ax = plt.subplots()
        XX = da['x'].values
        YY = da['y'].values

        lbl__XX = da['lbl__' + 'x'].values
        lbl__YY = da['lbl__' + 'y'].values

        ZZ = da.values

        da.plot(ax=ax, x='x', y='y', vmin=0, vmax=1, **kwargs)

        ax.set_xlabel('From')
        ax.set_xticks(XX)
        ax.set_xticklabels(lbl__XX, rotation=0)

        ax.set_ylabel('To')
        ax.set_yticks(YY)
        ax.set_yticklabels(lbl__YY, rotation=0)

        ax.set_title(lblEvent)
        ax.set_aspect('equal')
        for i, j in zip(*ZZ.nonzero()):
            ax.text(j, i, ZZ[i, j], color='white', ha='center', va='center')

        return fig, ax

    def export_plot_events(self, outfilename_prefix):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        with PdfPages(outfilename_prefix + ".pdf") as pdf_file:
            fig, ax = plt.subplots()
            self.parms.plot_Graph(ax=ax)
            fig.tight_layout()
            pdf_file.savefig(fig)

            # GeneChoice, Insertion, Deletion, DinucMarkov
            for event_nickname in self.xdata.keys():
                event = self.parms.get_Event(event_nickname)
                if event.event_type == 'GeneChoice':
                    fig, ax = self.plot_event_GeneChoice(event_nickname)
                    fig.tight_layout()
                    pdf_file.savefig(fig)
                    del fig
                elif event.event_type == 'Insertion':
                    fig, ax = self.plot_event_Insertion(event_nickname)
                    fig.tight_layout()
                    pdf_file.savefig(fig)
                    del fig
                elif event.event_type == 'Deletion':
                    fig, ax = self.plot_event_Deletion(event_nickname)
                    fig.tight_layout()
                    pdf_file.savefig(fig)
                    del fig
                elif event.event_type == 'DinucMarkov':
                    fig, ax = self.plot_event_DinucMarkov(event_nickname)
                    fig.tight_layout()
                    pdf_file.savefig(fig)
                    del fig
                else:
                    print("ERROR: EVENT NOT RECOGNIZE", event_nickname)


    def export_plot_Pmarginals(self, outfilename_prefix):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        with PdfPages(outfilename_prefix+".pdf") as pdf_file:
            fig, ax = plt.subplots()
            self.parms.plot_Graph(ax=ax)
            fig.tight_layout()
            pdf_file.savefig(fig)

            for event_nickname in self.Pmarginal.keys():
                fig, ax = plt.subplots(figsize=(20, 10))
                self.plot_Event_Marginal(event_nickname, ax=ax)
                fig.tight_layout()
                # flnOutput = flnPrefix + "_" + event_nickname + ".pdf"
                pdf_file.savefig(fig)
                # fig.savefig(flnOutput)

    def plot_Event_Marginal(self, event_nickname:str, ax=None, **kwargs):
        """
        Plot marginals of model events by nickname
        """
        event = self.parms.get_Event(event_nickname, by_nickname=True)
        da = self.Pmarginal[event_nickname] # real marginal DataArray

        lblEvent = event_nickname.replace("_", " ")
        xEtiqueta = lblEvent
        yEtiqueta = "P"

        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()

        ax.set_xlabel(xEtiqueta)
        ax.set_ylabel(yEtiqueta, rotation=0)

        if event.event_type == 'GeneChoice' :
            # Bar plot
            XX = da[event_nickname].values
            YY = da.values

            ax.bar(XX, YY, **kwargs)

            lbl_XX = da['lbl__' + event_nickname].values
            ax.set_xticks(XX)
            ax.set_xticklabels(v_genLabel(lbl_XX), rotation=90)
            #return ax

        elif event.event_type == 'Insertion':
            # Use labels as a coordinate.
            # Insertions are in principle independent,
            # FIXME: but if not what to do.
            XX = da['lbl__' + event_nickname].values
            YY = da.values
            if not 'marker' in kwargs.keys():
                kwargs['marker'] = 'o'
            ax.plot(XX, YY, **kwargs)

        elif event.event_type == 'Deletion':
            #YY = self.xdata[event_nickname].values
            #XX = self.xdata[event_nickname]['lbl__' + event_nickname].values
            #ax.plot(XX, YY)
            XX = da['lbl__' + event_nickname].values
            YY = da.values
            if not 'marker' in kwargs.keys():
                kwargs['marker'] = 's'
            ax.plot(XX, YY, **kwargs)

        elif event.event_type == 'DinucMarkov':
            XX = da['x'].values
            YY = da['y'].values

            lbl__XX = da['lbl__' + 'x'].values
            lbl__YY = da['lbl__' + 'y'].values

            ZZ = da.values

            da.plot(ax=ax, x='x', y='y', vmin=0, vmax=1, cmap='gnuplot2_r', **kwargs)

            ax.set_xlabel('From')
            ax.set_xticks(XX)
            ax.set_xticklabels(lbl__XX, rotation=0)

            ax.set_ylabel('To')
            ax.set_yticks(YY)
            ax.set_yticklabels(lbl__YY, rotation=0)

            ax.set_title(lblEvent)
            ax.set_aspect('equal')
            for i, j in zip(*ZZ.nonzero()):
                ax.text(j, i, ZZ[i, j], color='white', ha='center', va='center')


        else:
            print("Event nickname "+event_nickname+" is not present in this model.")
            print("Accepted Events nicknames are : "+str(self.get_events_nicknames_list()))

        #return self.get_Event_Marginal(nickname)
        # ax.set_title(lblEvent)
        return ax

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

    # PLOTS:
    def plot_Bayes_network(self, filename=None):
        if filename is None:
            return self.parms.plot_Graph()
        else:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax_ = self.parms.plot_Graph(ax=ax)
            fig.savefig(filename)
            return ax_

    def plot(self, event_nickname:str, ax=None):
        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()

        da = self.xdata[event_nickname]


        return ax

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

    def set_genomic_dataframe_dict(self, dataframe_dict):
        self.genomic_dataframe_dict = dataframe_dict

    def scenario_from_database(self, scenarios_list):
        scen = scenarios_list[0]
        scen_dict = scen.realizations_ids_dict

        for event in self.parms.Event_list:
            if not (event.event_type == 'DinucMarkov'):
                scen_dict[event.nickname] = scen_dict.pop('id_' + event.nickname)

    # FIXME:
    def export_model(self, model_parms_file=None, model_marginals_file=None):
        self.parms.write_model_parms(filename=model_parms_file)
        self.marginals
        self.xdata

        print("Exporting model to ")

    def get_event_realizations_DataFrame(self, event_nickname):
        return self.parms.Event_dict[event_nickname]

    def get_event_realization_of_event(self, event_nickname, event_id):
        if type(event_id) is list:
            return list( map( lambda x: self.parms.get_Event(event_nickname).realizations[x], event_id) )
        else:
            return self.parms.get_Event(event_nickname).realizations[event_id]

    def get_realizations_dict_from_scenario_dict(self, scenario_realization_dict:dict):
        realization_dict = dict()
        # print(scenario_realization_dict)
        for event_nickname, event_id in scenario_realization_dict.items():
            if not ( event_nickname == 'mismatcheslen' or event_nickname == 'mismatches' or event_nickname == 'errors')  :
                realization_dict[event_nickname] = self.get_event_realization_of_event(event_nickname, event_id)

        return realization_dict

    def set_realization_event_from_DataFrame(self, event_nickname, new_df):
        self.parms.set_event_realizations_from_DataFrame(event_nickname, new_df)
        self.marginals.initialize_uniform_from_model_parms(self.parms)
        self.generate_xdata()

    def write_model(self, fln_model_parms, fln_model_marginals):
        self.parms.write_model_parms(filename=fln_model_parms)
        self.marginals.write_model_marginals(filename=fln_model_marginals, model_parms=self.parms)

    # FIXME: DEPRECATED
    # FIXME: Find a better way to get the order to construct a sequence.
    def generate_sequence_construction_list(self):
        """ Generate the list of events to reconstruct a sequence from an scenario self.sequence_construction_event_list """

        sequence_arrengement_dict = dict()
        # 1. 'V_gene'
        V_gene_list = [event for event in self.parms.Event_list if event.seq_type == 'V_gene']

        # 2. 'D_gene'
        D_gene_list = [event for event in self.parms.Event_list if event.seq_type == 'D_gene']

        # 3. 'J_gene'
        J_gene_list = [event for event in self.parms.Event_list if event.seq_type == 'J_gene']

        V_gene_list = sorted(V_gene_list,
                             key=lambda event: 100 * event.priority - len(self.xdata[event.nickname].attrs['parents']),
                             reverse=True)

        J_gene_list = sorted(J_gene_list,
                             key=lambda event: 100 * event.priority - len(self.xdata[event.nickname].attrs['parents']),
                             reverse=True)

        sequence_arrengement_dict['V_gene'] = V_gene_list
        sequence_arrengement_dict['J_gene'] = J_gene_list
        # since d_3_del and d_5_del have the same priority then
        arrengement_list = list()
        if len(D_gene_list) == 0:
            VJ_gene_list = [event for event in self.parms.Event_list if event.seq_type == 'VJ_gene']
            VJ_gene_list = sorted(VJ_gene_list,
                                 key=lambda event: 100 * event.priority - len(self.xdata[event.nickname].attrs['parents']),
                                 reverse=True)
            sequence_arrengement_dict['VJ_gene'] = VJ_gene_list
            arrengement_list = V_gene_list + VJ_gene_list + J_gene_list
        else:
            D_gene_list = sorted(D_gene_list,
                                 key=lambda event: 100 * event.priority - len(self.xdata[event.nickname].attrs['parents']),
                                 reverse=True)
            sequence_arrengement_dict['D_gene'] = D_gene_list
            VD_gene_list = [event for event in self.parms.Event_list if event.seq_type == 'VD_genes']
            VD_gene_list = sorted(VD_gene_list,
                                  key=lambda event: 100 * event.priority - len(self.xdata[event.nickname].attrs['parents']),
                                  reverse=True)
            sequence_arrengement_dict['VD_gene'] = VD_gene_list

            DJ_gene_list = [event for event in self.parms.Event_list if event.seq_type == 'DJ_gene']
            DJ_gene_list = sorted(DJ_gene_list,
                                  key=lambda event: 100 * event.priority - len(self.xdata[event.nickname].attrs['parents']),
                                  reverse=True)
            sequence_arrengement_dict['DJ_gene'] = VD_gene_list

            arrengement_list = V_gene_list + VD_gene_list + D_gene_list + DJ_gene_list + J_gene_list




        self.sequence_construction_event_list = arrengement_list
        return sequence_arrengement_dict # arrengement_list


    def construct_sequence_VDJ_from_realization_dict(self, scen_realization_dict):
        """return VDJ gene segment, which are the gene with the deletions of palindromic insertions"""
        # print("scen_realization_dict : ", scen_realization_dict)
        V_segment_dict = get_gene_segment(scen_realization_dict['v_choice'].value,
                                     int_gene_3_del=scen_realization_dict['v_3_del'].value)

        D_segment_dict = get_gene_segment(scen_realization_dict['d_gene'].value,
                                     int_gene_5_del=scen_realization_dict['d_5_del'].value,
                                     int_gene_3_del=scen_realization_dict['d_3_del'].value)

        J_segment_dict = get_gene_segment(scen_realization_dict['j_choice'].value,
                                     int_gene_5_del=scen_realization_dict['j_5_del'].value)

        VD_segment_dict = collections.OrderedDict()
        DJ_segment_dict = collections.OrderedDict()
        VD_segment_dict['gene_segment'] = "".join([realiz.value for realiz in scen_realization_dict['vd_dinucl']])
        DJ_segment_dict['gene_segment'] = "".join([realiz.value for realiz in scen_realization_dict['dj_dinucl']][::-1])


        return V_segment_dict, VD_segment_dict, D_segment_dict, DJ_segment_dict, J_segment_dict

    def construct_sequence_VJ_from_realization_dict(self, scen_realization_dict):
        """return VJ sequence, which are the gene with the deletions of palindromic insertions"""

        # print("scen_realization_dict : ", scen_realization_dict)
        V_segment_dict = get_gene_segment(scen_realization_dict['v_choice'].value,
                                     int_gene_3_del=scen_realization_dict['v_3_del'].value)
        J_segment_dict = get_gene_segment(scen_realization_dict['j_choice'].value,
                                     int_gene_5_del=scen_realization_dict['j_5_del'].value)

        VJ_segment_dict = collections.OrderedDict()
        VJ_segment_dict['gene_segment'] = "".join([realiz.value for realiz in scen_realization_dict['vj_dinucl']])




        return V_segment_dict, VJ_segment_dict, J_segment_dict


    # TODO: MAKE A METHOD TO EXPORT A LINE FROM AN SCENARIO
    def get_AIRR_VDJ_rearragement_dict_from_scenario(self, scenario, str_sequence, v_offset=0, pgen=None, junction=None, junction_aa=None):
        # get_AIRR_VDJ_rearragement_dict_from_scenario(scenario, indexed_seq.seq_index, indexed_seq.sequence)
        # airr_dict = dict()

        from .AIRR import AIRR_VDJ_rearrangement

        realizations_ids_dict = scenario.realizations_ids_dict
        realization_dict = self.get_realizations_dict_from_scenario_dict(realizations_ids_dict)

        v_segment, vd_segment, d_segment, dj_segment, j_segment = self.construct_sequence_VDJ_from_realization_dict(realization_dict)

        airr_vdj = AIRR_VDJ_rearrangement(sequence_id=scenario.seq_index, sequence=str_sequence)

        airr_vdj.v_data.call = realization_dict['v_choice'].name
        airr_vdj.d_data.call = realization_dict['d_gene'].name
        airr_vdj.j_data.call = realization_dict['j_choice'].name

        airr_vdj.sequence_alignment = v_segment['gene_segment'] + vd_segment['gene_segment'] + d_segment['gene_segment'] + dj_segment['gene_segment'] + j_segment['gene_segment']

        airr_vdj.np1 = v_segment['palindrome_3_end'] + vd_segment['gene_segment']
        airr_vdj.np2 = dj_segment['gene_segment'] + j_segment['palindrome_5_end']

        airr_vdj.pgen = pgen

        airr_vdj.junction = junction
        airr_vdj.junction_aa = junction_aa

        airr_vdj.rev_comp = False

        # FIXME: CORRECT CIGAR FORMAT TEMPORARY SOLUTION
        airr_vdj.v_data.cigar = str(len(v_segment['gene_cut']))+"M"
        airr_vdj.d_data.cigar = str(len(d_segment['gene_cut'])) + "M"
        airr_vdj.j_data.cigar = str(len(j_segment['gene_cut'])) + "M"

        airr_vdj.v_data.score = 5 * len(v_segment['gene_cut'])
        airr_vdj.d_data.score = 5 * len(d_segment['gene_cut'])
        airr_vdj.j_data.score = 5 * len(j_segment['gene_cut'])



        # V
        airr_vdj.v_data.sequence_start = 1
        airr_vdj.v_data.sequence_end = len(v_segment['palindrome_5_end']) + len(v_segment['gene_cut'])
        airr_vdj.v_data.germline_start = airr_vdj.v_data.sequence_start - v_offset - 1
        airr_vdj.v_data.germline_end = airr_vdj.v_data.sequence_end - airr_vdj.v_data.sequence_start - 1
        # = airr_vdj.v_data.germline_start + len(v_segment['palindrome_5_end']) + len(v_segment['gene_cut'])
        airr_vdj.p3v_length = len(v_segment['palindrome_3_end'])

        # VD
        airr_vdj.n1_length = realization_dict['vd_ins'].value
        airr_vdj.np1_length = airr_vdj.p3v_length + airr_vdj.n1_length
        airr_vdj.np1 = vd_segment['gene_segment'] # This include the palindromic insertions

        # D
        airr_vdj.p5d_length = len(d_segment['palindrome_5_end'])
        airr_vdj.d_data.germline_start = d_segment['gene_ini'] + 1
        airr_vdj.d_data.germline_end = d_segment['gene_end'] + 1
        airr_vdj.d_data.sequence_start = airr_vdj.np1_length + (airr_vdj.v_data.sequence_end - airr_vdj.v_data.sequence_start - 1 )
        airr_vdj.d_data.sequence_end = airr_vdj.d_data.sequence_start + len(d_segment['gene_cut']) - 1
        airr_vdj.p3d_length = len(d_segment['palindrome_3_end'])

        # DJ
        airr_vdj.n2_length = realization_dict['dj_ins'].value
        airr_vdj.np2_length =  airr_vdj.p5d_length + airr_vdj.n2_length + airr_vdj.p3d_length
        airr_vdj.np2 = dj_segment['gene_segment'] # This include the palindromic insertions

        # J
        airr_vdj.p5j_length = len(j_segment['palindrome_5_end'])
        airr_vdj.j_data.germline_start = j_segment['gene_ini'] + 1
        airr_vdj.j_data.germline_end = j_segment['gene_end'] + 1
        airr_vdj.j_data.sequence_start = airr_vdj.np2_length + (airr_vdj.d_data.sequence_end - airr_vdj.d_data.sequence_start - 1)
        airr_vdj.j_data.sequence_end = airr_vdj.j_data.sequence_start + len(j_segment['gene_cut']) - 1

        return airr_vdj.to_dict()


    def get_AIRR_VJ_rearragement_dict_from_scenario(self, scenario, str_sequence, v_offset=0, pgen=None, junction=None, junction_aa=None):
        """
        Return airr rearragement from scenario.
        """
        # get_AIRR_VDJ_rearragement_dict_from_scenario(scenario, indexed_seq.seq_index, indexed_seq.sequence)
        # airr_dict = dict()

        from .AIRR import AIRR_VDJ_rearrangement

        realizations_ids_dict = scenario.realizations_ids_dict
        realization_dict = self.get_realizations_dict_from_scenario_dict(realizations_ids_dict)

        # FIXME: HERE
        v_segment, vj_segment, j_segment = self.construct_sequence_VJ_from_realization_dict(realization_dict)

        airr_vj = AIRR_VDJ_rearrangement(sequence_id=scenario.seq_index, sequence=str_sequence)

        airr_vj.v_data.call = realization_dict['v_choice'].name
        airr_vj.d_data.call = None #realization_dict['d_gene'].name
        airr_vj.j_data.call = realization_dict['j_choice'].name

        airr_vj.sequence_alignment = v_segment['gene_segment'] + vj_segment['gene_segment'] + j_segment['gene_segment']

        airr_vj.np1 = v_segment['palindrome_3_end'] + vj_segment['gene_segment'] + j_segment['palindrome_5_end']
        airr_vj.np2 = None

        airr_vj.pgen = pgen

        airr_vj.junction = junction
        airr_vj.junction_aa = junction_aa

        airr_vj.rev_comp = False

        # FIXME: CORRECT CIGAR FORMAT TEMPORARY SOLUTION
        airr_vj.v_data.cigar = str(len(v_segment['gene_cut']))+"M"
        airr_vj.d_data.cigar = None #str(len(d_segment['gene_cut'])) + "M"
        airr_vj.j_data.cigar = str(len(j_segment['gene_cut'])) + "M"

        airr_vj.v_data.score = 5 * len(v_segment['gene_cut'])
        airr_vj.d_data.score = None #5 * len(d_segment['gene_cut'])
        airr_vj.j_data.score = 5 * len(j_segment['gene_cut'])



        # V
        airr_vj.v_data.sequence_start = 1
        airr_vj.v_data.sequence_end = len(v_segment['palindrome_5_end']) + len(v_segment['gene_cut'])
        airr_vj.v_data.germline_start = airr_vj.v_data.sequence_start - v_offset - 1
        airr_vj.v_data.germline_end = airr_vj.v_data.sequence_end - airr_vj.v_data.sequence_start - 1
        # = airr_vdj.v_data.germline_start + len(v_segment['palindrome_5_end']) + len(v_segment['gene_cut'])
        airr_vj.p3v_length = len(v_segment['palindrome_3_end'])

        # FIXME: WHY i NEED TO PUT IT FIRST?
        airr_vj.p5j_length = len(j_segment['palindrome_5_end'])

        # VJ
        airr_vj.n1_length = realization_dict['vj_ins'].value
        airr_vj.np1_length = airr_vj.p3v_length + airr_vj.n1_length + airr_vj.p5j_length
        airr_vj.np1 = vj_segment['gene_segment']  # This include the palindromic insertions

        # J

        airr_vj.j_data.germline_start = j_segment['gene_ini'] + 1
        airr_vj.j_data.germline_end = j_segment['gene_end'] + 1
        airr_vj.j_data.sequence_start = airr_vj.np1_length + (
                    airr_vj.v_data.sequence_end - airr_vj.v_data.sequence_start - 1)
        airr_vj.j_data.sequence_end = airr_vj.j_data.sequence_start + len(j_segment['gene_cut']) - 1

        return airr_vj.to_dict()




class IgorModel_Parms:
    """
    Class to get a list of Events directly from the *_parms.txt
    :param model_parms_file: Igor parms file path.
    """
    def __init__(self, model_parms_file=None):
        ## Parms file representation
        self.Event_list = list() # list of Rec_event
        self.Edges      = list()
        self.ErrorRate_dict  = dict()

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
        tmpstr = "{ 'len Event_list': " + str(len(self.Event_list)) \
                 +", 'len Egdes': " + str(len(self.Edges)) \
                 +", 'len ErrorRate': " + str(len(self.ErrorRate_dict)) + " }"
        return tmpstr
        #return "{ Event_list, Egdes, ErrorRate}"

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

    @classmethod
    def from_database(cls, db):
        print("Loading Model Parms from database.")

    @classmethod
    def make_default_VJ(cls, df_genomicVs, df_genomicJs, lims_deletions=None, lims_insertions=None):
        """Create a default VJ model from V and J genes dataframes
            lims_deletions tuple with min and maximum value for deletions, e.g. (-4,20)
            lims_insertions tuple with min and maximum value for deletions, e.g. (0,30)
        """
        cls = IgorModel_Parms()

        if lims_deletions is None:
            lims_deletions = (-4, 17)

        if lims_insertions is None:
            lims_insertions = (0, 41)

        # Add events to Event_list
        for event_nickname in Igor_VJ_default_nickname_list:
            event_dict = IgorRec_Event_default_dict[event_nickname]
            if event_nickname == 'j_choice':
                event_dict["priority"] = 6

            event = IgorRec_Event.from_dict(event_dict)
            cls.Event_list.append(event)
            if event.event_type == 'DinucMarkov':
                value_list = ['A', 'C', 'G', 'T']
                name_list = ['' for val in value_list]
                event_df = pd.DataFrame.from_dict({'name': name_list, 'value': value_list})
                event_df.index.name = 'id'
                cls.set_event_realizations_from_DataFrame(event_nickname, event_df)
            elif event.event_type == 'Deletion':
                value_list = list(range(*lims_deletions))
                name_list = ['' for val in value_list]
                event_df = pd.DataFrame.from_dict({'name': name_list, 'value': value_list})
                event_df.index.name = 'id'
                cls.set_event_realizations_from_DataFrame(event_nickname, event_df)
            elif event.event_type == 'Insertion':
                value_list = list(range(*lims_insertions))
                name_list = ['' for val in value_list]
                event_df = pd.DataFrame.from_dict({'name': name_list, 'value': value_list})
                event_df.index.name = 'id'
                cls.set_event_realizations_from_DataFrame(event_nickname, event_df)
            elif event.event_type == 'GeneChoice':
                if event_nickname == 'v_choice':
                    cls.set_event_realizations_from_DataFrame(event_nickname, df_genomicVs)
                elif event_nickname == 'j_choice':
                    cls.set_event_realizations_from_DataFrame(event_nickname, df_genomicJs)
            else:
                print("Unrecognized type of event. There are only 4 types of events:")
                print(" - GeneChoice")
                print(" - Deletions")
                print(" - Insertions")
                print(" - DinucMarkov")

        # Update names
        # for event in self.Event_list:
        #     event.name ="jodete"
        # Now edges
        cls.set_Edges_from_dict(Igor_VJ_default_parents_dict)

        # Error Rate
        cls.ErrorRate_dict = {'error_type': 'SingleErrorRate', 'error_values': '0.000396072'}

        return cls

    @classmethod
    def make_default_VDJ(cls, df_genomicVs, df_genomicDs, df_genomicJs, lims_deletions=None, lims_insertions=None):
        """Create a default VJ model from V and J genes dataframes
            lims_deletions tuple with min and maximum value for deletions, e.g. (-4,20)
            lims_insertions tuple with min and maximum value for deletions, e.g. (0,30)
        """
        cls = IgorModel_Parms()

        if lims_deletions is None:
            lims_deletions = (-4, 17)

        if lims_insertions is None:
            lims_insertions = (0, 41)

        for event_nickname in Igor_VDJ_default_nickname_list:
            event_dict = IgorRec_Event_default_dict[event_nickname]
            event = IgorRec_Event.from_dict(event_dict)
            cls.Event_list.append(event)

            if event.event_type == 'DinucMarkov':
                value_list = ['A', 'C', 'G', 'T']
                name_list = ['' for val in value_list]
                event_df = pd.DataFrame.from_dict({'name': name_list, 'value': value_list})
                event_df.index.name = 'id'
                cls.set_event_realizations_from_DataFrame(event_nickname, event_df)
            elif event.event_type == 'Deletion':
                value_list = list(range(*lims_deletions))
                name_list = ['' for val in value_list]
                event_df = pd.DataFrame.from_dict({'name': name_list, 'value': value_list})
                event_df.index.name = 'id'
                cls.set_event_realizations_from_DataFrame(event_nickname, event_df)
            elif event.event_type == 'Insertion':
                value_list = list(range(*lims_insertions))
                name_list = ['' for val in value_list]
                event_df = pd.DataFrame.from_dict({'name': name_list, 'value': value_list})
                event_df.index.name = 'id'
                cls.set_event_realizations_from_DataFrame(event_nickname, event_df)
            elif event.event_type == 'GeneChoice':
                if event_nickname == 'v_choice':
                    cls.set_event_realizations_from_DataFrame(event_nickname, df_genomicVs)
                elif event_nickname == 'd_gene':
                    cls.set_event_realizations_from_DataFrame(event_nickname, df_genomicDs)
                elif event_nickname == 'j_choice':
                    cls.set_event_realizations_from_DataFrame(event_nickname, df_genomicJs)
                else:
                    print("ERROR: GeneChoice event "+event.nickname+" is not a default nickname.")

            else:
                print("ERROR: Unrecognized type of event. There are only 4 types of events:")
                print(" - GeneChoice")
                print(" - Deletions")
                print(" - Insertions")
                print(" - DinucMarkov")

        cls.update_events_name()

        # Now edges
        cls.set_Edges_from_dict(Igor_VDJ_default_parents_dict)

        # Error Rate
        cls.ErrorRate_dict = {'error_type': 'SingleErrorRate', 'error_values': '0.000396072'}

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
                event_realization.id = index
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
                event_realization.id = index
                event_realization.value = ndels
            print(event_nickname, " limits : ", limits)

    def load_Insertion_realizations_by_nickname(self, event_nickname: str, limits=(0, 24)):
        event = self.get_Event(event_nickname)
        if event.event_type == 'Insertion':
            start, end = limits
            # FIXME: VALIDATE FOR POSITIVE VALUES
            for index, nins in enumerate(range(start, end)):
                event_realization = IgorEvent_realization()
                event_realization.id = index
                event_realization.value = nins
            print(event_nickname, " limits : ", limits)

    def load_DinucMarkov_realizations_by_nickname(self, event_nickname: str):
        event = self.get_Event(event_nickname)
        if event.event_type == 'DinucMarkov':
            for index, nt_char in enumerate(['A', 'C', 'G', 'T']):
                event_realization = IgorEvent_realization()
                event_realization.id = index
                event_realization.value = nt_char

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

        self.gen_EventDict_DataFrame()

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
                    realization.id = int(realizData[2])
                elif event.event_type == "DinucMarkov":
                    realization.value = realizData[0]
                    realization.id = int(realizData[1])
                else:
                    realization.value = int(realizData[0])
                    realization.id = int(realizData[1])

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

    def add_Egde(self, parent_nickname, child_nickname):
        try:
            parent_name = self.dictNicknameName[parent_nickname]
            child_name = self.dictNicknameName[child_nickname]
            # TODO: CHECK IF EDGE exist!
            self.Edges.append([parent_name, child_name])
            self.getBayesGraph()
            #self.Edges_dict[child_nickname].append(parent_nickname)
        except Exception as e:
            print("Edge : ", parent_nickname, child_nickname, " couldn't be added.")
            print(e)
            pass

    def set_Edges_from_dict(self, parents_dict):
        try:
            for child_nickname, parents in parents_dict.items():
                for parent_nickname in parents:
                    parent_name = self.dictNicknameName[parent_nickname]
                    child_name = self.dictNicknameName[child_nickname]
                    self.Edges.append([parent_name, child_name])
            self.getBayesGraph()
        except Exception as e:
            print("set_Edges_from_dict : ", parent_nickname, child_nickname, " couldn't be added.")
            print(e)
            pass

    def remove_Edge(self, parent_nickname, child_nickname):
        try:
            parent_name = self.dictNicknameName[parent_nickname]
            child_name = self.dictNicknameName[child_nickname]
            # TODO: CHECK IF EDGE exist!
            new_Edges = [edge for edge in self.Edges if not (parent_name == edge[0] and child_name == edge[1])]
            self.Edges = new_Edges
            self.getBayesGraph()
            #self.Edges_dict[child_nickname].append(parent_nickname)
        except Exception as e:
            print("Edge : ", parent_nickname, child_nickname, " couldn't be added.")
            print(e)
            pass

    def set_event_realizations_from_DataFrame(self, event_nickname, df):
        # FIXME: unnecesary copy find a better way.
        new_Event_list = list()
        for event in self.Event_list:
            if event.nickname == event_nickname:
                event.update_realizations_from_dataframe(df)
            new_Event_list.append(event)
        self.Event_list = new_Event_list
        self.gen_EventDict_DataFrame()

    def read_ErrorRate(self, ofile):
        lastPos  = ofile.tell()
        line = ofile.readline()
        strip_line = line.rstrip('\n')  # Remove end of line character
        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        while strip_line[0] == '#':
            # TODO: SAVE THE FOLLOWING TEXT AFTER # AS ERROR TYPE
            self.ErrorRate_dict = dict()
            self.ErrorRate_dict['error_type'] = strip_line[1:]
            lastPos = ofile.tell()
            line = ofile.readline()
            strip_line = line.rstrip('\n').rstrip()
            error = strip_line
            self.ErrorRate_dict['error_values'] = error

            # if 'SingleErrorRate' == strip_line[1:] :
            #     lastPos  = ofile.tell()
            #     line = ofile.readline()
            #     strip_line = line.rstrip('\n').rstrip()
            #     error = strip_line
            #     self.ErrorRate = {"SingleErrorRate" : error }
        ofile.seek(lastPos)

    # FIXME: FINISH THIS METHOD
    def write_model_parms(self, filename=None):
        """Writes a model graph structure from a model params object.
        Note that for now this method does not read the error rate information.
        """
        if filename is None:
            filename = "tmp_mdl_parms.txt"
        # Sort events in list with the your specific preference.
        # FIXME: FIND ANOTHER WAY TO WRITE IN A CORRECT ORDER
        # igor_nickname_list = ["v_choice", "j_choice", "d_gene", "v_3_del"]
        # self.get_Event(nicknameList)
        # self.Event_list
        strSepChar = ";"
        try:
            import os
            #print("AAAAAAAAAAAAAAA:", os.path.dirname(filename), filename)
            os.makedirs(os.path.dirname(filename), exist_ok=True)
        except Exception as e:
            print("WARNING: write_model_parms path ", e)

        print("Writing model parms in file ", filename)
        with open(filename, "w") as ofile:
            # 1. Write events
            self.write_Event_list(ofile, delimiter=strSepChar)

            # 2. Write Edges
            self.write_Edges(ofile, delimiter=strSepChar)

            # 3. Write ErrorRate
            self.write_ErrorRate(ofile, delimiter=strSepChar)

    def write_Event_list(self, ofile, delimiter=None):
        if delimiter is None:
            strSepChar = ';'
        else:
            strSepChar = delimiter

        ofile.write("@Event_list\n")
        # for event in self.Event_list:
        # for nickname in Igor_nickname_list: # Igor_nicknameList is in IgorDefaults.py
        for event in self.Event_list:
            try:
                ## self.write_event(ofile, event:IgorRec_Event)
                # event = self.get_Event(nickname)
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
                str_realization_list = ""
                for strLine in str_df.split("\n"):
                    # Asi se tiene = ['id', 'value', 'name']
                    strLine_list = strLine.split(strSepChar)
                    # print(strLine, strLine_list)
                    if len(strLine_list) > 1:
                        if event.event_type == "GeneChoice":
                            str_id = strLine_list[0]
                            str_value = strLine_list[1]
                            str_name = strLine_list[2]
                            # Asi se quiere = ['name', 'value', 'id']
                            str_realization = str_name + strSepChar + str_value + strSepChar + str_id
                        else:
                            str_id = strLine_list[0]
                            str_value = strLine_list[1]
                            # Asi se quiere = ['value', 'id']
                            str_realization = str_value + strSepChar + str_id

                        str_realization_list = str_realization_list + "%" + str_realization + "\n"
                ofile.write(str_realization_list)

            except Exception as e:
                print("ERROR: write_Event_list, ", event.nickname)
                print(e)
                pass

    def write_Edges(self, ofile, delimiter=None):
        if delimiter is None:
            strSepChar = ';'
        else:
            strSepChar = delimiter
        ofile.write("@Edges\n")
        try:
            for edge in self.Edges:
                ofile.write("%"+edge[0]+strSepChar+edge[1]+"\n")
        except Exception as e:
            print("ERROR: write_Edges")
            print(e)
            pass

    def write_ErrorRate(self, ofile, delimiter=None):
        if delimiter is None:
            strSepChar = ';'
        else:
            strSepChar = delimiter
        ofile.write("@ErrorRate\n")
        ofile.write("#"+self.ErrorRate_dict['error_type']+"\n")
        ofile.write(self.ErrorRate_dict['error_values']+"\n")

    def get_EventsNickname_list(self):
        return [event.nickname for event in self.Event_list]

    def get_EventsName_list(self):
        return [event.name for event in self.Event_list]

    def get_Event(self, event_nickname_or_name, by_nickname=True):
        """Returns the RecEvent with corresponding name or nickname."""
        if by_nickname:
            for ev in self.Event_list:
                if ev.nickname == event_nickname_or_name:
                    return ev
            raise Exception(
                'RecEvent with nickname \"' + event_nickname_or_name + "\" not found.")
        else:
            for ev in self.Event_list:
                if ev.name == event_nickname_or_name:
                    return ev
            raise Exception(
                'RecEvent with name \"' + event_nickname_or_name + "\" not found.")

    def gen_EventDict_DataFrame(self):
        self.Event_dict = dict()
        #dictio = dict()
        for event in self.Event_list:
            #dictio[event.nickname] = event.get_realization_DataFrame()
            self.Event_dict[event.nickname] = event.get_realization_DataFrame()

        self.gen_NameNickname_dict()
        self.getBayesGraph()

    def gen_NameNickname_dict(self):
        self.dictNameNickname = dict()
        for event in self.Event_list:
            event.update_name()
            self.dictNameNickname[event.name] = event.nickname
        # return dictio
        self.dictNicknameName = {v: k for k, v in self.dictNameNickname.items()}

    def get_event_dict(self, str_key, str_value):
        """
        Return a python dictionary of the event_dict, like ('nickname', 'priority')
        {'v_choice:7, 'd_gene':6, ...}
        """
        dicto = dict()
        for event in self.Event_list:
            event_dict = event.to_dict()
            dicto[event_dict[str_key]] = event_dict[str_value]
        return dicto

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

    def plot_Graph(self, ax=None, **kwargs): # FIXME: ALLOW the possibility to pass an ax like ax=None):
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

        try:
            import hvplot.networkx as hvnx
            print("hvplot")
            graph = hvnx.draw(self.G, with_labels=True, FontSize=10, pos=pos, alpha=0.5,
                            arrowstyle='fancy', arrowsize=2000, node_size=1000, width=400, height=400)
                            ##, arrows=True, arrowsize=20, node_size=800, font_size=10, font_weight='bold')
            return graph

        except ImportError as e:
            try:
                if ax is None:
                    #print("matplotlib")
                    import matplotlib.pyplot as plt
                    fig, ax = plt.subplots()
                ax.set_aspect('equal')
                nx.draw(self.G, pos=pos, ax=ax, with_labels=True, arrows=True, arrowsize=20,
                        node_size=800, font_size=10, font_weight='bold')  # FIXME: make a better plot: cutting edges.

                return ax
            except Exception as e:
                print(e)
                raise

    def get_Event_dependencies(self, strEvent):
        print(strEvent)
        return list(self.G.predecessors(strEvent))

    def get_Event_list_sorted(self):
        # FIXME: GENERALIZE THIS PROCESS with the parents priority
        # Order events by priority.
        events_list_with_parents = list()
        for event in self.Event_list:
            events_list_with_parents.append((event, list(self.G.predecessors(event.nickname))))
        # print(events_list_with_parents)
        sorted_events_list_with_parents = sorted(events_list_with_parents,
                                                 key=lambda tupla: (tupla[0].priority, -len(tupla[1])), reverse=True)
        return [sorted_event_parent[0] for sorted_event_parent in sorted_events_list_with_parents]

    def from_scenario(self, scenario, strEvent):
        return self.Event_dict[strEvent].loc[scenario[strEvent]]

    # TODO: GIVEN A SCENARIO A DICT WITH REALIZATIONS
    def realiz_dict_from_scenario(self, scenario): #IgorScenario):
        realizations_dict = dict()
        for nickname_key in scenario.realizations_ids_dict:
            if nickname_key == 'mismatches':
                realizations_dict[nickname_key] = scenario[nickname_key]
            elif nickname_key == 'mismatcheslen':
                realizations_dict[nickname_key] = scenario[nickname_key]
            else:
                event = self.get_Event(nickname_key)
                print(nickname_key, scenario[nickname_key], event.event_type)
                if event.event_type == 'DinucMarkov':
                    realizations_dict[nickname_key] = list()
                    for realiz_id in scenario[nickname_key]:
                        realizations_dict[nickname_key].append(event.realizations[realiz_id])
                else:
                    realizations_dict[nickname_key] =  event.realizations[ scenario[nickname_key] ]
            # except Exception as e:
            #     print("ERROR: ", nickname_key, " while parsing to realizations.")
            #     print(e)
        return realizations_dict



    def update_events_name(self):
        for event in self.Event_list:
            event.update_name()
        self.gen_NameNickname_dict()

    # FIXME: SCENARIO FROM CSV LINE IN GENERATED SEQUENCES
    def get_scenario_from_line_CSV(self, str_line, file_header_list, sep=';'):
        dicto = dict()
        str_line_list = str_line.split(sep)
        for str_header in file_header_list:
            if str_header == 'seq_index':
                pass
        return dicto

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

    def __getitem__(self, item):
        return self.realizations[item]

    def __str__(self):
        return str(self.to_dict())

    def __lt__(self, other):
        # TODO: Less parents bigger event
        return self.priority < other.priority
        #        if ( self.priority < other.priority ):
        #            return True
        #        elif ( self.priority == other.priority  ):
        #            # FIXME: less dependencies should be on top

    def __gt__(self, other):
        # TODO: Less parents bigger event
        return self.priority > other.priority

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

    def add_realization(self):
        realization = IgorEvent_realization()
        self.realizations.append(realization)
        self.realizations = sorted(self.realizations)
        self.update_name()

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
        cls.update_name()
        return cls

    def update_realizations_from_fasta(self, flnGenomic):
        from Bio import SeqIO
        if self.event_type == 'GeneChoice':
            for index, record in enumerate(SeqIO.parse(flnGenomic, "fasta")):
                event_realization = IgorEvent_realization()
                event_realization.id = index
                event_realization.value = record.seq
                event_realization.name = record.description
                self.add_realization(event_realization)

    def export_realizations_to_fasta(self, flnGenomic):
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        sequences_list = list()
        for realization in self.realizations:
            record = SeqRecord(Seq(realization.value), realization.name, '', '')
            sequences_list.append(record)

        SeqIO.write(sequences_list, flnGenomic, "fasta")

    def update_realizations_from_dataframe(self, dataframe):
        """
        Update realizations with a dataframe (index, value, name)
        """
        self.realizations = list()
        for index, row in dataframe.iterrows():
            dict_realiz = row.to_dict()
            # print(index, dict_realiz)
            dict_realiz['index'] = index
            realiz = IgorEvent_realization.from_dict(dict_realiz)
            self.realizations.append(realiz)
            self.realizations = sorted(self.realizations)
        self.update_name()

    @classmethod
    def from_default_nickname(cls, nickname:str):
        cls = IgorRec_Event.to_dict(IgorRec_Event_default_dict[nickname])
        return cls

    # TODO:
    def add_realization(self, realization):
        """Add a realization to the RecEvent realizations list."""
        self.realizations.append(realization)
        self.realizations = sorted(self.realizations)
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
        return pd.DataFrame.from_records([realiz.to_dict() for realiz in self.realizations], index='id').sort_index()

class IgorEvent_realization:
    """A small class storing for each RecEvent realization its name, value and
    corresponding index.
    """
    __slots__ = ('id', 'name', 'value')
    def __init__(self):
        self.id = "" #index
        self.name  = "" #name
        self.value = "" #value

    def __lt__(self, other):
        return self.id < other.id

    def __gt__(self, other):
        return self.id > other.id

    def __str__(self):
        if self.name == "":
            return "{value};{id}".format(value=self.value, id=self.id)
            # return str(self.value)+";"+str(self.id)
        else:
            return "{name};{value};{id}".format(name=self.name, value=self.value, id=self.id)
            # return self.name+";"+str(self.value)+";"+str(self.id)

    def __repr__(self):
        return "Event_realization(" + str(self.id) + ")"

    def to_dict(self):
        return {
            'id': self.id,
            'value': self.value,
            'name': self.name
            }

    @classmethod
    def from_tuple(cls, id, value, name=""):
        cls = IgorEvent_realization()
        cls.id = id
        cls.value = value
        cls.name = name
        return cls

    @classmethod
    def from_dict(self, event_dict:dict):
        cls = IgorEvent_realization()
        cls.id = event_dict['index']
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
        self.model_marginals_file = ""
        if model_marginals_file is not None:
            self.read_model_marginals(model_marginals_file)

    #  @d_3_del
    #  $Dim[3,21,21]
    #  #[d_gene,0],[d_5_del,0]
    #  %0,0,0,1.6468e-08,0.00482319,1.08101e-09,0.0195311,0.0210679,0.0359338,0.0328678,2.25686e-05,4.97463e-07,0,9.31048e-08,1.01642e-05,0.000536761,0.0260845,0.0391021,0.319224,0.289631,0.211165
    #  #[d_gene,0],[d_5_del,1]
    #  %0,0,6.86291e-08,2.00464e-09,0.00163832,2.02919e-06,0.0306066,0.0126832,0.000872623,0.016518,0.00495292,0.000776747,4.45576e-05,0.000667902,0.00274004,0.00435049,0.300943,0.182499,0.13817,0.302534,0

    @classmethod
    def make_uniform_from_parms(cls, parms:IgorModel_Parms):
        cls = IgorModel_Marginals()
        cls.initialize_uniform_from_model_parms(parms)
        return cls

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

    def initialize_uniform_event_from_model_parms(self, event_nickname, parms:IgorModel_Parms):
        event = parms.get_Event(event_nickname)
        if event.event_type == 'DinucMarkov':
            # do something
            dimension = len(parms.get_Event(event.nickname).realizations)
            narr = np.ones(dimension * dimension) / (dimension * dimension)
            self.marginals_dict[event.nickname] = narr
        else:
            dimensions = [len(parms.get_Event(strEvent).realizations) for strEvent in
                          self.network_dict[event.nickname]]
            # print(dimensions[-1])
            narr = np.ones(dimensions) / dimensions[-1]

            self.marginals_dict[event.nickname] = narr

    def initialize_uniform_from_model_parms(self, parms:IgorModel_Parms):
        self.network_dict = dict()
        for key, value in parms.Edges_dict.items():
            self.network_dict[key] = value + [key]

        # Create a marginal for
        self.marginals_dict = dict()
        for event in parms.Event_list:
            self.initialize_uniform_event_from_model_parms(event.nickname, parms)
            # if event.event_type == 'DinucMarkov':
            #     # do something
            #     dimension = len(parms.get_Event(event.nickname).realizations)
            #     narr = np.ones(dimension*dimension) / (dimension*dimension)
            #     self.marginals_dict[event.nickname] = narr
            # else:
            #     dimensions = [len(parms.get_Event(strEvent).realizations) for strEvent in
            #                   self.network_dict[event.nickname]]
            #     # print(dimensions[-1])
            #     narr = np.ones(dimensions) / dimensions[-1]
            #
            #     self.marginals_dict[event.nickname] = narr

    def write_model_marginals(self, filename=None, model_parms=None):
        # self.marginals_dict = {}
        # self.network_dict = {}
        if filename is None:
            filename = "tmp_mdl_marginals.txt"

        if model_parms is None:
            print("model parms need it")
            raise
        parms = model_parms #IgorModel_Parms(model_parms_file=model_parms_file)

        if filename is None:
            filename = "tmp_mdl_marginals.txt"

        try:
            import os
            os.makedirs(os.path.dirname(filename), exist_ok=True)
        except Exception as e:
            print("WARNING: IgorModel_Marginals.write_model_marginals path ", e)

        print("Writing model marginals in file ", filename)
        with open(filename, "w") as fw:
            for event in parms.Event_list:
                strEvent = event.nickname
                # strEvent = "v_choice"
                self.write_event_probabilities(fw, strEvent)

    def write_event_probabilities(self, ofile, event_nickname):
        import itertools
        np_array = self.marginals_dict[event_nickname]
        parents_list = self.network_dict[event_nickname]
        parents_to_write = parents_list[:-1]
        # dims_list = tuple( map( lambda x: list(range(x)), array_to_write.shape) )
        dims_list = tuple(map(lambda x: list(range(x)), np_array.shape[:-1]))

        ofile.write("@" + event_nickname + "\n")
        str_shape = str(list(np_array.shape)).replace(" ", "").replace("[", "").replace("]", "")
        strDim = "$Dim[" + str_shape + "]\n"
        ofile.write(strDim)
        for elem in itertools.product(*dims_list):
            title = ""
            #     print(parents_to_write)
            title = str(list(zip(parents_to_write, elem)))
            title = title.replace("[", "#")
            title = title.replace("]", "\n")
            title = title.replace(" ", "")
            title = title.replace("(", "[")
            title = title.replace(")", "]")
            title = title.replace("\'", "")
            ofile.write(title)
            slice_index = tuple(list(elem) + [None])
            linea = str(list(np_array[slice_index].flat))
            linea = linea.replace(" ", "")
            linea = linea.replace(" ", "")
            linea = linea.replace("[", "%")
            linea = linea.replace("]", "\n")
            ofile.write(linea)

        # ofile.write()


class IgorAnchors:
    def __init__(self, flnVanchors, flnJanchors):
        self.flnVanchors = flnVanchors
        self.flnJanchors = flnJanchors
        self.df_Vanchors = pd.read_csv(flnVanchors, sep=';')
        self.df_Janchors = pd.read_csv(flnJanchors, sep=';')
        # rename indices.

class IgorScenario:
    def __init__(self):
        self.seq_index = -1
        self.scenario_rank = -1
        self.scenario_proba_cond_seq = -1
        #self.events_ordered_list = list()
        self.realizations_ids_dict = dict()
        # given a templated list with ids
        self.mdl = None
        # self.mdl.parms.Event_dict[strEv].loc[self.id_d_gene]['name']

    def __getitem__(self, key):
        return self.realizations_ids_dict[key]

    def to_dict(self):
        dictScenario = dict()
        dictScenario['seq_index'] = self.seq_index
        dictScenario['scenario_rank'] = self.scenario_rank
        dictScenario['scenario_proba_cond_seq'] = self.scenario_proba_cond_seq
        dictScenario.update(self.realizations_ids_dict)
        return dictScenario

    # TODO: This method should return a scenario in a fasta format with corresponding ID and events
    def get_scenario_fasta(self, mdl:IgorModel):
        str_fasta = ""
        # sort events to construct fasta sequence:
        mdl.parms.Event_list
        for key in mdl.xdata.keys():
            self.realizations_ids_dict[key]

        return str_fasta

    def set_model(self, mdl:IgorModel):
        """ Initiate scenario dictionary with a IgorModel """
        for key in mdl.xdata.keys():
            self.realizations_ids_dict[key] = -1

    # TODO: in DEV - FINISH THIS METHOD
    def set_model_from_headers(self, header_line:str):
        # seq_index;scenario_rank;scenario_proba_cond_seq;GeneChoice_V_gene_Undefined_side_prio7_size35;GeneChoice_J_gene_Undefined_side_prio7_size14;GeneChoice_D_gene_Undefined_side_prio6_size2;Deletion_V_gene_Three_prime_prio5_size21;Deletion_D_gene_Five_prime_prio5_size21;Deletion_D_gene_Three_prime_prio5_size21;Deletion_J_gene_Five_prime_prio5_size23;Insertion_VD_genes_Undefined_side_prio4_size31;DinucMarkov_VD_genes_Undefined_side_prio3_size16;Insertion_DJ_gene_Undefined_side_prio2_size31;DinucMarkov_DJ_gene_Undefined_side_prio1_size16;Mismatches
        header_line = "seq_index;scenario_rank;scenario_proba_cond_seq;GeneChoice_V_gene_Undefined_side_prio7_size35;GeneChoice_J_gene_Undefined_side_prio7_size14;GeneChoice_D_gene_Undefined_side_prio6_size2;Deletion_V_gene_Three_prime_prio5_size21;Deletion_D_gene_Five_prime_prio5_size21;Deletion_D_gene_Three_prime_prio5_size21;Deletion_J_gene_Five_prime_prio5_size23;Insertion_VD_genes_Undefined_side_prio4_size31;DinucMarkov_VD_genes_Undefined_side_prio3_size16;Insertion_DJ_gene_Undefined_side_prio2_size31;DinucMarkov_DJ_gene_Undefined_side_prio1_size16;Mismatches"
        header_fields = header_line.split(";")
        events_list = header_fields[3:]
        print("hoajs")

    # FIXME:
    @classmethod
    def load_FromLineBestScenario(cls, line, delimiter=";"):
        #seq_index;scenario_rank;scenario_proba_cond_seq;GeneChoice_V_gene_Undefined_side_prio7_size35;GeneChoice_J_gene_Undefined_side_prio7_size14;GeneChoice_D_gene_Undefined_side_prio6_size2;Deletion_V_gene_Three_prime_prio5_size21;Deletion_D_gene_Five_prime_prio5_size21;Deletion_D_gene_Three_prime_prio5_size21;Deletion_J_gene_Five_prime_prio5_size23;Insertion_VD_genes_Undefined_side_prio4_size31;DinucMarkov_VD_genes_Undefined_side_prio3_size16;Insertion_DJ_gene_Undefined_side_prio2_size31;DinucMarkov_DJ_gene_Undefined_side_prio1_size16;Mismatches
        cls = IgorScenario()
        linesplit = line.split(delimiter)
        for ii in range(len(linesplit)):
            # TODO: find a better way to do this, if is a list keep it as list
            if (ii in [ 11, 13, 14 ]):
                linesplit[ii] = linesplit[ii]
            else:
                linesplit[ii] = linesplit[ii].replace("(", "").replace(")", "")

    @classmethod
    def load_FromSQLRecord(cls, sqlRecordScenario:list, sql_scenario_name_type_list:list):
        cls = IgorScenario()
        for ii, (col_name, tipo) in enumerate(sql_scenario_name_type_list):
            if col_name == 'seq_index':
                cls.seq_index = int(sqlRecordScenario[ii])
            elif col_name == 'scenario_rank':
                cls.scenario_rank = int(sqlRecordScenario[ii])
            elif col_name == 'scenario_proba_cond_seq':
                cls.scenario_proba_cond_seq = float(sqlRecordScenario[ii])
            else:
                if tipo == 'integer':
                    cls.realizations_ids_dict[col_name] = int(sqlRecordScenario[ii])
                else:
                    cls.realizations_ids_dict[col_name] = eval(sqlRecordScenario[ii])

        return cls

    def export_to_AIRR_line(self, scenario_col_list:list, sep='\t'):
        str_line = ""
        self.seq_index = -1
        self.scenario_rank = -1
        self.scenario_proba_cond_seq = -1

        # n_d_5_del = self.mdlParms.Event_dict[strEv].loc[self.id_d_5_del]['value']
        # name_D = self.mdlParms.Event_dict[strEv].loc[self.id_d_gene]['name']
        # header_list=['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score', 'd_score', 'j_score'])
        # sequence_id	sequence	rev_comp	productive	v_call	d_call	j_call	c_call	sequence_alignment	germline_alignment	junction	junction_aa	v_score	v_cigar	d_score	d_cigar	j_score	j_cigar	c_score	c_cigar	vj_in_frame	stop_codon	v_identity	v_evalue	d_identity	d_evalue	j_identity	j_evalue	v_sequence_start	v_sequence_end	v_germline_start	v_germline_end	d_sequence_start	d_sequence_end	d_germline_start	d_germline_end	j_sequence_start	j_sequence_end	j_germline_start	j_germline_end	junction_length	np1_length	np2_length	duplicate_count	consensus_count
        airr_header_list = ["sequence_id", "sequence", "rev_comp", "productive", "v_call", "d_call", "j_call", "c_call",
                            "sequence_alignment", "germline_alignment", "junction", "junction_aa", "v_score", "v_cigar", "d_score", "d_cigar", "j_score", "j_cigar", "c_score", "c_cigar", "vj_in_frame", "stop_codon", "v_identity", "v_evalue", "d_identity", "d_evalue", "j_identity", "j_evalue", "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end", "d_sequence_start", "d_sequence_end", "d_germline_start", "d_germline_end", "j_sequence_start", "j_sequence_end", "j_germline_start", "j_germline_end", "junction_length", "np1_length", "np2_length", "duplicate_count", "consensus_count"]

        from pygor3 import IgorModel_Parms
        mdl_parms = IgorModel_Parms()
        # mdl_parms = self.mdl.parms
        # TODO: No general way, just select between VJ OR VDJ, SO RECHECK IN MODEL IF 'd_gene' is present and make arrangement.
        airr_line_list = list()
        for event_nickname in scenario_col_list:
            event_realization_id = self.realizations_ids_dict[event_nickname]
            event_realization_value = mdl_parms.Event_dict[event_nickname].loc[event_realization_id]['value']
            event_realization_name = mdl_parms.Event_dict[event_nickname].loc[event_realization_id]['name']
            airr_line_list.append(str(self.seq_index))

            bs_realiz = mdl_parms.realiz_dict_from_scenario(bs)
            mdl_parms.from_scenario()
            # GeneChoice
            self.realizations_ids_dict[event_nickname]
            # Deletions

            # Insertions

            # DinucMarkov

        str_line = sep.join([self.seq_index, self.scenario_rank, self.scenario_proba_cond_seq])


        return str_line

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

