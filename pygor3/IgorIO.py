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
_flag_verbose = False
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
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

from pathlib import Path
from typing import Union
import tempfile

### GENERIC FUNCTIONS
def genLabel(strName):
    """
    Generation of label for a simple identification of genomic template sequence.
    """
    aaa = strName.split("|")
    if len(aaa) > 1:
        return aaa[1]
    else:
        return strName


v_genLabel = np.vectorize(genLabel)


def command_from_dict_options(dicto: dict):
    """ Return igor options from dictionary"""
    dicto_copy = copy.deepcopy(dicto)
    cmd = ''
    for key in dicto_copy.keys():
        if dicto_copy[key]['active']:
            if dicto_copy[key]['active'] is None:
                cmd = cmd + " " + key + " "
            else:
                if dicto_copy[key]['value'] is None:
                    cmd = cmd + " " + key + " "
                else:
                    cmd = cmd + " " + key + " " + str( dicto_copy[key]['value'] )

            if 'dict_options' in dicto_copy[key]:
                if dicto_copy[key]['dict_options'] is not None:
                    # print(key, dicto[key]['dict_options'])
                    cmd = cmd + " " + command_from_dict_options(dicto_copy[key]['dict_options'])
    return cmd


def run_command(cmd):
    """from http://blog.kagesenshi.org/2008/02/teeing-python-subprocesspopen-output.html
    """
    # print(cmd)
    # p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # stdout = []
    # while True:
    #     line = p.stdout.readline()
    #     line = line.decode("utf-8")
    #     stdout.append(line)
    #     # print (line, end='')
    #     if line == '' and p.poll() != None:
    #         break
    # return ''.join(stdout)
    try:
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return p.stdout
    except TypeError as e:
        try:
            p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            return p.stdout
        except Exception as e:
            raise e
    except Exception as e:
        raise e


def execute_command_generator(cmd):
    popen = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    # popen = subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def run_command_print(cmd):
    try:
        std_output_str = ""
        for path in execute_command_generator(cmd):
            print(path, end="")
            std_output_str = std_output_str + '\n'

        return std_output_str
    except Exception as e:
        raise e


def run_command_no_output(cmd):
    """from http://blog.kagesenshi.org/2008/02/teeing-python-subprocesspopen-output.html
    """
    # p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    try:
        from subprocess import PIPE
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        return p
    except TypeError as e:
        p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        return p
    except Exception as e:
        raise e
    # p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    # return p


# FIXME: IT IS BETTER TO USE DECORATORS FOR VARIABLES LIKE igor_batchname and update the dependencies on that automatically?


### IGOR INPUT SEQUENCES  ####

class IgorIndexedSequence:
    """
    Return a IgorIndexedSequence instance
    """

    def __init__(self, seq_index=-1, sequence=''):
        self.seq_index = seq_index
        self.sequence = sequence

    def __str__(self):
        return str(self.to_dict())

    def to_dict(self):
        """
        Return a IgorIndexedSequence instance as a python dictionary.
        """
        dictIndexedSequence = {
            "seq_index": self.seq_index, \
            "sequence": self.sequence
        }

        return dictIndexedSequence

    @classmethod
    def load(cls, seq_index, sequence):
        cls = IgorIndexedSequence()
        try:
            cls.seq_index = seq_index
            cls.sequence = sequence
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
            cls.seq_index = int(csvsplit[0])
            cls.sequence = int(csvsplit[1])
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

class IgorAnchors:
    def __init__(self, path_ref_genome: Union[None, str, Path] = None,
                 flnVanchors: Union[None, str, Path] = None,
                 flnJanchors: Union[None, str, Path] = None,
                 df_Vanchors: Union[None,pd.DataFrame] = None,
                 df_Janchors: Union[None, pd.DataFrame] = None,
                 sep=';'):
        self.flnVanchors = flnVanchors
        self.flnJanchors = flnJanchors
        self.path_ref_genome = path_ref_genome
        self.update_default_filenames()
        self.df_Vanchors = df_Vanchors
        self.df_Janchors = df_Janchors
        try:
            self.load_dataframes()
        except:
            pass

    def load_dataframes(self):

        try:
            self.load_V_dataframe()
        except Exception as e:
            print("V anchors were not loaded.")
            pass

        try:
            self.load_J_dataframe()
        except Exception as e:
            print("J anchors were not loaded.")
            pass

    @classmethod
    def load_from_path(cls, path_ref_genome):
        cls = IgorAnchors()
        cls.path_ref_genome = path_ref_genome
        cls.update_default_filenames()
        try:
            cls.load_V_dataframe()
            cls.load_J_dataframe()
        except Exception as e:
            raise e
        else:
            return cls

    def update_default_filenames(self, path_ref_genome: Union[None, str, Path] = None):
        if path_ref_genome is not None:
            self.path_ref_genome = path_ref_genome

        if self.path_ref_genome is not None:
            self.flnVanchors = str(self.path_ref_genome) + "/" + "V_gene_CDR3_anchors.csv"
            self.flnJanchors = str(self.path_ref_genome) + "/" + "J_gene_CDR3_anchors.csv"

    def load_V_dataframe(self, flnVanchors: Union[None, str, Path] = None, sep=';'):
        try:
            if flnVanchors is not None:
                self.flnVanchors = flnVanchors
            self.df_Vanchors = pd.read_csv(self.flnVanchors, sep=sep).set_index('gene')
        except Exception as e:
            raise e

    def load_J_dataframe(self, flnJanchors: Union[None, str, Path] = None, sep=';'):
        try:
            if flnJanchors is not None:
                self.flnJanchors = flnJanchors
            self.df_Janchors = pd.read_csv(self.flnJanchors, sep=sep).set_index('gene')
        except Exception as e:
            raise e

        # rename indices.



class IgorRefGenome:
    def __init__(self, fln_genomicVs: Union[None, str, Path] = None,
                 fln_genomicDs: Union[None, str, Path] = None,
                 fln_genomicJs: Union[None, str, Path] = None,
                 fln_V_gene_CDR3_anchors: Union[None, str, Path] = None,
                 fln_J_gene_CDR3_anchors: Union[None, str, Path] = None,
                 path_ref_genome: Union[None, str, Path]=None):
        """Class to save genomic information"""
        # FIXME: find a better way to add a default value for this and also the "/" separator
        self.path_ref_genome = None
        self.fln_genomicVs = None
        self.fln_genomicDs = None
        self.fln_genomicJs = None
        self.fln_V_gene_CDR3_anchors = None
        self.fln_J_gene_CDR3_anchors = None

        if path_ref_genome is None:
            self.path_ref_genome = "."

        if fln_genomicVs is not None:
            self.fln_genomicVs = fln_genomicVs  # "genomicVs.fasta"
        if fln_genomicDs is not None:
            self.fln_genomicDs = fln_genomicDs  # "genomicDs.fasta"
        if fln_genomicJs is not None:
            self.fln_genomicJs = fln_genomicJs  # "genomicJs.fasta"

        if fln_V_gene_CDR3_anchors is not None:
            self.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors  # "V_gene_CDR3_anchors.csv"
        if fln_J_gene_CDR3_anchors is not None:
            self.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors  # "J_gene_CDR3_anchors.csv"

        self.df_genomicVs = None # ['id', 'name', 'value'] 'id' as index of dataframe
        self.df_genomicDs = None # ['id', 'name', 'value'] 'id' as index of dataframe
        self.df_genomicJs = None # ['id', 'name', 'value'] 'id' as index of dataframe

        self.df_V_anchors = None # ['gene', 'anchor_index'] 'gene' as index of dataframe
        self.df_J_anchors = None # ['gene', 'anchor_index'] 'gene' as index of dataframe

        # self.dict_genomicVs = None  # (self.df_genomicVs.set_index('name').to_dict())['value']
        # self.dict_genomicDs = None
        # self.dict_genomicJs = None

        # self.df_V_ref_genome = None
        # self.df_J_ref_genome = None

        try:
            if (self.fln_genomicVs is not None) or (self.fln_genomicJs is not None):
                self.load_dataframes_from_ref_genome_files()
        except:
            pass

    @property
    def dict_genomicVs(self):
        try:
            df_genomicVs_copy = self.df_genomicVs.copy()
            return (df_genomicVs_copy.set_index('name').to_dict())['value']
        except Exception as e:
            return None

    @dict_genomicVs.setter
    def dict_genomicVs(self, new_dict:dict):
        # TODO: IN DEV FINISH ME TO ADD NEW GENOMIC REFERENCES
        self.df_genomicVs

    @property
    def dict_genomicDs(self):
        try:
            df_genomicDs_copy = self.df_genomicDs.copy()
            return (df_genomicDs_copy.set_index('name').to_dict())['value']
        except Exception as e:
            return None

    @property
    def dict_genomicJs(self):
        try:
            df_genomicJs_copy = self.df_genomicJs.copy()
            return (df_genomicJs_copy.set_index('name').to_dict())['value']
        except Exception as e:
            return None

    @property
    def df_V_ref_genome(self):
        """
        Property that returns V dataframe of genomic templates with anchor's column if present.
        """
        try:
            if (self.df_V_anchors is None) or (self.df_V_anchors.empty):
                return self.df_genomicVs
            else:
                return get_join_genomics_anchors_dataframes(self.df_genomicVs, self.df_V_anchors)
            # return self.df_genomicVs.set_index('name').join(self.df_V_anchors.set_index('gene')).reset_index()
        except Exception as e:
            return None

    @property
    def df_J_ref_genome(self):
        """
        Property that returns J dataframe of genomic templates with anchor's column if present.
        """
        try:
            if (self.df_J_anchors is None) or (self.df_J_anchors.empty):
                return self.df_genomicJs
            else:
                return get_join_genomics_anchors_dataframes(self.df_genomicJs, self.df_J_anchors)
            # return self.df_genomicJs.set_index('name').join(self.df_J_anchors.set_index('gene')).reset_index()
        except Exception as e:
            return None

    @df_V_ref_genome.setter
    def df_V_ref_genome(self, new_V_ref_genome: pd.DataFrame):
        self.df_genomicVs = new_V_ref_genome[['name', 'value']].copy()
        self.df_V_anchors = new_V_ref_genome.drop(columns=['value']).set_index('name')


    @df_J_ref_genome.setter
    def df_J_ref_genome(self, new_J_ref_genome: pd.DataFrame):
        self.df_genomicJs = new_J_ref_genome[['name', 'value']].copy()
        self.df_J_anchors = new_J_ref_genome.drop(columns=['value']).set_index('name')

    @property
    def V(self):
        try:
            return self.df_V_ref_genome
        except Exception as e:
            return None

    @property
    def D(self):
        try:
            return self.df_genomicDs
        except Exception as e:
            return None

    @property
    def J(self):
        try:
            return self.df_J_ref_genome
        except Exception as e:
            return None

    def __repr__(self):
        return self.to_dict()

    def to_dict(self):
        dicto = dict()
        if self.df_V_ref_genome is not None:
            dicto['V'] = self.df_V_ref_genome

        if self.df_genomicDs is not None:
            dicto['D'] = self.df_genomicDs

        if self.df_J_ref_genome is not None:
            dicto['J'] = self.df_J_ref_genome

        return dicto

    def __getitem__(self, key):
        if key == 'V':
            return self.df_V_ref_genome
        elif key == 'D':
            return self.df_genomicDs
        elif key == 'J':
            return self.df_J_ref_genome
        else:
            return None


    @classmethod
    def load_default(cls, IgorSpecie, IgorChain, modelpath=None, ref_genome=None):
        """
        Return IgorRefGenome
        """

        ref_genome_fln_dict = get_default_fln_dict_ref_genomes_species_chain(IgorSpecie, IgorChain,
                                                                         modelspath=modelpath, ref_genome_path=ref_genome)
        print(ref_genome_fln_dict)
        try:
            cls = IgorRefGenome(**ref_genome_fln_dict)
            cls.specie = IgorSpecie
            cls.chain = IgorChain
        except Exception as e:
            raise e
        return cls

    @classmethod
    def load_FromSQLRecord_list(cls, sqlrecords_genomicVs=None, sqlrecords_genomicDs=None, sqlrecords_genomicJs=None,
                                sqlrecords_V_gene_CDR3_anchors=None, sqlrecords_J_gene_CDR3_anchors=None):
        """
        Return IgorRefGenome from database records.
        """
        cls = IgorRefGenome()
        # TODO: make query to database

        cls.df_genomicVs = pd.DataFrame.from_records(sqlrecords_genomicVs, columns=['id', 'name', 'value']).set_index(
            'id')

        # Fasta to dataframe
        try:
            # df_V_anchors = pd.read_csv(self.fln_V_gene_CDR3_anchors, sep=';')
            df_V_anchors = pd.DataFrame.from_records(sqlrecords_V_gene_CDR3_anchors,
                                                     columns=['id', 'gene', 'anchor_index']).set_index(('id'))

            cls.df_V_ref_genome = cls.df_genomicVs.set_index('name').join(df_V_anchors.set_index('gene')).reset_index()
            cls.dict_genomicVs = (cls.df_genomicVs.set_index('name').to_dict())['value']
        except Exception as e:
            print('No V genes were found.')
            print(e)
            pass

        # J genes
        cls.df_genomicJs = pd.DataFrame.from_records(sqlrecords_genomicJs, columns=['id', 'name', 'value']).set_index(
            'id')
        try:
            df_J_anchors = pd.DataFrame.from_records(sqlrecords_J_gene_CDR3_anchors,
                                                     columns=['id', 'gene', 'anchor_index']).set_index(('id'))
            cls.df_J_ref_genome = cls.df_genomicJs.set_index('name').join(df_J_anchors.set_index('gene')).reset_index()
            # cls.dict_genomicJs = (cls.df_genomicJs.set_index('name').to_dict())['value']
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
    def load_from_path(cls, path_ref_genome:Union[str, Path]):
        """
        Return IgorRefGenome from directory path with default names:
        genomicVs.fasta, genomicDs.fasta, genomicJs.fasta,
        V_gene_CDR3_anchors.csv and J_gene_CDR3_anchors.csv
        :param path_ref_genome: Path of directory
        :return : IgorRefGenome
        """
        try:
            cls = IgorRefGenome()
            cls.path_ref_genome = path_ref_genome
            cls.update_fln_names(path_ref_genome=cls.path_ref_genome)
            cls.load_dataframes_from_ref_genome_files()
            return cls
        except Exception as e:
            raise e

    @classmethod
    def load_from_dataframe_genomics_dict(cls, df_genomics_dict:dict):
        """
        Return IgorRefGenome from directory path with default names:
        genomicVs.fasta, genomicDs.fasta, genomicJs.fasta,
        V_gene_CDR3_anchors.csv and J_gene_CDR3_anchors.csv
        :param df_genomics_dict: dictionary with 'V', 'J' and/or 'D' keys with pandas dataframes.
        :return : IgorRefGenome
        """
        cls = IgorRefGenome()
        if 'V' in df_genomics_dict:
            # Check columns name convention in df_genomics_dict['V']
            df_genome = get_dataframe_with_ref_genome_column_names(df_genomics_dict['V'])
            # cls.df_genomicVs = df_genomics_dict['V']
            cls.df_V_ref_genome = df_genome.copy() #df_genomics_dict['V']

        if 'D' in df_genomics_dict:
            df_genome = get_dataframe_with_ref_genome_column_names(df_genomics_dict['D'])
            cls.df_genomicDs = df_genome.copy() #df_genomics_dict['D']

        if 'J' in df_genomics_dict:
            df_genome = get_dataframe_with_ref_genome_column_names(df_genomics_dict['J'])
            cls.df_J_ref_genome = df_genome.copy()  # df_genomics_dict['D']
            # cls.df_genomicJs = df_genomics_dict['J']
            # cls.df_J_ref_genome = df_genomics_dict['J']

        return cls


    @staticmethod
    def get_imgt_list_species():
        from .imgt import get_species_list
        return get_species_list()

    @classmethod
    def load_VJ_from_IMGT_website(cls, imgt_species, imgt_chain, **kwargs):
        """
        Return IgorRefGenome from IMGT website:
        :param imgt_species: species in IMGT format
        :param imgt_chain: chain in IMGT format
        :param modelspath: (Optional) If specified will not be deleted.
        """
        try:
            from .imgt import download_ref_genome_VJ

            flag_temporal_dir = False
            if not ('modelspath' in kwargs):
                import tempfile
                with tempfile.TemporaryDirectory() as tmp_path:
                    kwargs['modelspath'] = tmp_path
                    flag_temporal_dir = True
                    df_genes_dict = download_ref_genome_VJ(imgt_species, imgt_chain, **kwargs)
                    ref_genome_path = kwargs['modelspath'] + "/" + imgt_species + "/" + imgt_chain + "/" + "ref_genome"
                    cls = IgorRefGenome.load_from_path(ref_genome_path)
            else:
                if kwargs['modelspath'] is None:
                    import tempfile
                    with tempfile.TemporaryDirectory() as tmp_path:
                        kwargs['modelspath'] = tmp_path
                        flag_temporal_dir = True
                        df_genes_dict = download_ref_genome_VJ(imgt_species, imgt_chain, **kwargs)
                        ref_genome_path = kwargs['modelspath'] + "/" + imgt_species + "/" + imgt_chain + "/" + "ref_genome"
                        cls = IgorRefGenome.load_from_path(ref_genome_path)
                else:
                    df_genes_dict = download_ref_genome_VJ(imgt_species, imgt_chain, **kwargs)
                    cls = IgorRefGenome.load_from_path(kwargs['modelspath'])

        except Exception as e:
            raise e
        else:
            return cls

    @classmethod
    def load_VDJ_from_IMGT_website(cls, imgt_species, imgt_chain, **kwargs):
        """
        Return IgorRefGenome from IMGT website:
        :param imgt_species: species in IMGT format
        :param imgt_chain: chain in IMGT format
        :param modelspath: (Optional) If specified will not be deleted.
        """
        try:
            from .imgt import download_ref_genome_VDJ

            flag_temporal_dir = False
            if not ('modelspath' in kwargs):
                import tempfile
                with tempfile.TemporaryDirectory() as tmp_path:
                    kwargs['modelspath'] = tmp_path
                    flag_temporal_dir = True
                    df_genes_dict = download_ref_genome_VDJ(imgt_species, imgt_chain, **kwargs)
                    ref_genome_path = kwargs['modelspath'] + "/" + imgt_species + "/"+ imgt_chain + "/" + "ref_genome"
                    cls = IgorRefGenome.load_from_path(ref_genome_path)
            else:
                if kwargs['modelspath'] is None:
                    import tempfile
                    with tempfile.TemporaryDirectory() as tmp_path:
                        kwargs['modelspath'] = tmp_path
                        flag_temporal_dir = True
                        df_genes_dict = download_ref_genome_VDJ(imgt_species, imgt_chain, **kwargs)
                        ref_genome_path = kwargs['modelspath'] + "/" + imgt_species + "/" + imgt_chain + "/" + "ref_genome"
                        cls = IgorRefGenome.load_from_path(ref_genome_path)
                else:
                    df_genes_dict = download_ref_genome_VDJ(imgt_species, imgt_chain, **kwargs)
                    cls = IgorRefGenome.load_from_path(kwargs['modelspath'])

        except Exception as e:
            raise e
        else:
            return cls

    def update_fln_names(self, path_ref_genome: Union[None, str] = None,
                         fln_genomicVs: Union[None, str] = None,
                         fln_genomicDs: Union[None, str] = None,
                         fln_genomicJs: Union[None, str] = None,
                         fln_V_gene_CDR3_anchors: Union[None, str] = None,
                         fln_J_gene_CDR3_anchors: Union[None, str] = None):
        """Update genomic filenames
            :param fln_genomicVs: Path of fasta file for V genomic templates,
            :param fln_genomicDs: Path of fasta file for D genomic templates,
            :param fln_genomicJs: Path of fasta file for J genomic templates,
            :param fln_V_gene_CDR3_anchors: Path of csv anchor file for V genes,
            :param fln_J_gene_CDR3_anchors: Path of csv anchor file for J genes
        """
        try:
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
        except Exception as e:
            e_message = "IgorRefGenome.update_fln_names : path_ref_genome " + str(self.path_ref_genome)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    # TODO: LOAD INSTANCE FROM DEFINED FILES, what is the difference btwn load_dataframes?
    def load_dataframes_from_dict(self, df_genomics_dict):
        # FIXME: IN DEV
        self.fln_genomicVs = None
        self.fln_genomicDs = None
        self.fln_genomicJs = None
        self.fln_V_gene_CDR3_anchors = None
        self.fln_J_gene_CDR3_anchors = None

        self.load_dataframes_from_ref_genome_files()

    def load_dataframes_from_ref_genome_files(self,
                                              fln_genomicVs: Union[None, str, Path] = None,
                                              fln_genomicDs: Union[None, str, Path] = None,
                                              fln_genomicJs: Union[None, str, Path] = None,
                                              fln_V_gene_CDR3_anchors: Union[None, str, Path] = None,
                                              fln_J_gene_CDR3_anchors: Union[None, str, Path] = None,
                                              sep=';'
                                              ):

        if fln_genomicVs is not None:
            self.fln_genomicVs = fln_genomicVs  # "genomicVs.fasta"
        if fln_genomicDs is not None:
            self.fln_genomicDs = fln_genomicDs  # "genomicDs.fasta"
        if fln_genomicJs is not None:
            self.fln_genomicJs = fln_genomicJs  # "genomicJs.fasta"

        if fln_V_gene_CDR3_anchors is not None:
            self.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors  # "V_gene_CDR3_anchors.csv"
        if fln_J_gene_CDR3_anchors is not None:
            self.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors  # "J_gene_CDR3_anchors.csv"

        # Fasta to dataframe
        # V genes
        try:
            self.load_genomicVs_from_file(self.fln_genomicVs)
            try:
                self.load_V_anchors_from_file(self.fln_V_gene_CDR3_anchors)
            except Exception as e:
                print(e)
                pass
        except Exception as e:
            print(e)
            pass

        # J genes
        try:
            self.load_genomicJs_from_file(self.fln_genomicJs)
            try:
                self.load_J_anchors_from_file(self.fln_J_gene_CDR3_anchors)
            except Exception as e:
                print(e)
                pass
        except Exception as e:
            # print("WARNING: separator", sep, " do not work for ", self.fln_genomicJs)
            print(e)
            pass

        # D genes
        try:
            self.load_genomicDs_from_file(self.fln_genomicDs)
            # self.df_genomicDs = get_dataframe_from_fasta(self.fln_genomicDs)
            # self.dict_genomicDs = (self.df_genomicDs.set_index('name').to_dict())['value']
        except Exception as e:
            self.fln_genomicDs = None
            print(e)
            pass

        # # Add anchors
        # try:
        #     df_V_anchors = get_anchors_dataframe_from_csv(self.fln_V_gene_CDR3_anchors, sep=sep)
        #     self.df_V_ref_genome = self.df_genomicVs.set_index('name').join(df_V_anchors).reset_index()
        #     self.df_V_ref_genome.index.name = 'id'
        #     self.dict_genomicVs = (self.df_genomicVs.set_index('name').to_dict())['value']
        # except Exception as e:
        #     print('No V genes were found.')
        #     print(e)
        #     pass
        #
        # # J genes
        # self.df_genomicJs = get_dataframe_from_fasta(self.fln_genomicJs)
        # try:
        #     df_J_anchors = pd.read_csv(self.fln_J_gene_CDR3_anchors, sep=';')
        #     self.df_J_ref_genome = self.df_genomicJs.set_index('name').join(
        #         df_J_anchors.set_index('gene')).reset_index()
        #     self.df_J_ref_genome.index.name = 'id'
        #     self.dict_genomicJs = (self.df_genomicJs.set_index('name').to_dict())['value']
        # except Exception as e:
        #     print('No J genes were found.')
        #     print(e)
        #     pass


        # return df_V_ref_genome, df_J_ref_genome

    def get_anchors_dict(self):
        dict_anchor_index = dict()
        dict_anchor_index['V'] = self.df_V_ref_genome.set_index('name')['anchor_index'].to_dict()
        dict_anchor_index['J'] = self.df_J_ref_genome.set_index('name')['anchor_index'].to_dict()
        return dict_anchor_index

    def write_ref_genome(self,
                         fln_genomicVs: Union[None, str, Path] = None,
                         fln_genomicDs: Union[None, str, Path] = None,
                         fln_genomicJs: Union[None, str, Path] = None,
                         fln_V_gene_CDR3_anchors: Union[None, str, Path] = None,
                         fln_J_gene_CDR3_anchors: Union[None, str, Path] = None, sep=';'):
        """Save genomes in files
        :param fln_genomicVs: Output V gene fasta genomic file.
        :param fln_genomicDs: Output V gene fasta genomic file.
        :param fln_genomicJs: Output V gene fasta genomic file.
        :param fln_V_gene_CDR3_anchors: Output csv anchor file for V gene.
        :param fln_J_gene_CDR3_anchors: Output csv anchor file for J gene.
        """

        if fln_genomicVs is None:
            fln_genomicVs = self.fln_genomicVs

        if self.D is not None:
            if fln_genomicDs is None:
                fln_genomicDs = self.fln_genomicDs
        if fln_genomicJs is None:
            fln_genomicJs = self.fln_genomicJs

        if fln_V_gene_CDR3_anchors is None:
            fln_V_gene_CDR3_anchors = self.fln_V_gene_CDR3_anchors
        if fln_J_gene_CDR3_anchors is None:
            fln_J_gene_CDR3_anchors = self.fln_J_gene_CDR3_anchors

        try:
            write_ref_genome_files_from_dataframe(self.df_V_ref_genome, fln_genomicVs,
                                                  fln_V_gene_CDR3_anchors, sep=sep)

            if self.df_genomicDs is not None:
                write_ref_genome_files_from_dataframe(self.df_genomicDs, fln_genomicDs, sep=sep)

            write_ref_genome_files_from_dataframe(self.df_J_ref_genome, fln_genomicJs,
                                                  fln_J_gene_CDR3_anchors, sep=sep)

            # write_genetemplate_dataframe_to_fasta(self.fln_genomicVs, self.df_genomicVs)
            # write_genetemplate_dataframe_to_fasta(self.fln_genomicJs, self.df_genomicJs)
        except Exception as e:
            raise e

        # try:
        #     write_genetemplate_dataframe_to_fasta(self.fln_genomicDs, self.df_genomicDs)
        # except Exception as e:
        #     pass
        #
        # try:
        #     write_geneanchors_dataframe_to_csv(self.fln_V_gene_CDR3_anchors, self.df_V_ref_genome)
        #     write_geneanchors_dataframe_to_csv(self.fln_J_gene_CDR3_anchors, self.df_J_ref_genome)
        # except Exception as e:
        #     raise e


    def write_ref_genome_dir(self, ref_genome_dir_path, sep=';'):
        """Write ref_genome directory in path
        :param ref_genome_dir_path: Path to directory to save ref_genomic files.
        :param sep: default = ';' to save anchors files.
        """
        fln_dict = get_default_ref_genome_fln_paths(ref_genome_path=ref_genome_dir_path)
        # TODO: CHECK FOR D GENES
        if (self.D is None) and 'fln_genomicDs' in fln_dict.keys():
            fln_dict['fln_genomicDs'] = None
        # print(fln_dict)
        if self.df_genomicVs is None:
            del fln_dict['fln_genomicVs']

        if self.df_genomicDs is None:  # ['id', 'name', 'value'] 'id' as index of dataframe
            del fln_dict['fln_genomicDs']

        if self.df_genomicJs is None:  # ['id', 'name', 'value'] 'id' as index of dataframe
            del fln_dict['fln_genomicJs']

        if (self.df_V_anchors is None) or (self.df_V_anchors.empty):
            del fln_dict['fln_V_gene_CDR3_anchors']

        if (self.df_J_anchors is None) or (self.df_J_anchors.empty):
            del fln_dict['fln_J_gene_CDR3_anchors']

        self.write_ref_genome(sep=sep, **fln_dict)

    def clean_empty_anchors(self):
        """
        Remove genes without anchors
        """
        try:
            tmp_df = self.df_V_ref_genome[self.df_V_ref_genome['anchor_index'].notna()].copy()
            tmp_df.reset_index(inplace=True)
            tmp_df['anchor_index'] = tmp_df['anchor_index'].apply(lambda x: int(x))
            tmp_df = tmp_df.drop(columns=['id'])
            tmp_df.index.name = 'id'
            self.df_V_ref_genome = tmp_df.copy()

            tmp_df = self.df_J_ref_genome[self.df_J_ref_genome['anchor_index'].notna()].copy()
            tmp_df.reset_index(inplace=True)
            tmp_df['anchor_index'] = tmp_df['anchor_index'].apply(lambda x: int(x))
            tmp_df = tmp_df.drop(columns=['id'])
            tmp_df.index.name = 'id'
            self.df_J_ref_genome = tmp_df.copy()

        except Exception as e:
            raise e

    def load_genomicVs_from_file(self, fln_genomicVs):
        """
        Load V genes dataframe (df_genomicVs) to IgorRefGenome object.
        :param fln_genomicVs: Filename of fasta gene templates for V gene
        """
        try:

            self.df_genomicVs = get_dataframe_from_fasta(fln_genomicVs)
            # self.dict_genomicVs = (self.df_genomicVs.set_index('name').to_dict())['value']
            self.fln_genomicVs = fln_genomicVs
            if _flag_verbose:
                print("Loaded genomic V templates from file ", self.fln_genomicVs)
        except Exception as e:
            e_message = "load_genomicVs_from_file " + str(fln_genomicVs)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def load_genomicDs_from_file(self, fln_genomicDs):
        """
        Load D genes dataframe (df_genomicDs) to IgorRefGenome object.
        :param fln_genomicDs: Filename of fasta gene templates for D gene
        """
        try:
            self.df_genomicDs = get_dataframe_from_fasta(fln_genomicDs)
            # self.dict_genomicDs = (self.df_genomicDs.set_index('name').to_dict())['value']
            self.fln_genomicDs = fln_genomicDs
            if _flag_verbose:
                print("Loaded genomic D templates from file ", self.fln_genomicDs)
        except Exception as e:
            e_message = "load_genomicDs_from_file " + str(fln_genomicDs)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def load_genomicJs_from_file(self, fln_genomicJs):
        """
        Load J genes dataframe (df_genomicJs) to IgorRefGenome object.
        :param fln_genomicJs: Filename of fasta gene templates for J gene
        """
        try:
            self.df_genomicJs = get_dataframe_from_fasta(fln_genomicJs)
            # self.dict_genomicJs = (self.df_genomicJs.set_index('name').to_dict())['value']
            self.fln_genomicJs = fln_genomicJs
            if _flag_verbose:
                print("Loaded genomic J templates from file ", self.fln_genomicJs)
        except Exception as e:
            e_message = "load_genomicJs_from_file " + str(fln_genomicJs)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def load_V_anchors_from_file(self, fln_V_gene_CDR3_anchors, sep=';'):
        """
        Load CDR3 V anchors dataframe (df_V_anchors) to IgorRefGenome object.
        :param fln_V_gene_CDR3_anchors: Filename of csv anchors file templates for V gene
        """
        try:
            self.df_V_anchors = get_anchors_dataframe_from_csv(fln_V_gene_CDR3_anchors, sep=sep)
            # _df_V_anchors = get_anchors_dataframe_from_csv(fln_V_gene_CDR3_anchors, sep=sep)
            # self.df_V_ref_genome = get_join_genomics_anchors_dataframes(self.df_V_ref_genome, _df_V_anchors)
            self.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors
            if _flag_verbose:
                print("Loaded genomic V CDR3 anchors from file ", self.fln_V_gene_CDR3_anchors)
        except Exception as e:
            e_message = "load_V_anchors_from_file " + str(fln_V_gene_CDR3_anchors)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def load_J_anchors_from_file(self, fln_J_gene_CDR3_anchors, sep=';'):
        """
        Load V genes dataframe (df_genomicVs) to IgorRefGenome object.
        :param fln_J_gene_CDR3_anchors: Filename of fasta gene templates for V gene
        """
        try:
            self.df_J_anchors = get_anchors_dataframe_from_csv(fln_J_gene_CDR3_anchors, sep=sep)
            self.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors
            if _flag_verbose:
                print("Loaded genomic J CDR3 anchors from file ", self.fln_J_gene_CDR3_anchors)
        except Exception as e:
            e_message = "load_J_anchors_from_file " + str(fln_J_gene_CDR3_anchors)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])


class IgorAlignment_data:
    def __init__(self):
        self.seq_index = -1
        self.gene_id = -1
        self.score = -1
        self.offset = 0
        self.insertions = list()
        self.deletions = list()
        self.mismatches = list()
        self.length = 0
        self.offset_5_p = 0
        self.offset_3_p = 0

        self.strGene_name = ""
        self.strGene_class = ""
        self.strGene_seq = ""

        self.anchor_in_read = None

    def __str__(self):
        return str(self.to_dict())

    def to_dict(self):
        dictAlignment_data = {
            "seq_index": self.seq_index, \
            "gene_id": self.gene_id, \
            "score": self.score, \
            "offset": self.offset, \
            "insertions": self.insertions, \
            "deletions": self.deletions, \
            "mismatches": self.mismatches, \
            "length": self.length, \
            "offset_5_p": self.offset_5_p, \
            "offset_3_p": self.offset_3_p, \
            "strGene_name": self.strGene_name, \
            "strGene_class": self.strGene_class, \
            "strGene_seq": self.strGene_seq
        }

        return dictAlignment_data

    @classmethod
    def load_FromCSVLine(cls, csvline, strGene_name="", delimiter=";"):
        # seq_index;gene_name;score;offset;insertions;deletions;mismatches;length;5_p_align_offset;3_p_align_offset
        cls = IgorAlignment_data()
        csvsplit = csvline.replace("\n", "").split(";")
        try:
            cls.seq_index = int(csvsplit[0])
            cls.strGene_name = str(csvsplit[1])
            cls.score = float(csvsplit[2])
            cls.offset = int(csvsplit[3])
            cls.insertions = eval(csvsplit[4].replace("{", "[").replace("}", "]"))
            cls.deletions = eval(csvsplit[5].replace("{", "[").replace("}", "]"))
            cls.mismatches = eval(csvsplit[6].replace("{", "[").replace("}", "]"))
            cls.length = int(csvsplit[7])
            cls.offset_5_p = int(csvsplit[8])
            cls.offset_3_p = int(csvsplit[9])
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
            cls.seq_index = int(sqlRecordAlign[0])
            cls.gene_id = int(sqlRecordAlign[1])
            cls.score = float(sqlRecordAlign[2])
            cls.offset = int(sqlRecordAlign[3])
            cls.insertions = eval(sqlRecordAlign[4])
            cls.deletions = eval(sqlRecordAlign[5])
            cls.mismatches = eval(sqlRecordAlign[6])
            cls.length = int(sqlRecordAlign[7])
            cls.offset_5_p = int(sqlRecordAlign[8])
            cls.offset_3_p = int(sqlRecordAlign[9])
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
class IgorEvent_realization:
    """A small class storing for each RecEvent realization its name, value and
    corresponding index.
    """
    __slots__ = ('id', 'name', 'value')

    def __init__(self):
        self.id = ""  # index
        self.name = ""  # name
        self.value = ""  # value

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
        return str( self.to_dict() )
        # return "Event_realization(" + str(self.id) + ")"

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
    def from_pandas(cls, df:pd.DataFrame):
        return df.to_records()

    @classmethod
    def from_dict(self, event_dict: dict):
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

        self._pd_realizations = None
        try:
            self.update_pd_realizations_from_realizations()
        except:
            pass

    @property
    def pd_realizations(self):
        return self._pd_realizations

    @pd_realizations.setter
    def pd_realizations(self, value):
        self.update_realizations_from_dataframe(value)
        self._pd_realizations = self.get_realization_DataFrame()

    @pd_realizations.deleter
    def pd_realizations(self):
        del self._pd_realizations

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
        self.update_pd_realizations_from_realizations()


    def get_realization(self, index:Union[int,list]) -> IgorEvent_realization:
        """get realization object by index
        :param index: Id of realization
        :return : IgorEvent_realization
        """
        try:
            # TODO: CHANGE THIS TO USE PANDAS DATAFRAME
            tmp_list = [realization for realization in self.realizations if realization.id == index]
            return tmp_list[0]
        except Exception as e:
            raise e

    @classmethod
    def from_dict(cls, dict_IgorRec_Event: dict):
        """Returns a IgorRec_Event based on dictionary
        """
        cls = IgorRec_Event(dict_IgorRec_Event["event_type"], dict_IgorRec_Event["seq_type"],
                            dict_IgorRec_Event["seq_side"], dict_IgorRec_Event["priority"],
                            dict_IgorRec_Event["nickname"])
        # 'event_type', 'seq_type', 'seq_side', 'priority', and 'nickname'
        # FIXME: Is better to make this class as an extension of a dictionary container?
        # Given the nickname complete the events with
        # cls.nickname = dict_IgorRec_Event["nickname"]
        # cls.event_type = dict_IgorRec_Event["event_type"]
        # cls.seq_type = dict_IgorRec_Event["seq_type"]
        # cls.seq_side = dict_IgorRec_Event["seq_side"]
        # cls.priority = dict_IgorRec_Event["priority"]
        cls.realizations = dict_IgorRec_Event["realizations"]  # TODO: CREATE FUNCTION TO GENERATE realizations vector
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

    def update_pd_realizations_from_realizations(self, realizations:Union[None,list]=None):
        if realizations is not None:
            self.realizations = realizations
        self.pd_realizations = self.get_realization_DataFrame()

    @classmethod
    def from_default_nickname(cls, nickname: str):
        cls = IgorRec_Event.to_dict(IgorRec_Event_default_dict[nickname])
        return cls

    # TODO:
    def add_realization(self, realization):
        """Add a realization to the RecEvent realizations list."""
        self.realizations.append(realization)
        self.realizations = sorted(self.realizations)
        self.update_name()
        self.pd_realizations = self.get_realization_DataFrame()

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

    def set_realization_vector_GeneChoice(self, flnGenomic: str):
        """
        Sets a realization vector from a filename
        :param flnGenomic: fasta file with the genomic template IMGT or other template.
        """
        # FIXME: FINISH IT
        # TODO: Add realizations from fasta file.
        from Bio import SeqIO
        # for record in list(SeqIO.parse(flnGenomic, "fasta")):

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
        """ Return an Event realizations as a pandas DataFrame with id, value and name columns
        and attributes
        - event_type
        - seq_type
        - seq_side
        - priority
        - NOT event name, because conflicts with pandas dataframe name
        - nickname
        """
        try:
            if len(self.realizations) == 0:
                df_event = pd.DataFrame(columns=['id', 'value', 'name']).set_index('id')
            else:
                df_event = pd.DataFrame.from_records([realiz.to_dict() for realiz in self.realizations], index='id').sort_index()
            df_event.event_type = self.event_type
            df_event.seq_type = self.seq_type
            df_event.seq_side = self.seq_side
            df_event.priority = self.priority
            df_event.nickname = self.nickname
            return df_event
        # except KeyError:
        #
        #     df_event.event_type = self.event_type
        #     df_event.seq_type = self.seq_type
        #     df_event.seq_side = self.seq_side
        #     df_event.priority = self.priority
        #     df_event.nickname = self.nickname
        #     return df_event
        except Exception as e:
            raise e



class IgorModel_Parms:
    """
    Class to get a list of Events directly from the *_parms.txt
    :param model_parms_file: Igor parms file path.
    """

    def __init__(self, model_parms_file=None,
                 fln_V_gene_CDR3_anchors=None,
                 fln_J_gene_CDR3_anchors=None):
        ## Parms file representation
        # IgorRec_Event
        self.Event_list = list()  # list of IgorRec_Event
        self.Edges = list()
        self.ErrorRate_dict = dict()

        # Now there are properties
        # self.df_V_ref_genome = None
        # self.df_D_ref_genome = None
        # self.df_J_ref_genome = None

        self.df_V_anchors = None
        self.df_J_anchors = None

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
            if not (fln_V_gene_CDR3_anchors is None):
                try:
                    self.attach_V_anchors_from_file(fln_V_gene_CDR3_anchors)
                except:
                    pass
            if not (fln_J_gene_CDR3_anchors is None):
                try:
                    self.attach_J_anchors_from_file(fln_J_gene_CDR3_anchors)
                except:
                    pass

            # self.get_EventDict_DataFrame()

    def __getitem__(self, item):
        return self.Event_dict[item]

    @property
    def event_GeneChoice_list(self):
        return [event for event in self.Event_list if event.event_type == 'GeneChoice']

    @property
    def event_Insertion_list(self):
        return [event for event in self.Event_list if event.event_type == 'Insertion']

    @property
    def event_Deletion_list(self):
        return [event for event in self.Event_list if event.event_type == 'Deletion']

    @property
    def event_DinucMarkov_list(self):
        return [event for event in self.Event_list if event.event_type == 'DinucMarkov']

    @property
    def event_GeneChoice_V(self) -> Union[IgorRec_Event, None]:
        """
        Return IgorRec_Event GeneChoice and V_gene event from self.Event_list (usual nickname 'v_choice')
        """
        try:
            GeneChoice_list = self.event_GeneChoice_list # [event for event in self.Event_list if event.event_type == 'GeneChoice']
            event_GeneChoice_V = [event for event in GeneChoice_list if event.seq_type == 'V_gene'][0]
        except IndexError:
            print("No V genes event found!")
            return None
        except Exception as e:
            raise e
        else:
            return event_GeneChoice_V

    @property
    def event_GeneChoice_D(self) -> Union[IgorRec_Event, None]:
        """
        Return IgorRec_Event GeneChoice and D_gene event from self.Event_list (usual nickname 'd_gene')
        """
        try:
            GeneChoice_list = self.event_GeneChoice_list #[event for event in self.Event_list if event.event_type == 'GeneChoice']
            event_GeneChoice_D = [event for event in GeneChoice_list if event.seq_type == 'D_gene'][0]
        except IndexError:
            # print("No D genes event found!")
            return None
        except Exception as e:
            raise e
        else:
            return event_GeneChoice_D

    @property
    def event_GeneChoice_J(self) -> Union[IgorRec_Event, None]:
        """
        Return IgorRec_Event GeneChoice and D_gene event from self.Event_list (usual nickname 'd_gene')
        """
        try:
            GeneChoice_list = self.event_GeneChoice_list #[event for event in self.Event_list if event.event_type == 'GeneChoice']
            event_GeneChoice_J = [event for event in GeneChoice_list if event.seq_type == 'J_gene'][0]
        except IndexError:
            print("No J genes event found!")
            return None
        except Exception as e:
            raise e
        else:
            return event_GeneChoice_J

    @property
    def df_V_ref_genome(self) -> Union[pd.DataFrame, None]:
        """
        Return pandas dataframe with anexed anchors if anchors are available.
        """
        try:
            # Use the function to join dataframes
            if self.df_V_anchors is None:
                return self.event_GeneChoice_V.get_realization_DataFrame()
            else:
                if isinstance(self.df_V_anchors, pd.DataFrame) :
                    if self.df_V_anchors.empty:
                        return self.event_GeneChoice_V.get_realization_DataFrame()
                    else:
                        return get_join_genomics_anchors_dataframes(
                            self.event_GeneChoice_V.get_realization_DataFrame(),
                            self.df_V_anchors)
        except Exception as e:
            raise e

    @property
    def df_D_ref_genome(self) -> Union[pd.DataFrame, None]:
        """
        Return pandas dataframe with anexed anchors if anchors are available.
        """
        try:
            if self.event_GeneChoice_D is None:
                return None
            else:
                df_all = self.event_GeneChoice_D.get_realization_DataFrame()
                columnas = df_all.columns.to_list()
                ini_cols = ['name', 'value']
                other_cols = list()
                for col in columnas:
                    if not col in ini_cols:
                        other_cols.append(col)
                new_order = ini_cols + other_cols
                return df_all[new_order]

        except Exception as e:
            raise e

    @property
    def df_J_ref_genome(self) -> Union[pd.DataFrame, None]:
        """
        Return pandas dataframe with anexed anchors if anchors are available.
        """
        try:
            # Use the function to join dataframes
            if self.df_J_anchors is None:
                return self.event_GeneChoice_J.get_realization_DataFrame()
            else:
                if isinstance(self.df_J_anchors, pd.DataFrame) :
                    if self.df_J_anchors.empty:
                        return self.event_GeneChoice_J.get_realization_DataFrame()
                    else:
                        return get_join_genomics_anchors_dataframes(
                            self.event_GeneChoice_J.get_realization_DataFrame(), self.df_J_anchors)
        except Exception as e:
            raise e


    def __str__(self):
        tmpstr = "{ 'len Event_list': " + str(len(self.Event_list)) \
                 + ", 'len Egdes': " + str(len(self.Edges)) \
                 + ", 'len ErrorRate': " + str(len(self.ErrorRate_dict)) + " }"
        return tmpstr
        # return "{ Event_list, Egdes, ErrorRate}"

    @classmethod
    def from_network_dict(cls, network_dict: dict):
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
        # mdl0.parms.Event_list[0].realizations[0])

        # load events from default dictionary.
        #
        return cls

    @classmethod
    def from_database(cls, db):
        print("Loading Model Parms from database.")

    @classmethod
    def make_default_VJ(cls, df_V_ref_genome, df_J_ref_genome, lims_deletions=None, lims_insertions=None):
        """Create a default VJ model from V and J genes dataframes
        :param df_V_ref_genome: Pandas Dataframe of Genome reference for V gene with CDR3 anchors
        :param df_J_ref_genome: Pandas Dataframe of Genome reference for J gene with CDR3 anchors
        :param lims_deletions: Tuple with min and maximum value for deletions, e.g. (-4,20). Negative numbers are palidromic insertions
        :param lims_insertions: Tuple with min and maximum value for deletions, e.g. (0,30)
        """
        cls = IgorModel_Parms()
        # df_genomicVs, df_genomicJs
        try:
            genomic_cols = ['name', 'value']
            df_genomicVs = df_V_ref_genome[genomic_cols]
            df_genomicJs = df_J_ref_genome[genomic_cols]
        except KeyError as e:
            print("ERROR: gene name column name should be 'name' and sequence column name should be 'value'")
            raise e
        except Exception as e:
            raise e

        if lims_deletions is None:
            lims_deletions = (-4, 17)

        if lims_insertions is None:
            lims_insertions = (0, 41)

        # Add events to Event_list
        for event_nickname in Igor_VJ_default_nickname_list:
            event_dict = IgorRec_Event_default_dict[event_nickname].copy()
            if event_nickname == 'j_choice':
                event_dict["priority"] = 6

            event = IgorRec_Event.from_dict(event_dict)
            cls.Event_list.append(event)


        cls.gen_NameNickname_dict()
        for edge_parent_child in Igor_VJ_default_Edges_parent_child_tuples:
            cls.add_Edge(*edge_parent_child)

        for event in cls.Event_list:
            event_nickname = event.nickname

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

        cls.Event_list = cls.get_Event_list_sorted()
        cls.update_events_name()

        cls.set_Edges_from_dict(Igor_VJ_default_parents_dict)

        # Error Rate
        cls.ErrorRate_dict = {'error_type': 'SingleErrorRate', 'error_values': '0.000396072'}

        try:
            # Attach anchors
            cls.df_V_anchors = get_df_anchors_from_df_ref_genome(cls.df_V_ref_genome)
            cls.df_J_anchors = get_df_anchors_from_df_ref_genome(cls.df_J_ref_genome)

        except Exception as e:
            print(e)
            pass
        return cls

        return cls

    @classmethod
    def make_default_VDJ(cls, df_V_ref_genome, df_D_ref_genome, df_J_ref_genome, lims_deletions=None, lims_insertions=None):
        """Create a default VJ model from V and J genes dataframes
        :param df_V_ref_genome: Pandas Dataframe of Genome reference for V gene with CDR3 anchors
        :param df_D_ref_genome: Pandas Dataframe of Genome reference for D gene
        :param df_J_ref_genome: Pandas Dataframe of Genome reference for J gene with CDR3 anchors
        :param lims_deletions: Tuple with min and maximum value for deletions, e.g. (-4,20). Negative numbers are palidromic insertions
        :param lims_insertions: Tuple with min and maximum value for deletions, e.g. (0,30)
        """
        cls = IgorModel_Parms()
        # df_genomicVs, df_genomicDs, df_genomicJs
        # FIXME: CREATE EVENTS AND EDGES

        try:
            genomic_cols = ['name', 'value']
            df_genomicVs = df_V_ref_genome[genomic_cols].copy()
            df_genomicDs = df_D_ref_genome[genomic_cols].copy()
            df_genomicJs = df_J_ref_genome[genomic_cols].copy()
        except KeyError as e:
            print("ERROR: gene name column name should be 'name' and sequence column name should be 'value'")
            raise e
        except Exception as e:
            raise e


        if lims_deletions is None:
            lims_deletions = (-4, 17)

        if lims_insertions is None:
            lims_insertions = (0, 41)

        for event_nickname in Igor_VDJ_default_nickname_list:
            event_dict = IgorRec_Event_default_dict[event_nickname].copy()
            event = IgorRec_Event.from_dict(event_dict)
            cls.Event_list.append(event)

        cls.gen_NameNickname_dict()
        # TODO: ADD EDGES # parent, child
        for edge_parent_child in Igor_VDJ_default_Edges_parent_child_tuples:
            cls.add_Edge(*edge_parent_child)

        for event in cls.Event_list:
            event_nickname = event.nickname

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
                if event.seq_type == 'V_gene':
                #if event_nickname == 'v_choice':
                    cls.set_event_realizations_from_DataFrame(event_nickname, df_genomicVs)
                elif event.seq_type == 'D_gene':
                #elif event_nickname == 'd_gene':
                    cls.set_event_realizations_from_DataFrame(event_nickname, df_genomicDs)
                elif event.seq_type == 'J_gene':
                # elif event_nickname == 'j_choice':
                    cls.set_event_realizations_from_DataFrame(event_nickname, df_genomicJs)
                else:
                    print("ERROR: GeneChoice event " + event.nickname + " is not a default nickname.")

            else:
                print("ERROR: Unrecognized type of event. There are only 4 types of events:")
                print(" - GeneChoice")
                print(" - Deletions")
                print(" - Insertions")
                print(" - DinucMarkov")

        cls.Event_list = cls.get_Event_list_sorted()
        cls.update_events_name()

        # Now edges
        cls.set_Edges_from_dict(Igor_VDJ_default_parents_dict)

        # Error Rate
        cls.ErrorRate_dict = {'error_type': 'SingleErrorRate', 'error_values': '0.000396072'}

        try:
            # Attach anchors
            cls.df_V_anchors = get_df_anchors_from_df_ref_genome(df_V_ref_genome)
            cls.df_J_anchors = get_df_anchors_from_df_ref_genome(df_J_ref_genome)

        except Exception as e:
            print(e)
            pass
        return cls

    @classmethod
    def make_default_VDJ_from_IgorRefGenome(cls, ref_genome:IgorRefGenome, lims_deletions=None, lims_insertions=None):
        """
        Return IgorModel_Parms from IgorRefGenome
        """
        cls = IgorModel_Parms.make_default_VDJ(ref_genome.df_genomicVs, ref_genome.df_genomicDs, ref_genome.df_genomicJs,
                                               lims_deletions=lims_deletions, lims_insertions=lims_insertions)
        cls.attach_anchors_from_files()
        return cls

    @classmethod
    def load_default(cls, IgorSpecie, IgorChain, modelpath=None, ref_genome_path=None):  # rcParams['paths.igor_models']):
        """
        Return IGoR default model parms for species and chain specified.
        """
        flnModelParms, flnModelMargs = get_default_models_paths_species_chain(IgorSpecie, IgorChain, modelpath=modelpath)

        fln_dict = get_default_fln_dict_ref_genomes_species_chain(IgorSpecie, IgorChain, modelspath=modelpath, ref_genome_path=ref_genome_path)

        cls = IgorModel_Parms(model_parms_file=flnModelParms,
                              fln_V_gene_CDR3_anchors=fln_dict['fln_V_gene_CDR3_anchors'],
                              fln_J_gene_CDR3_anchors=fln_dict['fln_J_gene_CDR3_anchors'])

        return cls

    def load_events_from_dict(self, dicto):
        print(dicto)

    # TODO: Check how the imgt functions return data
    def load_GeneChoice_realizations_by_nickname(self, event_nickname: str, flnGenomic):
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
        try:
            print("Reading Parms filename from: ", filename)
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
                if strip_line == "@ErrorRate":
                    self.read_ErrorRate(ofile)
            self.model_parms_file = filename

            self.gen_EventDict_DataFrame()
        except Exception as e:
            raise e

    # save in Event_list
    def read_Event_list(self, ofile):
        lastPos = ofile.tell()
        line = ofile.readline()
        strip_line = line.rstrip('\n')  # Remove end of line character
        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        # event = Rec_Event()

        while strip_line[0] == '#':
            # get the metadata of the event list
            event_metadata = strip_line[1:].split(";")  # GeneChoice;V_gene;Undefined_side;7;v_choice
            event_metadata[3] = int(event_metadata[3])  # change priority to integer
            event = IgorRec_Event(*event_metadata)
            # self.G.add_node(event.nickname)
            # Now read the realizations (or possibilities)
            lastPos = ofile.tell()
            line = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)

            while strip_line[0] == '%':
                realization = IgorEvent_realization()
                realizData = strip_line[1:].split(";")
                if event.event_type == "GeneChoice":
                    realization.name = realizData[0]
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
                lastPos = ofile.tell()
                line = ofile.readline()
                strip_line = line.rstrip('\n')  # Remove end of line character
                strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)

            self.Event_list.append(event)
        ofile.seek(lastPos)

    def read_Edges(self, ofile):
        # print "read_Edges"
        lastPos = ofile.tell()
        line = ofile.readline()
        strip_line = line.rstrip('\n')  # Remove end of line character
        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        while strip_line[0] == '%':
            edge = strip_line[1:].split(';')
            self.Edges.append(edge)
            # read nextline
            lastPos = ofile.tell()
            line = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        ofile.seek(lastPos)

    def add_Edge(self, parent_nickname, child_nickname):
        try:
            parent_name = self.dictNicknameName[parent_nickname]
            child_name = self.dictNicknameName[child_nickname]
            # TODO: CHECK IF EDGE exist!
            if not ([parent_name, child_name] in self.Edges):
                self.Edges.append([parent_name, child_name])
            self.getBayesGraph()
            # self.Edges_dict[child_nickname].append(parent_nickname)
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
                    if not ([parent_name, child_name] in self.Edges):
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
            # self.Edges_dict[child_nickname].append(parent_nickname)
        except Exception as e:
            print("Edge : ", parent_nickname, child_nickname, " couldn't be added.")
            print(e)
            pass

    def set_event_realizations_from_DataFrame(self, event_nickname, df):
        """
        Set realizations of a defined event from a pandas dataframe.
        :param event_nickname: Event nickname to set the realizations
        :param df: Pandas dataframe with 'id', 'value', 'name' columns (id as index)
        """
        # FIXME: unnecesary copy find a better way.
        #  if GeneChoice if anchors present attach it to
        #  self.df_V_anchors or
        #  self.df_J_anchors

        original_event_name = self.dictNicknameName[event_nickname]
        new_Event_list = list()
        for event in self.Event_list:
            if event.nickname == event_nickname:
                event.update_realizations_from_dataframe(df)
            new_Event_list.append(event)
        self.Event_list = new_Event_list
        # FIXME: UPDATE EDGES NAMES
        self.gen_NameNickname_dict()
        new_event_name = self.dictNicknameName[event_nickname]
        for edge in self.Edges:
            if edge[0] == original_event_name:
                edge[0] = new_event_name
            if edge[1] == original_event_name:
                edge[1] = new_event_name

        event = self.get_Event(event_nickname)
        self.Event_dict[event_nickname] = event.get_realization_DataFrame()
        self.gen_NameNickname_dict()
        self.getBayesGraph()
        # self.gen_EventDict_DataFrame()

        # IF Anchors in dataframe attach it.
        if event.nickname == self.event_GeneChoice_V:
            try:
                # Attach anchors
                self.df_V_anchors = get_df_anchors_from_df_ref_genome(df)
            except Exception as e:
                print(e)
                pass
        elif event.nickname == self.event_GeneChoice_J:
            try:
                # Attach anchors
                self.df_V_anchors = get_df_anchors_from_df_ref_genome(df)
            except Exception as e:
                print(e)
                pass



    def attach_anchors_from_files(self, fln_V_gene_CDR3_anchors=None, fln_J_gene_CDR3_anchors=None, sep=';'):
        """
        Add anchors to IgorModel_Parms from file, pandas dataframe or dictionary
        1. Get a dataframe from parms.Event_dict
        """

        try:
            self.attach_V_anchors_from_file(fln_V_gene_CDR3_anchors, sep=sep)
        except Exception as e:
            raise e

        try:
            self.attach_J_anchors_from_file(fln_J_gene_CDR3_anchors, sep=sep)
        except Exception as e:
            raise e


    def attach_V_anchors_from_file(self, fln_V_gene_CDR3_anchors, sep=';'):
        """
        Attach V anchors from file
        :param fln_V_gene_CDR3_anchors: IGoR's V anchors file
        """
        try:
            self.df_V_anchors = pd.read_csv(fln_V_gene_CDR3_anchors, sep=sep).set_index('gene')
        except Exception as e:
            raise e

    def attach_J_anchors_from_file(self, fln_J_gene_CDR3_anchors, sep=';'):
        """
        Attach J anchors from file
        :param fln_J_gene_CDR3_anchors: IGoR's J anchors file
        """
        try:
            self.df_J_anchors = pd.read_csv(fln_J_gene_CDR3_anchors, sep=sep).set_index('gene')
        except Exception as e:
            raise e

    def attach_V_anchors_from_Dataframe(self, df_V_anchors):
        # TODO: IN DEV
        #  if .set_index('gene')
        self.df_V_anchors = df_V_anchors

    def get_IgorRefGenome(self)->IgorRefGenome:
        """
        Return IgorRefGenome instance from events and df_V_anchors
        """
        df_genomics_dict = dict()
        if self.df_V_ref_genome is not None:
            df_genomics_dict['V'] = self.df_V_ref_genome

        if self.df_D_ref_genome is not None:
            df_genomics_dict['D'] = self.df_D_ref_genome

        if self.df_J_ref_genome is not None:
            df_genomics_dict['J'] = self.df_J_ref_genome

        ref_genome = IgorRefGenome.load_from_dataframe_genomics_dict(df_genomics_dict)
        return ref_genome


    # # FIXME: deprecated method
    # def get_ref_genome(self)->IgorRefGenome:
    #     # FIXME: IN DEV
    #     """Return IgorRefGenome genomes"""
    #     # Get all gene choice events from self (IgorModel_Parms
    #     tmp_dir = tempfile.TemporaryDirectory(prefix='ref_genome', dir='.')
    #     try:
    #         fln_genomicVs = tmp_dir.name + "/" + "genomicVs.fasta"
    #         fln_genomicDs = tmp_dir.name + "/" + "genomicDs.fasta"
    #         fln_genomicJs = tmp_dir.name + "/" + "genomicJs.fasta"
    #         fln_V_gene_CDR3_anchors = tmp_dir.name + "/" + "V_gene_CDR3_anchors.csv"
    #         fln_J_gene_CDR3_anchors = tmp_dir.name + "/" + "J_gene_CDR3_anchors.csv"
    #
    #         GeneChoice_list = [event for event in self.Event_list if event.event_type == 'GeneChoice']
    #         event_V = [event for event in GeneChoice_list if event.seq_type == 'V_gene'][0]
    #
    #         df_genetemplates = self.Event_dict[event_V.nickname].copy()
    #         df_genetemplates = df_genetemplates[['name', 'value']]
    #         df_genetemplates['id'] = df_genetemplates.index.get_level_values('id')
    #         df_genetemplates.set_index('name', inplace=True)
    #
    #         event_V.export_realizations_to_fasta(fln_genomicVs)
    #         event_J = [event for event in GeneChoice_list if event.seq_type == 'J_gene'][0]
    #         event_J.export_realizations_to_fasta(fln_genomicJs)
    #         try:
    #             event_D = [event for event in GeneChoice_list if event.seq_type == 'D_gene'][0]
    #             event_D.export_realizations_to_fasta(fln_genomicDs)
    #         except Exception as e:
    #             pass
    #
    #         # V Anchors
    #         try:
    #             self.df_V_ref_genome.to_csv(fln_V_gene_CDR3_anchors, columns=['name', 'anchor_index', 'function'],
    #                                         sep=';', index=False)
    #         except KeyError:
    #             self.df_V_ref_genome.to_csv(fln_V_gene_CDR3_anchors, columns=['name', 'anchor_index'],
    #                                         sep=';', index=False)
    #         except Exception as e:
    #             raise e
    #
    #         # J Anchors
    #         try:
    #             self.df_J_ref_genome.to_csv(fln_J_gene_CDR3_anchors, columns=['name', 'anchor_index', 'function'],
    #                                         sep=';', index=False)
    #         except KeyError:
    #             self.df_J_ref_genome.to_csv(fln_J_gene_CDR3_anchors, columns=['name', 'anchor_index'],
    #                                         sep=';', index=False)
    #         except Exception as e:
    #             raise e
    #
    #         # FIXME: IN DEV ADD ANCHORS
    #         # gene;anchor_index;function
    #
    #         ref_genome = IgorRefGenome.load_from_path(tmp_dir.name)
    #
    #     except Exception as e:
    #         raise e
    #     else:
    #         return ref_genome
    #     finally:
    #         tmp_dir.cleanup()
    #
    #
    #     return ref_genome

    def read_ErrorRate(self, ofile):
        lastPos = ofile.tell()
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

    def write_ref_genome_dir(self, ref_genome_dir_path):
        try:
            os.system("mkdir -p " + ref_genome_dir_path)
            ref_genome = self.get_IgorRefGenome()
            ref_genome.write_ref_genome_dir(ref_genome_dir_path)
        except Exception as e:
            raise e

    # FIXME: FINISH THIS METHOD
    def write_model_parms(self, filename=None, sep=';'):
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
        strSepChar = sep
        try:
            import os
            # print("AAAAAAAAAAAAAAA:", os.path.dirname(filename), filename)
            os.makedirs(os.path.dirname(filename), exist_ok=True)
        except Exception as e:
            pass
            # print("WARNING: write_model_parms path ", e)

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
        # TODO: WRITE FIRST EVENTS WITH HIGHEST PRIORITY AND LESS PARENTS
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
                df = df[['value', 'name']]
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
                ofile.write("%" + edge[0] + strSepChar + edge[1] + "\n")
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
        ofile.write("#" + self.ErrorRate_dict['error_type'] + "\n")
        ofile.write(self.ErrorRate_dict['error_values'] + "\n")

    def get_EventsNickname_list(self):
        return [event.nickname for event in self.Event_list]

    def get_EventsName_list(self):
        return [event.name for event in self.Event_list]

    def get_Event(self, event_nickname_or_name, by_nickname=True) -> IgorRec_Event:
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

    def get_Event_realization(self, event_nickname: str, index:Union[int, list]) -> IgorRec_Event:
        """Return event realization by event_nickname and index
        :param event_nickname: Nickname of event to get realization.
        :param index: Id of realization in event.
        :return: IgorRec_Event with nickname 'event_nickname' and id 'index'.
        """
        ps_realiz = self.parms.Event_dict[event_nickname].loc[index]
        return IgorEvent_realization.from_tuple(index, ps_realiz.value, ps_realiz.name)
        # return self.get_Event(event_nickname).get_realization(index)
        # elif isinstance(index, list):
        #     return [self.get_Event(event_nickname).get_realization(id) for id in index]

    def gen_EventDict_DataFrame(self):
        self.Event_dict = dict()
        # dictio = dict()
        for event in self.Event_list:
            # dictio[event.nickname] = event.get_realization_DataFrame()
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
        # self.G = self.G.reverse()

    def genPreMarginalDF(self):
        data = []
        for event in self.Event_list:
            # print (parms.dictNameNickname[event.name])
            # parms.Edges
            # tmpDict = dict()
            lista = []
            for edge in self.Edges:
                if edge[1] == event.name:
                    # print(parms.dictNameNickname[edge[0]])
                    # print(edge[0])
                    lista.append(self.dictNameNickname[edge[0]])
            tmpDict = {'event': event.nickname, 'priority': event.priority, 'Edges': lista}
            data.append(tmpDict)

        self.preMarginalDF = pd.DataFrame(data)  # .set_index('event')
        self.preMarginalDF['nEdges'] = self.preMarginalDF['Edges'].map(len)
        self.preMarginalDF.sort_values(['priority', 'nEdges'], ascending=[False, True])

    def genMarginalFile(self, model_marginals_file=None):
        self.genPreMarginalDF()
        # self.preMarginalDF
        if model_marginals_file == None:
            model_marginals_file = "model_marginals.txt"
        ofile = open(model_marginals_file, "w")
        for index, row in self.preMarginalDF.iterrows():
            nickname = row['event']
            ofile.write("@" + nickname + "\n")
            # DimEvent = len(parms.Event_dict[event.nickname])
            # DimEdges = len(parms.Edges_dict[event.nickname])

            DimEvent = len(self.Event_dict[nickname])
            strDimLine = "$Dim["
            DimList = []
            if row['nEdges'] == 0:
                strDimLine = strDimLine + str(DimEvent)
                strDimLine = strDimLine + "]"
            else:
                for evNick in row['Edges']:  # parms.Edges_dict[event.nickname]:
                    Dim = len(self.Event_dict[evNick])
                    strDimLine = strDimLine + str(Dim) + ","
                    DimList.append(Dim)
                strDimLine = strDimLine + str(DimEvent)
                strDimLine = strDimLine + "]"
            ofile.write(strDimLine + "\n")

            lista = row['Edges']  # self.Event_dict[nickname]
            for indices in np.ndindex(tuple(DimList)):
                # print indices
                strTmp = "#"
                for ii in range(len(lista)):
                    strTmp = strTmp + "[" + lista[ii] + "," + str(indices[ii]) + "]"
                    if not (ii == len(lista) - 1):
                        strTmp = strTmp + ","
                ofile.write(strTmp + "\n")
                ofile.write("%")
                unifProb = (1. / DimEvent)
                for jj in range(DimEvent):
                    ofile.write(str(unifProb))
                    if not (jj == DimEvent - 1):
                        ofile.write(",")
                ofile.write("\n")

        ofile.close()

    def plot_Graph(self, ax=None, **kwargs):  # FIXME: ALLOW the possibility to pass an ax like ax=None):
        """Return a plot of the bayesian network """
        # if ax is None:

        pos = nx.spring_layout(self.G)
        # priorities up
        prio_dict = dict()
        for event in self.Event_list:
            if not (event.priority in prio_dict):
                prio_dict[event.priority] = list()
            prio_dict[event.priority].append(event)
        # print(str(prio_dict))
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
                    # print("matplotlib")
                    import matplotlib.pyplot as plt
                    fig, ax = plt.subplots()
                ax.set_aspect('equal')
                dict_nickname_seq_type = self.get_event_dict('nickname', 'seq_type')
                colors_list = list(
                    map(lambda x: Igor_seq_type_color_dict[dict_nickname_seq_type[x]], list(self.G.nodes())))

                nx.draw(self.G, pos=pos, ax=ax, with_labels=True, arrows=True, arrowsize=20,
                        node_size=800, alpha=0.5, font_size=10, font_weight='bold',
                        nodelist=self.G.nodes(), node_color=colors_list)  # FIXME: make a better plot: cutting edges.

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
    def realiz_dict_from_scenario(self, scenario):  # IgorScenario):
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
                    realizations_dict[nickname_key] = event.realizations[scenario[nickname_key]]
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

    def __getitem__(self, item):
        return self.marginals_dict[item]
    #  @d_3_del
    #  $Dim[3,21,21]
    #  #[d_gene,0],[d_5_del,0]
    #  %0,0,0,1.6468e-08,0.00482319,1.08101e-09,0.0195311,0.0210679,0.0359338,0.0328678,2.25686e-05,4.97463e-07,0,9.31048e-08,1.01642e-05,0.000536761,0.0260845,0.0391021,0.319224,0.289631,0.211165
    #  #[d_gene,0],[d_5_del,1]
    #  %0,0,6.86291e-08,2.00464e-09,0.00163832,2.02919e-06,0.0306066,0.0126832,0.000872623,0.016518,0.00495292,0.000776747,4.45576e-05,0.000667902,0.00274004,0.00435049,0.300943,0.182499,0.13817,0.302534,0

    @classmethod
    def make_uniform_from_parms(cls, parms: IgorModel_Parms):
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
            print("Reading Marginals filename from: ", filename)
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

        # return marginals_dict, network_dict

    def initialize_uniform_event_from_model_parms(self, event_nickname, parms: IgorModel_Parms):
        event = parms.get_Event(event_nickname)
        if event.event_type == 'DinucMarkov':
            # do something
            dimension = len(parms.get_Event(event.nickname).realizations)
            narr = np.ones(dimension * dimension) / (dimension) # * dimension)
            self.marginals_dict[event.nickname] = narr
        else:
            dimensions = [len(parms.get_Event(strEvent).realizations) for strEvent in
                          self.network_dict[event.nickname]]
            # print(dimensions[-1])
            narr = np.ones(dimensions) / dimensions[-1]

            self.marginals_dict[event.nickname] = narr

    def update_network_dict_from_model_parms(self, parms: IgorModel_Parms):
        self.network_dict = dict()
        for key, value in parms.Edges_dict.items():
            self.network_dict[key] = value + [key]

    def initialize_uniform_from_model_parms(self, parms: IgorModel_Parms):
        """
        Update with uniform distribution IgorModel_Marginals from IgorModel_Parms object
        :param parms: IgorModel_Parms object
        """
        self.update_network_dict_from_model_parms(parms)

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
        parms = model_parms  # IgorModel_Parms(model_parms_file=model_parms_file)

        if filename is None:
            filename = "tmp_mdl_marginals.txt"

        try:
            import os
            os.makedirs(os.path.dirname(filename), exist_ok=True)
        except Exception as e:
            pass
            # print("WARNING: IgorModel_Marginals.write_model_marginals path ", e)

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


class IgorModel:
    """
    :class: IgorModel,
    :param model_parms_file: Path of IGoR's model parms file
    :type model_parms_file: Union[None, str, Path]=None
    :param model_marginals_file: Path of IGoR's model marginals file
    :type model_marginals_file: Union[None, str, Path]=None
    :param parms: IgorModel_Parms instance
    :type parms: Union[None, IgorModel_Parms] = None
    :param marginals: IgorModel_Marginals instance
    :type marginals: Union[None, IgorModel_Marginals] = None
    """
    def __init__(self, model_parms_file: Union[None, str, Path]=None,
                 model_marginals_file: Union[None, str, Path]=None,
                 parms: Union[None, IgorModel_Parms] = None,
                 marginals: Union[None, IgorModel_Marginals] = None,
                 fln_V_gene_CDR3_anchors:Union[None, str, Path]=None,
                 fln_J_gene_CDR3_anchors:Union[None, str, Path]=None):
        """Constructor method
        """
        self.parms = None
        self.marginals = None

        if parms is not None:
            self.parms = parms

        if marginals is not None:
            self.marginals = marginals

        if self.parms is None:
            self.parms = IgorModel_Parms()
        if self.marginals is None:
            self.marginals = IgorModel_Marginals()

        self.genomic_dataframe_dict = dict() # FIXME: CHANGE THIS AS PROPERTY
        self.xdata = dict()
        self.factors = list()
        self.metadata = dict()
        self.specie = ""
        self.chain = ""

        self.anchors_CDR3_V = None
        self.anchors_CDR3_J = None

        self.Pmarginal = dict()

        # FIXME: But since DB is in refactor keep it for the moment
        self.BestScenariosHeaderList = list()  # This is a ordered list store the nicknames of events in the header of the file
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

        if (not (fln_V_gene_CDR3_anchors is None)) and (not (fln_J_gene_CDR3_anchors is None)):
            try:
                self.parms.attach_anchors_from_files(fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                                                 fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors, sep=';')
            except KeyError as e:
                # OLGA
                self.parms.attach_anchors_from_files(fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                                                     fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors, sep=',')
                pass
            except Exception as e:
                raise e
            if self.parms.event_GeneChoice_V is not None:
                self.genomic_dataframe_dict['V'] = self.parms.df_V_ref_genome

            if self.parms.event_GeneChoice_D is not None:
                self.genomic_dataframe_dict['D'] = self.parms.df_D_ref_genome

            if self.parms.event_GeneChoice_J is not None:
                self.genomic_dataframe_dict['J'] = self.parms.df_J_ref_genome


    def __getitem__(self, key):
        return self.xdata[key]

    def __str__(self):
        return ".xdata" + str(self.get_events_nicknames_list())

    def get_Event_value(self, event_nickname, index):
        return self.parms[event_nickname].value.loc[index]

    def get_df_events_rearrangement(self):
        """
        Return event dataframe sorted for V(D)J rearrangement
        """
        try:
            import pandas as pd

            list_events_properties = list()
            for ev_nickname in self.xdata.keys():
                list_events_properties.append(self.xdata[ev_nickname].attrs)

            df_events = pd.DataFrame(list_events_properties)
            df_events.sort_values(by='seq_type').reset_index()

            Igor_dict_order_seq_type = {'V_gene': 0, 'VD_genes': 1, 'D_gene': 2, 'DJ_gene': 3, 'VJ_gene': 4, 'J_gene': 5,
                                        'VDJ_genes': 6, 'Undefined_gene': 7}
            Igor_dict_order_event_type = {'GeneChoice': 0, 'Deletion': 1, 'Insertion': 2, 'DinucMarkov': 3, 'Undefined': 4}
            Igor_dict_order_seq_side = {'Five_prime': 0, 'Three_prime': 1, 'Undefined_side': 2}

            df_events['n_seq_type'] = df_events['seq_type'].map(Igor_dict_order_seq_type)
            df_events['n_event_type'] = df_events['event_type'].map(Igor_dict_order_event_type)
            df_events['n_seq_side'] = df_events['seq_side'].map(Igor_dict_order_seq_side)
            df_rearrengement_events = df_events.sort_values(
                by=['n_seq_type', 'n_event_type', 'n_seq_side']).reset_index().drop(columns=['index'])
            df_rearrengement_events.drop(columns=['n_seq_type', 'n_event_type', 'n_seq_side'], inplace=True)
            return df_rearrengement_events
        except Exception as e:
            raise e

    @property
    def event_GeneChoice_nickname_list(self):
        try:
            return list(map(lambda x: x.nickname, self.parms.event_GeneChoice_list))
        except Exception as e:
            raise e

    @property
    def event_Insertion_nickname_list(self):
        try:
            return list(map(lambda x: x.nickname, self.parms.event_Insertion_list))
        except Exception as e:
            raise e

    @property
    def event_Deletion_nickname_list(self):
        try:
            return list(map(lambda x: x.nickname, self.parms.event_Deletion_list))
        except Exception as e:
            raise e

    @property
    def event_DinucMarkov_nickname_list(self):
        try:
            return list(map(lambda x: x.nickname, self.parms.event_DinucMarkov_list))
        except Exception as e:
            raise e

    @property
    def event_GeneChoice_V_nickname(self):
        return self.parms.event_GeneChoice_V.nickname

    @property
    def event_GeneChoice_D_nickname(self):
        return self.parms.event_GeneChoice_D.nickname

    @property
    def event_GeneChoice_J_nickname(self):
        return self.parms.event_GeneChoice_J.nickname

    @property
    def ErrorRate_dict(self):
        try:
            return self.parms.ErrorRate_dict
        except AttributeError:
            return None
        except Exception as e:
            raise e

    @property
    def Pconditionals(self):
        return self.xdata


    @property
    def V_anchors(self):
        # return self.__V_anchors
        return self.genomic_dataframe_dict['V']['anchor_index'].dropna().astype(int)

    @property
    def J_anchors(self):
        return self.genomic_dataframe_dict['J']['anchor_index'].dropna().astype(int)


    def V_anchor(self, id: int):
        try:
            return self.V_anchors.loc[id]
        except Exception as e:
            raise e


    def J_anchor(self, id: int):
        try:
            return self.J_anchors.loc[id]
        except Exception as e:
            raise e

    @classmethod
    def make_default_model_from_IgorRefGenome(cls, genomes:IgorRefGenome):
        cls = IgorModel()
        try:
            if genomes.df_genomicDs is None:
                cls = IgorModel.make_model_default_VJ_from_dataframes(genomes.df_V_ref_genome,
                                                          genomes.df_J_ref_genome)
            else:
                cls = IgorModel.make_model_default_VDJ_from_dataframes(genomes.df_V_ref_genome,
                                                           genomes.df_genomicDs,
                                                           genomes.df_J_ref_genome)
            # Define anchors
            # cls.anchors_CDR3_V = genomes.df_V_ref_genome['anchor_index']
            # cls.anchors_CDR3_J = genomes.df_J_ref_genome['anchor_index']
        except Exception as e:
            raise e
        else:
            return cls

    @classmethod
    def make_model_default_VDJ_from_dataframes(cls,
                                               df_V_ref_genome: Union[pd.DataFrame],
                                               df_D_ref_genome: Union[pd.DataFrame],
                                               df_J_ref_genome: Union[pd.DataFrame],
                                               lims_deletions=None, lims_insertions=None):
        """
        Returns IgorModel with uniform ditribution with for the default
        :param df_V_ref_genome:Union[pd.DataFrame],
        :param df_D_ref_genome:Union[pd.DataFrame],
        :param df_J_ref_genome:Union[pd.DataFrame]
        :return: IgorModel object
        """
        try:
            cls = IgorModel()
            cls.parms = IgorModel_Parms.make_default_VDJ(df_V_ref_genome,
                                                         df_D_ref_genome,
                                                         df_J_ref_genome,
                                                         lims_deletions=lims_deletions,
                                                         lims_insertions=lims_insertions)
            cls.parms.Event_list = cls.parms.get_Event_list_sorted()
            cls.marginals.initialize_uniform_from_model_parms(cls.parms)
            cls.generate_xdata()
            return cls
        except Exception as e:
            raise e

    @classmethod
    def make_model_default_VJ_from_dataframes(cls,
                                               df_V_ref_genome: Union[pd.DataFrame],
                                               df_J_ref_genome: Union[pd.DataFrame],
                                               lims_deletions=None, lims_insertions=None):
        """
        Returns IgorModel with uniform ditribution with for the default
        :param df_V_ref_genome:Union[pd.DataFrame],
        :param df_D_ref_genome:Union[pd.DataFrame],
        :param df_J_ref_genome:Union[pd.DataFrame]
        :return: IgorModel object
        """
        cls = IgorModel()
        cls.parms.make_default_VJ(df_V_ref_genome, df_J_ref_genome, lims_deletions=lims_deletions,
                                   lims_insertions=lims_insertions)
        cls.parms.Event_list = cls.parms.get_Event_list_sorted()
        cls.marginals.initialize_uniform_from_model_parms(cls.parms)
        cls.generate_xdata()
        return cls



    # TODO: finish this method to load model with default installed igor.
    @classmethod
    def load_default(cls, IgorSpecie, IgorChain, modelpath=None, ref_genome_path=None):  # rcParams['paths.igor_models']):
        """
        :return: IgorModel loaded with the default location for specie and chain
        """

        if IgorChain in igor_option_path_dict:
            IgorChain = igor_option_path_dict[IgorChain]

        flnModelParms, flnModelMargs = get_default_models_paths_species_chain(IgorSpecie, IgorChain, modelpath=modelpath)

        fln_dict = get_default_fln_dict_ref_genomes_species_chain(IgorSpecie, IgorChain, modelspath=modelpath, ref_genome_path=ref_genome_path)

        # if modelpath is None:
        #     try:
        #         modelpath = run_igor_datadir() + "/models"
        #     except Exception as e:
        #         print("ERROR: getting default igor datadir.", e)
        #
        #
        #
        #
        # IgorModelPath = modelpath+"/"+IgorSpecie+"/"+IgorChain+"/"
        # print("Loading default IGoR model from path : ", IgorModelPath)
        # # FIXME: FIND A WAY TO GENERALIZE THIS WITH SOMEKIND OF STANDARD NAME
        # flnModelParms = IgorModelPath + "models/model_parms.txt"
        # flnModelMargs = IgorModelPath + "models/model_marginals.txt"
        # print("Parms filename: ", flnModelParms)
        # print("Margs filename: ", flnModelMargs)
        # print("-"*50)

        #        IgorRefGenomePath = IgorModelPath+"ref_genome/"
        #        flnVGeneTemplate = IgorRefGenomePath+"genomicVs.fasta"
        #        flnDGeneTemplate = IgorRefGenomePath+"genomicDs.fasta"
        #        flnJGeneTemplate = IgorRefGenomePath+"genomicJs.fasta"
        #
        #        flnVGeneCDR3Anchors = IgorRefGenomePath+"V_gene_CDR3_anchors.csv"
        #        flnJGeneCDR3Anchors = IgorRefGenomePath+"J_gene_CDR3_anchors.csv"
        cls = IgorModel(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs,
                        fln_V_gene_CDR3_anchors=fln_dict['fln_V_gene_CDR3_anchors'],
                        fln_J_gene_CDR3_anchors=fln_dict['fln_J_gene_CDR3_anchors'])

        try:
            fln_ref_genomes_dict = get_default_fln_dict_ref_genomes_species_chain(IgorSpecie, IgorChain, modelspath=modelpath,
                                                                                  ref_genome_path=ref_genome_path)
            fln_V_gene_CDR3_anchors = fln_ref_genomes_dict['fln_V_gene_CDR3_anchors']
            fln_J_gene_CDR3_anchors = fln_ref_genomes_dict['fln_J_gene_CDR3_anchors']
            cls.parms.attach_anchors_from_files(fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors)
        except:
            print("No anchors attached.")
            pass
        cls.specie = IgorSpecie
        cls.chain = IgorChain

        return cls

    @classmethod
    def load_from_txt(cls, model_parms_file: str, model_marginals_file: Union[None, str] = None):
        """
        load model from txt path files (model_parms.txt and model_marginals.txt
        if model marginals is not specified a uniform distribution is loaded.
        """
        try:
            cls = IgorModel()
            cls.read_model_from_txt(model_parms_file, model_marginals_file)
            return cls
        except Exception as e:
            e_message = "IgorModel.load_from_txt : model_parms_file = " + str(model_parms_file)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    # FIXME: THIS COULD BE CONFUSING WITH igor_model_dir_path
    @classmethod
    def load_from_directory(cls, model_files_dir):
        """
        return a IgorModel from directory with default names 'model_parms.txt' and 'model_marginals.txt'
        """
        try:
            cls = IgorModel()
            cls.read_model_from_directory(model_files_dir)
            return cls
        except Exception as e:
            e_message = "IgorModel.load_from_directory : model_files_dir = " + str(model_files_dir)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    @classmethod
    def load_from_parms_marginals_object(cls, mdl_parms: IgorModel_Parms,
                                         mdl_marginals: Union[None, IgorModel_Marginals] = None):
        """
        Load IgorModel from IgorModel_Parms and IgorModel_Marginals instances.
        If IgorModel_Marginals not provided, a uniform distribution is provided for marginals.
        :return : IgorModel instance
        """
        try:
            cls = IgorModel()
            cls.parms = mdl_parms
            if mdl_marginals is None:
                mdl_marginals = IgorModel_Marginals()
                mdl_marginals.initialize_uniform_event_from_model_parms(cls.parms)

            cls.marginals = mdl_marginals
            cls.generate_xdata()
            return cls
        except Exception as e:
            raise e

    # FIXME:
    @classmethod
    def load_from_networkx(cls, IgorSpecie, IgorChain):
        """
        :return IgorModel loaded with the default location for specie and chain
        """
        cls = IgorModel(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
        return cls

    def generate_xdata(self):
        """Load model dictionary with xarray structures from parms and marginals instances"""
        # TODO: CHANGE TO A QUERY IN PARMS WITH HIGHT PRIORITY FIRST AND LESS NUMBER OF PARENTS

        try:
            Event_Genechoice_List = ['v_choice', 'j_choice', 'd_gene']
            Event_Dinucl_List = ['vd_dinucl', 'dj_dinucl', 'vj_dinucl']
            Event_Insertion_List = ['vd_ins', 'dj_ins', 'vj_ins']
            Event_Deletion_List = ['v_3_del', 'j_5_del', 'd_3_del', 'd_5_del']

            for key in self.marginals.marginals_dict:
                event = self.parms.get_Event(key)

                if event.event_type == 'DinucMarkov':
                    # if key in Event_Dinucl_List:
                    # FIXME: reshape by dimension in marginal
                    self.xdata[key] = xr.DataArray(self.marginals.marginals_dict[key].reshape(4, 4), \
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
                    # print "key: ", key, self.xdata[key].dims

                    for strDim in self.xdata[key].dims:
                        self.xdata[key][strDim] = range(len(self.xdata[key][strDim]))
                        if strDim in Event_Genechoice_List:
                            # print strDim
                            # labels = self.parms.Event_dict[strDim]['name'].map(genLabel).values # FIXME: use the exact name defined in model_parms
                            labels = self.parms.Event_dict[strDim]['name'].values
                            strCoord = 'lbl__' + strDim
                            self.xdata[key][strCoord] = (strDim, labels)  # range(len(self.xdata[key][coord]))

                            sequences = self.parms.Event_dict[strDim]['value'].values
                            strCoord = 'seq__' + strDim
                            self.xdata[key][strCoord] = (strDim, sequences)

                        elif not (strDim in Event_Dinucl_List):
                            labels = self.parms.Event_dict[strDim]['value'].values
                            strCoord = 'lbl__' + strDim
                            self.xdata[key][strCoord] = (strDim, labels)  # range(len(self.xdata[key][coord]))
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

            # TODO: ADD genomic_dataframe_dict when generate_xdata is call
            if self.parms.event_GeneChoice_V is not None:
                self.genomic_dataframe_dict['V'] = self.parms.df_V_ref_genome
            if self.parms.event_GeneChoice_D is not None:
                self.genomic_dataframe_dict['D'] = self.parms.df_D_ref_genome
            if self.parms.event_GeneChoice_J is not None:
                self.genomic_dataframe_dict['J'] = self.parms.df_J_ref_genome

            self.generate_Pmarginals()
        except Exception as e:
            raise e

    def read_model_from_txt(self, model_parms_file: str,
                            model_marginals_file: Union[None, str] = None):
        """
        Read model from model_parms.txt and model_marginals.txt.
        :param model_parms_file: Path to model parms txt file.
        :param model_marginals_file: Path to model marginals txt file.
        """

        try:
            self.parms = IgorModel_Parms()
            self.parms.read_model_parms(model_parms_file)

            try:
                self.marginals = IgorModel_Marginals()
                if model_marginals_file is None:
                    # make a marginals uniform from parms and dont write it
                    self.marginals.initialize_uniform_from_model_parms(parms=self.parms)
                else:
                    self.marginals.read_model_marginals(model_marginals_file)

            except Exception as e:
                e_message = "IgorModel.read_model_from_txt : " + str(model_marginals_file)
                import sys
                raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

            self.generate_xdata()

        except Exception as e:
            e_message = "IgorModel.read_model_from_txt : " + str(model_parms_file)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def read_model_from_directory(self, model_files_dir):
        """
        Read model from model files directory.
        :param model_files_dir: Path to model directory.
        """
        try:
            model_parms_file = model_files_dir + "/models/model_parms.txt"
            model_marginals_file = model_files_dir + "/models/model_marginals.txt"
            self.read_model_from_txt(model_parms_file, model_marginals_file)
        except IOError as e:
            print("model_parms.txt not found in ", model_files_dir)
            try:
                print("looking for model_params.txt in ", model_files_dir, "(OLGA name)")
                model_parms_file = model_files_dir + "/model_params.txt"
                model_marginals_file = model_files_dir + "/model_marginals.txt"
                self.read_model_from_txt(model_parms_file, model_marginals_file)
            except Exception as e:
                raise e

        # TODO: LOOK FOR ANCHORS
        try:
            fln_V_gene_CDR3_anchors = model_files_dir + '/ref_genome/V_gene_CDR3_anchors.csv'
            fln_J_gene_CDR3_anchors = model_files_dir + '/ref_genome/J_gene_CDR3_anchors.csv'
            self.parms.attach_anchors_from_files(fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                                                 fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors, sep=';')
        except IOError as e:
            fln_V_gene_CDR3_anchors = model_files_dir + '/V_gene_CDR3_anchors.csv'
            fln_J_gene_CDR3_anchors = model_files_dir + '/J_gene_CDR3_anchors.csv'
            self.parms.attach_anchors_from_files(fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                                                 fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors, sep=',')

        except Exception as e:
            print(fln_V_gene_CDR3_anchors + " not found!")
            print(fln_J_gene_CDR3_anchors + " not found!")
            try:
                fln_V_gene_CDR3_anchors = model_files_dir + '/V_gene_CDR3_anchors.csv'
                fln_J_gene_CDR3_anchors = model_files_dir + '/J_gene_CDR3_anchors.csv'
                self.parms.attach_anchors_from_files(fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                                                     fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors, sep=',')
            except Exception as ee:
                print('neither: ')
                print(fln_V_gene_CDR3_anchors + "not found!")
                print(fln_J_gene_CDR3_anchors + "not found!")
                print("WARNING: Not anchors for model, use IgorParms.attach_anchor_from_files to attach correct anchors!")
                pass
        else:
            print("Anchors loaded from " + fln_V_gene_CDR3_anchors + " and " + fln_J_gene_CDR3_anchors)

        if self.parms.event_GeneChoice_V is not None:
            self.genomic_dataframe_dict['V'] = self.parms.df_V_ref_genome

        if self.parms.event_GeneChoice_D is not None:
            self.genomic_dataframe_dict['D'] = self.parms.df_D_ref_genome

        if self.parms.event_GeneChoice_J is not None:
            self.genomic_dataframe_dict['J'] = self.parms.df_J_ref_genome

    def get_zero_xarray_from_list(self, strEvents_list: list):
        """Get xarray with labels and dimensions for strEvents_list
        :param strEvents_list: list of events nickname.
        :return: xarray with dimensions and coordinates with zero as values.
        """
        # strEvents_list = ['v_choice', 'j_choice']
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

    def get_observable_xarray_from_function(self, observable_func, variables_tuple_list):
        """
        observable_func(x,y,z)
        variables_tuple_list = [('v_choice', 'id'), ('vd_ins', 'value'), ('v_3_del', 'value')]
        """

        # 1. Get the nicknames only
        strEvents_list = [var[0] for var in variables_tuple_list]
        strEvents_tuple = tuple(strEvents_list)

        # Use model parms to create xarray with values
        da_shape_list = [len(self.parms.Event_dict[str_event_nickname]) for str_event_nickname in strEvents_list]
        da_shape_tuple = tuple(da_shape_list)

        data_arr = np.zeros(da_shape_tuple)

        # create meshgrid
        tmp_list_2_broadcast = list()
        for event_nickname, choose_type in variables_tuple_list:
            if choose_type == 'id':
                tmp_list_2_broadcast.append(self.parms.Event_dict[event_nickname].index.values)
            else:
                tmp_list_2_broadcast.append(self.parms.Event_dict[event_nickname][choose_type].values)

        if len(strEvents_list) > 1:
            variables_mesh = np.meshgrid(*tmp_list_2_broadcast)
            vect_observable_func = np.vectorize(observable_func)
            data_arr = vect_observable_func(*variables_mesh)

        da = xr.DataArray(data_arr, dims=strEvents_tuple)
        # da.assign_coords()
        for strDim in strEvents_list:
            da[strDim] = range(len(da[strDim]))
        # self.xdata[key][strDim] = range(len(self.xdata[key][strDim]))

        # for event_nickname in strEvents_list:
        #     da[event_nickname] = self.parms.Event_dict[event_nickname].index.values
        #     labels = self.parms.Event_dict[event_nickname]['name'].values
        #     strCoord = 'lbl__' + event_nickname
        #     da[strCoord] = (event_nickname, labels)

        return da

    def get_observable_from_df_scenarios(self, observable_function, df_scenarios:pd.DataFrame):
        """
        Return a pandas series with the calculated observable over the df_scenarios dataframe.
        :param observable_function: This function should use the varibles with self.realization
        :param df_scenarios: Scenarios dataframe loaded with self.get_dataframe_scenarios.
        """
        return df_scenarios.apply(lambda row: observable_function(row), axis=1)


    def get_realization_value_from_df_scenarios(self, df_scenarios, event_nickname):
        return self.get_observable_from_df_scenarios(lambda x: self.realization(x, event_nickname).value, df_scenarios)

    def get_ones_xarray_from_list(self, strEvents_list: list):
        """Get xarray with labels and dimensions for strEvents_list
        :param strEvents_list: list of events nickname.
        :return: xarray with dimensions and coordinates with one as values.
        """
        return xr.ones_like(self.get_zero_xarray_from_list(strEvents_list))

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
        """
        Variable elimination to get probabily marginals of event strEvent
        :parm strEvent: event nickname
        """
        # FIXME: use xdata instead of self.parms
        try:
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
        except Exception as e:
            raise e

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

    def get_P_joint(self, not_sum_out_nickname_list:list):
        """
        Return xarray DataArray of the joint probability of nickname event list
        :param not_sum_out_nickname_list: list of nickname events to get
        the probability joint distribution. DinucMarkov events not accepted.
        """
        try:
            # TODO: CHECK THAT NO DINUCMARKOV EVENT PRESENT IN not_not_sum_out_nickname_list
            # use all or any
            sorted_events = self.parms.get_Event_list_sorted()
            # print(list(map(lambda x: x.nickname, sorted_events)))
            # not_sum_out_nickname_list = ['j_choice', 'd_gene']  # Events to get joint distribution
            sorted_nicknames_to_sum_out = [event.nickname for event in sorted_events if (
                        not event.event_type == "DinucMarkov" and not event.nickname in not_sum_out_nickname_list)]

            factors = self.VE_get_Pmarginals_initial_factors()
            for event_nickname_to_eliminate in sorted_nicknames_to_sum_out:
                factors = self.VE_get_factors_by_sum_out_variable(event_nickname_to_eliminate, factors)

            # Now multiply the remaining factors to get the marginal.
            P_joint = 1
            for factor in factors:
                P_joint = P_joint * factor

            return P_joint
        except Exception as e:
            raise e


    def get_mutual_information_events(self, event_nickname1, event_nickname2):
        """ Return xarray with
        """
        try:
            if not nx.d_separated(self.parms.G, {event_nickname1}, {event_nickname2}, {}):
                da_P_x_y = self.get_P_joint([event_nickname1, event_nickname2])
                da_P_x = self.Pmarginal[event_nickname1]
                da_P_y = self.Pmarginal[event_nickname2]

                da_P_x_times_P_y = (da_P_x * da_P_y)
                da_log_P_ratio = xr.zeros_like(da_P_x_y)
                da_log_P_ratio.values = np.nan_to_num(
                    np.log2(da_P_x_y / da_P_x_times_P_y), nan=0.0, neginf=0.0
                )
                return xr.dot(da_P_x_y, da_log_P_ratio)
                # return get_D_KL_from_xarray(da_P_x_y, da_P_x, da_P_y)
            else:
                return 0
        except Exception as e:
            raise e

    def get_mutual_information(self):
        """
        Return xarray with mutual information
        """
        # sorted_events = self.parms.get_Event_list_sorted()
        # GeneChoice_list = [event for event in sorted_events if event.event_type == 'GeneChoice']
        # Insertion_list = [event for event in sorted_events if event.event_type == 'Insertion']
        # Deletion_list = [event for event in sorted_events if event.event_type == 'Deletion']
        # DinucMarkov_list = [event for event in sorted_events if event.event_type == 'DinucMarkov']
        # events_no_DinucMarkov = [event for event in sorted_events if not event.event_type == 'DinucMarkov']

        try:
            np.seterr(divide='ignore')

            dict_nickname_event_type = self.parms.get_event_dict('nickname', 'event_type')
            dict_events = {key: val for key, val in dict_nickname_event_type.items() if val != 'DinucMarkov'}
            event_lista_nicknames = list(dict_events.keys())
            data_0 = np.zeros((len(event_lista_nicknames), len(event_lista_nicknames)))
            da_mi = xr.DataArray(data_0, dims=('x', 'y'), coords={'x': event_lista_nicknames, 'y': event_lista_nicknames})
            da_mi.name = 'mutual_information'


            import itertools
            for event_nickname_x, event_nickname_y in itertools.combinations_with_replacement(event_lista_nicknames, 2):
            # for event_nickname_x, event_nickname_y in itertools.product(event_lista_nicknames, event_lista_nicknames):
                if event_nickname_x != event_nickname_y:
                    mi = self.get_mutual_information_events(event_nickname_x, event_nickname_y)
                    da_mi.loc[{"x": event_nickname_x, "y": event_nickname_y}] = mi
                    da_mi.loc[{"x": event_nickname_y, "y": event_nickname_x}] = mi
                else:
                    mi = 0
                    da_mi.loc[{"x": event_nickname_x, "y": event_nickname_y}] = mi
                # print(event_nickname_x, event_nickname_y, mi)
            return da_mi
        except Exception as e:
            raise e
        finally:
            np.seterr(divide='warn')


    def get_entropy_event(self, event_nickname):
        """
        Returns entropy of event. If event has parents it returns the conditional
        entropy of the event with this parents
        H(X|Y) = - \sum_{x,y} p(x,y) \log{p(x|y)}
        :param event_nickname: Event nickname to calculate entropy
        """
        try:

            da_CP_event = self.Pconditionals[event_nickname]

            da_log2_event = xr.zeros_like(da_CP_event)
            da_log2_event.values = np.nan_to_num(np.log2(da_CP_event.values), neginf=0, nan=0)

            if da_CP_event.attrs['event_type'] == 'DinucMarkov':
                # FIXME: IN DEV
                H_event = -(da_CP_event * da_log2_event).sum(dim='y')
                return H_event

            elif (len(da_CP_event.attrs['parents']) > 0):
                H_event_given_parents = -(da_CP_event * da_log2_event).sum(dim=event_nickname)
                P_joint = self.get_P_joint(da_CP_event.attrs['parents'])
                H_event_conditional = (P_joint*H_event_given_parents).sum()
                return H_event_conditional
            else:
                H_event = -(da_CP_event * da_log2_event).sum(dim=event_nickname)
                return H_event
                # log2_Pmarginal = np.log2(self.Pconditionals[event_nickname])
                # log2_Pmarginal.values = np.nan_to_num(log2_Pmarginal.values, neginf=0, nan=0)
                # da_entropy = - xr.dot(self.Pmarginal[event_nickname], log2_Pmarginal)
                # return da_entropy
        except Exception as e:
            raise e

    def get_conditional_entropy_dinucl_function_l_ins(self, event_nickname_dinucl:str):
        """
        Return a function that depends on insertion length H(P_{m_i|l}) = H(p_{ss}) - (l-1) \sum_m p_{ss}(m) \sum_n T(n|m) \log2{T(n|m)}
        where p_{ss} is the stationary state calculated for T (eigenvector for eigenvalue 1)
        T(n|m) jump matrix from m to n, where m and n are nucleotides.
        """
        da_dinucl = self.Pconditionals[event_nickname_dinucl]
        if da_dinucl.attrs['event_type'] == 'DinucMarkov':
            # 1. Calculate the stationary distribution
            p_ss = np.real(get_P_stationary_state_dinucl(da_dinucl))
            H_dinucl = self.get_entropy_event(event_nickname_dinucl)

            # 2. Calculate the entropy of the stationary distribution
            # get_entropy_event(event_nickname_dinucl)
            H_p_ss = -(p_ss * np.log2(p_ss)).sum()

            # 3. Calculate the entropic contributions of lenght l to the entropy
            matmul_H_dinucl_p_ss = np.matmul(H_dinucl.values, p_ss)
            def tmp_function(l_ins):
                # 4. Return  H(P_{m_i}|l)
                if l_ins > 0:
                    return H_p_ss + (l_ins - 1) * matmul_H_dinucl_p_ss
                elif l_ins == 0:
                    return H_p_ss
                else:
                    return None

            return np.vectorize(tmp_function)
        else:
            return None


    def get_df_Insertion_entropy_contribution(self):
        """
        H(P({m_i})) = H(P_{ins}) - \sum_l P_{ins}(l) \sum_{m_i |l} H( P(m_i |l) )
        Return entropy contributions of insertions considering the diversity of nucleotides.
        """
        # 1. Find the Insertions events
        list_Insertion_events_nickname = list() # vd_ins, dj_ins
        list_DinucMarkov_events_nickname = list()
        for event_nickname in self.Pconditionals.keys():
            if self.Pconditionals[event_nickname].attrs['event_type'] == 'Insertion':
                list_Insertion_events_nickname.append(event_nickname)
            elif self.Pconditionals[event_nickname].attrs['event_type'] == 'DinucMarkov':
                list_DinucMarkov_events_nickname.append(event_nickname)

        # 2. Associate DinucMarkov with Insertion event by seq_type
        dict_Insertion_DinucMarkov = dict() # dict_Insertion_DinucMarkov['vd_ins'] =  vd_dinucl
        for insertion_event_nickname in list_Insertion_events_nickname:
            insertion_seq_type = self.Pconditionals[insertion_event_nickname].attrs['seq_type']
            for dinucl_event_nickname in list_DinucMarkov_events_nickname:
                if insertion_seq_type == self.Pconditionals[dinucl_event_nickname].attrs['seq_type']:
                    dict_Insertion_DinucMarkov[insertion_event_nickname] = dinucl_event_nickname

        # 3. Calculate the stationary state of DinucMarkov events
        # H(P({m_i})) = H(P_{ins}) - \sum_l P_{ins}(l) \sum_{m_i |l} H( P(m_i |l) )

        # 4. Finally Get the entropic contribution
        list_Insertion_entropy = list()
        for insertion_event_nickname in list_Insertion_events_nickname:
            vf_H_dinucl_given_l = self.get_conditional_entropy_dinucl_function_l_ins(
                dict_Insertion_DinucMarkov[insertion_event_nickname])
            H_P_mi_l = vf_H_dinucl_given_l(self.parms[insertion_event_nickname]['value'].values)
            # FIXME: WHICH ONE IS THE BEST Pmarginal or Pconditional?
            # mdl.get_entropy_event('vd_ins') + np.dot(H_P_mi_l, mdl['vd_ins'])
            entropy_tmp = self.get_entropy_event(insertion_event_nickname)
            entropy_tmp = entropy_tmp.values + np.dot(H_P_mi_l, self.Pmarginal[insertion_event_nickname].values)
            list_Insertion_entropy.append(entropy_tmp)

        data_tmp = {
            'event_nickname': list_Insertion_events_nickname,
            'event_type': ['Insertion' for ev_nick in list_Insertion_events_nickname],
            'seq_type': [ self.Pconditionals[ev_nick].attrs['seq_type'] for ev_nick in self.event_Insertion_nickname_list],
            'entropy': list_Insertion_entropy
        }
        # print("list_Insertion_entropy: ", list_Insertion_entropy, type(list_Insertion_entropy[0]))
        df_Insertion_entropy = pd.DataFrame(data_tmp)
        df_Insertion_entropy['entropy'] = df_Insertion_entropy['entropy'].astype(float)

        return df_Insertion_entropy
        # self.Pconditionals.attrs['seq_type']
        # event_nickname_insertion = 'vd_ins'
        # H_P_ins = self.get_entropy_event(event_nickname_insertion)

    def get_df_GeneChoice_entropy_contribution(self):
        """
        Return pandas dataframe entropy decomposition for GeneChoice
        """
        list_GeneChoice_entropy = list()
        for event_GeneChoice_nickname in self.event_GeneChoice_nickname_list:
            list_GeneChoice_entropy.append(self.get_entropy_event(event_GeneChoice_nickname))
        data_tmp = {
            'event_nickname': self.event_GeneChoice_nickname_list,
            'event_type': ['GeneChoice' for ev_nick in self.event_GeneChoice_nickname_list],
            'seq_type': [self.Pconditionals[ev_nick].attrs['seq_type'] for ev_nick in self.event_GeneChoice_nickname_list],
            'entropy': list_GeneChoice_entropy
        }
        df_GeneChoice_entropy = pd.DataFrame(data_tmp)
        df_GeneChoice_entropy['entropy'] = df_GeneChoice_entropy['entropy'].astype(float)
        return df_GeneChoice_entropy

    def get_df_Deletion_entropy_contribution(self):
        """
        Return pandas dataframe entropy decomposition for Deletions
        """
        list_Deletion_entropy = list()
        for event_Deletion_nickname in self.event_Deletion_nickname_list:
            list_Deletion_entropy.append(self.get_entropy_event(event_Deletion_nickname))

        data_tmp = {
            'event_nickname': self.event_Deletion_nickname_list,
            'event_type': ['Deletion' for ev_nick in self.event_Deletion_nickname_list],
            'seq_type': [self.Pconditionals[ev_nick].attrs['seq_type'] for ev_nick in self.event_Deletion_nickname_list],
            'entropy': list_Deletion_entropy
        }
        df_Deletion_entropy = pd.DataFrame(data_tmp)
        df_Deletion_entropy['entropy'] = df_Deletion_entropy['entropy'].astype(float)
        return df_Deletion_entropy


    def get_df_entropy_decomposition(self):
        """
        Return entropy decomposition of events in a Dataframe
        """
        try:
            np.seterr(divide='ignore')
            # calculate entropy of V, D and J Genechoice
            df_GeneChoice_entropy = self.get_df_GeneChoice_entropy_contribution()

            # calculate deletion entropy
            df_Deletion_entropy = self.get_df_Deletion_entropy_contribution()

            # calculate insertion entropy
            df_Insertion_entropy = self.get_df_Insertion_entropy_contribution()

            np.seterr(divide='warn')

            df_entropy_decomposition = pd.concat([df_GeneChoice_entropy, df_Insertion_entropy, df_Deletion_entropy])
            df_entropy_decomposition.reset_index(inplace=True)

            return df_entropy_decomposition.drop(columns='index')
        except Exception as e:
            raise e

    def plot_recombination_entropy(self, ax=None, df_entropy=None):
        if df_entropy is None:
            df_entropy = self.get_df_entropy_decomposition()

        if ax is None:
            # print("matplotlib")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(25, 12))

        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle
        edgecolor = 'g'


        # ax.set_xlim(0, 20)
        # ax.set_ylim(0, 20)
        # ax.add_patch(Rectangle( (0,0), 5, 5 , edgecolor='red') )
        ypos = 0
        left = 0
        height = 0.10 * df_entropy['entropy'].sum()
        for index, row in df_entropy.iterrows():  # ['entropy']:
            h_entropy = row['entropy']
            color = Igor_event_type_color_dict[row['event_type']]
            # print(aaa)
            hatch = Igor_seq_type_hatch_pattern_dict[row['seq_type']]
            ax.barh(ypos, h_entropy, height=height, left=left, color=color, edgecolor=edgecolor, hatch=hatch)
            # label
            seq_type = row['seq_type'].split('_')[0]
            xpos = left + h_entropy / 2
            label = "{} \n{:.1f} \nbits".format(seq_type, h_entropy)
            # label = "{}".format(seq_type)
            ax.text(xpos, ypos - 0.25 * height, label, ha='center', color='k', weight='bold', fontsize='xx-large')
            left += h_entropy

        height = 0.5 * height
        ypos += height
        left = 0
        for event_type in df_entropy[
            'event_type'].unique():  # df_entropy.groupby('event_type').apply(lambda g: g['entropy'].sum()):
            h_event_type = df_entropy[df_entropy['event_type'] == event_type]['entropy'].sum()
            color = Igor_event_type_color_dict[event_type]
            ax.barh(ypos, h_event_type, height=height, left=left, color=color, edgecolor=edgecolor)
            xpos = left + h_event_type / 2
            label = "{}: {:.1f} bits".format(event_type, h_event_type)
            ax.text(xpos, ypos, label, ha='center', color='k', weight='bold', fontsize='xx-large')
            left += h_event_type

        # height = 0.5*height
        ypos += height
        left = 0
        h_recombination = df_entropy['entropy'].sum()
        color = '#f68294'
        xpos = left + h_recombination / 2
        label = "Recombination Events: {:.1f} bits".format(h_recombination)
        ax.text(xpos, ypos, label, ha='center', color='k', weight='bold', fontsize='xx-large')
        ax.barh(ypos, h_recombination, height=height, left=left, color=color, edgecolor=edgecolor)
        ax.set_aspect('equal')
        ax.axis('off')
        return ax

    @staticmethod
    def get_cross_entropy(self):
        pass

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
        """
        Export model events in different csv files for event.
        :param fln_prefix: filename prefix to save events files
        :param sep: csv field separator
        """
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
        # print("********", dependencias, strEvent)
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
        else:
            print("Recombination event " + strEvent + " has an export problem!")

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
            # print(parents)
            evento = self.parms.get_Event(strEvent)
            # print(evento.event_type)
            # print(evento.seq_type)
            dependencias = list(self.xdata[strEvent].dims)
            # print("********", dependencias, strEvent)
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
                        title = "P(" + da["lbl__" + strEvent].values[ii] + "| " + dependencias[0] + "," + dependencias[
                            1] + ")"
                        ofile.write("\n" + title + "\n")
                        da_ii = da[{strEvent: ii}]
                        df = pd.DataFrame(data=da_ii.values, index=da['lbl__' + dependencias[0]].values,
                                          columns=da['lbl__' + dependencias[1]].values)

                        df.to_csv(ofile, mode='a', sep=sep)  # , index_label=evento.seq_type)
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
                        title = "P(" + da["lbl__" + strEvent].values[ii] + "| " + dependencias[0] + "," + dependencias[
                            1] + ")"
                        ofile.write("\n" + title + "\n")
                        da_ii = da[{strEvent: ii}]
                        df = pd.DataFrame(data=da_ii.values, index=da['lbl__' + dependencias[0]].values,
                                          columns=da['lbl__' + dependencias[1]].values)

                        df.to_csv(ofile, mode='a', sep=sep)  # , index_label=evento.seq_type)
            else:
                print("Recombination event " + strEvent + " has an export problem!")

            # return df

            ## P(D3, D5 | D) = P( D3| D5,D) x P (D5,D)
            #### Deletions in D
            da = self.xdata['d_3_del'] * self.xdata['d_5_del']

            ### DELETIONS
            strEvent = 'd_gene'
            da = self.xdata[strEvent]
            dependencias = list(da.dims)
            # print("********", dependencias, strEvent)
            dependencias.remove(strEvent)
            dependencias_dim = [da[dep].shape[0] for dep in dependencias]

            lbl_file = fln_prefix + "P__" + strEvent + "__deletions" + ".csv"
            with open(lbl_file, 'w') as ofile:
                for ii in da[strEvent].values:
                    da_ii = da[{strEvent: ii}]
                    lbl_event_realization = da['lbl__' + strEvent].values[ii]
                    title = "_P(" + dependencias[0] + "," + dependencias[
                        1] + "| " + strEvent + " = " + lbl_event_realization + ")"
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
                                 columns=["P(" + strEvent + ")"])  # da['lbl__' + strEvent].values
            strEvent = 'dj_ins'
            da = self.xdata[strEvent]
            df_dj = pd.DataFrame(data=da.values, index=da['lbl__' + strEvent].values,
                                 columns=["P(" + strEvent + ")"])  # da['lbl__' + strEvent].values

            df = df_vd.merge(df_dj, left_index=True, right_index=True)
            lbl_file = fln_prefix + "P__" + "insertions" + ".csv"
            df.to_csv(lbl_file, index_label="Insertions", sep=sep)

            ### DINUCL
            strEvent = 'vd_dinucl'
            da = self.xdata[strEvent]
            # print(da)
            df = pd.DataFrame(data=da.values, index=da['lbl__x'].values,
                              columns=da['lbl__y'].values)
            lbl_file = fln_prefix + "P__" + strEvent + ".csv"
            df.to_csv(lbl_file, index_label="From\To", sep=sep)

            strEvent = 'dj_dinucl'
            da = self.xdata[strEvent]
            # print(da)
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

    def export_Pmarginal_to_csv(self, event_nickname: str, *args, **kwargs):

        if kwargs.get('sep') is None:
            kwargs['sep'] = ';'

        event = self.parms.get_Event(event_nickname, by_nickname=True)
        da = self.xdata[event_nickname]
        if event.event_type == 'GeneChoice':
            df = da.to_dataframe(name="P")  # .drop('priority', 1)
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
            df.to_csv(*args, **kwargs)  # , index_label=, sep=sep)

        else:
            print("Event nickname " + event_nickname + " is not present in this model.")
            print("Accepted Events nicknames are : " + str(self.get_events_nicknames_list()))

    # FIXME: CHANGE EVENT MARGINAL!!!
    def get_Event_Marginal(self, event_nickname: str):
        """Returns an xarray with the marginal probability of the event given the nickname"""
        # FIXME: add new way to make the recursion.
        # FIXME: FIRST without recursion, return
        if event_nickname in self.parms.get_EventsNickname_list():
            da_event = self.xdata[event_nickname]
            dependencies = self.parms.Edges_dict[event_nickname]
            event_parents = list(self.parms.G.predecessors(event_nickname))
            # 1. Sort the dependencies by priority, then by dependencie
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
                print('&' * 20)
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
                da_marginal = da_event_parent * da_event_parent2
                return da_marginal

        else:
            print("Event nickname : " + event_nickname + " is not an event in this IGoR model.")
            return list()

    def plot_event_GeneChoice(self, event_nickname: str, **kwargs):
        """ Return GeneChoice plot """
        # Default values in plot

        import numpy as np
        v_genLabel = np.vectorize(genLabel)
        import matplotlib.pyplot as plt
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
            fig.tight_layout()
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
            fig.tight_layout()
            return fig, ax
        elif len(parents_list) == 2:
            if not 'cmap' in kwargs.keys():
                kwargs['cmap'] = 'gnuplot2_r'
            # da = self.xdata[event_nickname]
            fig, ax = plt.subplots(*da[event_nickname].shape, figsize=(10, 20))
            for ii, ev_realiz in enumerate(da[event_nickname]):
                # print(ev_realiz.values)
                da[{event_nickname: ev_realiz.values}].plot(ax=ax[ii], cmap='gnuplot2_r')
                lbl_ev_realiz = str(ev_realiz["lbl__" + event_nickname].values)
                lbl_parents = str(",".join(da.attrs['parents']))
                titulo = "$P($" + event_nickname + "$ = $ " + lbl_ev_realiz + " $|$" + lbl_parents + "$)$"
                ax[ii].set_title(titulo)
            fig.tight_layout()
            return fig, ax
        else:
            fig, ax = plt.subplots()
            ax.set_title("Dimensionality not supportted for event : ", event_nickname)
            fig.tight_layout()
            return fig, ax

    def plot_event_Insertion(self, event_nickname: str, **kwargs):
        """ Return Insertion plot """

        import numpy as np
        v_genLabel = np.vectorize(genLabel)
        import matplotlib.pyplot as plt
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
            fig.tight_layout()
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
            fig.tight_layout()
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
            fig.tight_layout()
            return fig, ax
        else:
            fig, ax = plt.subplots()
            ax.set_title("Dimensionality not supportted for event : ", event_nickname)
            fig.tight_layout()
            return fig, ax

    def plot_event_Deletion(self, event_nickname: str, **kwargs):
        """ Return GeneChoice plot """
        # FIXME: I THINK THAT THE BEST WAY SHOULD BE ONLY RETURN AX,
        #  and to save the figure in a pdf use the matplotlib function getcf() get current figure.
        import numpy as np
        v_genLabel = np.vectorize(genLabel)
        import matplotlib.pyplot as plt
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
            fig.tight_layout()
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
            fig.tight_layout()
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
            fig.tight_layout()
            return fig, ax
        else:
            fig, ax = plt.subplots()
            ax.set_title("Dimensionality not supportted for event : ", event_nickname)
            fig.tight_layout()
            return fig, ax

    def plot_event_DinucMarkov(self, event_nickname: str, **kwargs):
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

        fig.tight_layout()
        return fig, ax

    def export_plot_Pconditionals(self, outfilename_prefix):
        """
        Create a pdf file with preliminary plots of conditional probabilities
        :param outfilename_prefix: Prefix for pdf file
        """
        self.export_plot_events(outfilename_prefix)

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

    def plot_Event(self, event_nickname: str, ax=None, **kwargs):
        event = self.parms.get_Event(event_nickname)
        if event.event_type == 'GeneChoice':
            return self.plot_event_GeneChoice(event_nickname)
        elif event.event_type == 'Insertion':
            return self.plot_event_Insertion(event_nickname)
        elif event.event_type == 'Deletion':
            fig, ax = self.plot_event_Deletion(event_nickname)
            fig.tight_layout()
            return fig, ax
        elif event.event_type == 'DinucMarkov':
            fig, ax = self.plot_event_DinucMarkov(event_nickname)
            fig.tight_layout()
            return fig, ax
        else:
            print("ERROR: EVENT NOT RECOGNIZE", event_nickname)

    def export_plot_Pmarginals(self, outfilename_prefix):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        with PdfPages(outfilename_prefix + ".pdf") as pdf_file:
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

    def plot_Event_Marginal(self, event_nickname: str, ax=None, **kwargs):
        """
        Plot marginals of model events by nickname
        """
        event = self.parms.get_Event(event_nickname, by_nickname=True)
        da = self.Pmarginal[event_nickname]  # real marginal DataArray

        lblEvent = event_nickname.replace("_", " ")
        xEtiqueta = lblEvent
        yEtiqueta = "P"

        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()

        ax.set_xlabel(xEtiqueta)
        ax.set_ylabel(yEtiqueta, rotation=0)

        if event.event_type == 'GeneChoice':
            # Bar plot
            XX = da[event_nickname].values
            YY = da.values

            ax.bar(XX, YY, **kwargs)

            lbl_XX = da['lbl__' + event_nickname].values
            ax.set_xticks(XX)
            ax.set_xticklabels(v_genLabel(lbl_XX), rotation=90)
            # return ax

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
            # YY = self.xdata[event_nickname].values
            # XX = self.xdata[event_nickname]['lbl__' + event_nickname].values
            # ax.plot(XX, YY)
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
            print("Event nickname " + event_nickname + " is not present in this model.")
            print("Accepted Events nicknames are : " + str(self.get_events_nicknames_list()))

        # return self.get_Event_Marginal(nickname)
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

    def get_sorted_events_nicknames_list(self):
        return [ event.nickname for event in self.parms.get_Event_list_sorted()]

    # PLOTS:
    def plot_Bayes_network(self, ax=None, filename=None):
        if filename is None:
            if ax is None:
                return self.parms.plot_Graph()
            else:
                ax_ = self.parms.plot_Graph(ax=ax)
                return ax_

        else:
            if ax is None:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots()
                ax_ = self.parms.plot_Graph(ax=ax)
                fig.savefig(filename)
                return ax_
            else:
                ax_ = self.parms.plot_Graph(ax=ax)
                return ax_


    def plot(self, event_nickname: str, ax=None):
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


    def plot_GeneChoice_Pmarginal(self):
        try:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(18, 20))
            grid = plt.GridSpec(2, 3, hspace=0.2, wspace=0.2)
            grid
            ax_v = fig.add_subplot(grid[0, :])
            ax_j = fig.add_subplot(grid[1, 0:2])
            ax_d = fig.add_subplot(grid[1, 2])

            self.plot_Event_Marginal('v_choice', ax=ax_v)
            self.plot_Event_Marginal('j_choice', ax=ax_j)
            self.plot_Event_Marginal('d_gene', ax=ax_d)
            return fig
        except Exception as e:
            raise e

    def plot_InsertionsDeletions_Pmarginal(self):
        try:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 2, figsize=(16, 6))
            self.plot_Event_Marginal('v_3_del', ax=ax[0], marker='o', label='V3')
            self.plot_Event_Marginal('j_5_del', ax=ax[0], marker='x', label='J5')
            self.plot_Event_Marginal('d_3_del', ax=ax[0], marker='s', label='D3')
            self.plot_Event_Marginal('d_5_del', ax=ax[0], marker='s', label='D5')
            ax[0].set_xlabel("Deletions")
            ax[0].legend()

            self.plot_Event_Marginal('vd_ins', ax=ax[1], marker='o', label='VD')
            self.plot_Event_Marginal('dj_ins', ax=ax[1], marker='s', label='DJ')
            ax[1].set_xlabel("Insertions")
            ax[1].legend()
            return fig
        except Exception as e:
            raise e
    # def plot_dumm_report(self, strEvent):
    #     # strEvent = 'd_gene'
    #     import matplotlib.pyplot as plt
    #     fig, ax = plt.subplots()
    #     dependencias = list(self.xdata[strEvent].dims)
    #     dependencias.remove(strEvent)
    #     dependencias_dim = [self.xdata[strEvent][dep].shape[0] for dep in dependencias]
    #     # eventos = eventos.remove(strEvent)
    #     lista = list()
    #     import numpy as np
    #     # np.ndindex()
    #     for index in np.ndindex(*dependencias_dim):
    #         dictionary = dict(zip(dependencias, index))
    #         # TODO: PLOT EACH DAMM FIGURE
    #         self.xdata[strEvent][dictionary].plot()
    #         aaa = [str(key) + "__" + str(dictionary[key]) for key in dictionary.keys()]
    #         lbl_file = "___".join(aaa)
    #         df = self.xdata[strEvent][dictionary].to_dataframe("P").drop('priority', 1)
    #         df.plot.bar(x="lbl__" + strEvent, y='P', ax=ax)
    #         print("*" * 10)
    #         print(lbl_file)
    #         print(df)
    #         # df.to_csv(lbl_file+".csv")
    #         # fig.savefig(lbl_file+".png")
    #         # ax.clear()
    #     return fig

    def set_genomic_dataframe_dict(self, dataframe_dict):
        # TODO: UPDATE THE df_V_ref_genome in parms and generate_xdata()
        try:
            self.genomic_dataframe_dict = dataframe_dict
            self.parms.set_event_realizations_from_DataFrame(self.parms.event_GeneChoice_V.nickname,
                                                             dataframe_dict['V'])
            if self.parms.event_GeneChoice_D is not None:
                self.parms.set_event_realizations_from_DataFrame(self.parms.event_GeneChoice_D.nickname,
                                                                 dataframe_dict['D'])
            self.parms.set_event_realizations_from_DataFrame(self.parms.event_GeneChoice_J.nickname,
                                                             dataframe_dict['J'])

            try:
                if 'anchor_index' in dataframe_dict['V'].columns:
                    self.parms.df_V_anchors = get_df_anchors_from_df_ref_genome(dataframe_dict['V'])
            except Exception as e:
                print(e)
                pass

            try:
                if 'anchor_index' in dataframe_dict['J'].columns:
                    self.parms.df_J_anchors = get_df_anchors_from_df_ref_genome(dataframe_dict['J'])
            except Exception as e:
                print(e)
                pass

            self.generate_xdata()

        except Exception as e:
            raise e


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
            return list(map(lambda x: self.parms.get_Event(event_nickname).realizations[x], event_id))
        else:
            return self.parms.get_Event(event_nickname).realizations[event_id]

    def get_realizations_dict_from_scenario_dict(self, scenario_realization_dict: dict):
        realization_dict = dict()
        # print(scenario_realization_dict)
        for event_nickname, event_id in scenario_realization_dict.items():
            if not (event_nickname == 'mismatcheslen' or event_nickname == 'mismatches' or event_nickname == 'errors'):
                realization_dict[event_nickname] = self.get_event_realization_of_event(event_nickname, event_id)

        return realization_dict

    def set_realization_event_from_DataFrame(self, event_nickname, new_df):
        self.parms.set_event_realizations_from_DataFrame(event_nickname, new_df)
        self.marginals.initialize_uniform_from_model_parms(self.parms)
        self.generate_xdata()

    def write_model(self, fln_model_parms, fln_model_marginals, fln_V_gene_CDR3_anchors=None, fln_J_gene_CDR3_anchors=None):
        """
        Write model parms and marginals(conditional probabilities) in IGoR txt format files.
        :param fln_model_parms: Filename for model parameters.
        :param fln_model_marginals: Filename for model marginals (conditional probabilities).
        :param fln_V_gene_CDR3_anchors: Filename of CDR3 anchors for V gene(optional).
        :param fln_J_gene_CDR3_anchors: Filename of CDR3 anchors for J gene(optional).
        """
        self.parms.Event_list = self.parms.get_Event_list_sorted()
        self.parms.write_model_parms(filename=fln_model_parms)
        self.marginals.write_model_marginals(filename=fln_model_marginals, model_parms=self.parms)
        if fln_V_gene_CDR3_anchors is not None:
            try:
                write_geneanchors_dataframe_to_csv(fln_V_gene_CDR3_anchors, self.genomic_dataframe_dict['V'])
            except Exception as e:
                raise e

        if fln_J_gene_CDR3_anchors is not None:
            try:
                write_geneanchors_dataframe_to_csv(fln_J_gene_CDR3_anchors, self.genomic_dataframe_dict['J'])
            except Exception as e:
                raise e


    def write_ref_genome(self, fln_genomicVs=None, fln_genomicDs=None, fln_genomicJs=None,
                         fln_V_gene_CDR3_anchors=None, fln_J_gene_CDR3_anchors=None):
        """
        Write ref_genome from genomic_dataframe_dict
        :param fln_genomicVs: V fasta file
        :param fln_genomicDs: D fasta file
        :param fln_genomicJs: J fasta file
        :param fln_V_gene_CDR3_anchors: V csv anchors file
        :param fln_J_gene_CDR3_anchors: J csv anchors file
        """
        try:
            if 'V' in self.genomic_dataframe_dict:
                # Write fasta file
                # Write anchors
                try:
                    write_ref_genome_files_from_dataframe(self.genomic_dataframe_dict['V'], fln_genomicVs,
                                                          fln_V_gene_CDR3_anchors)
                except Exception as e:
                    print("ERROR: write_ref_genome ", self.genomic_dataframe_dict)
                    raise e

            if 'J' in self.genomic_dataframe_dict:
                # Write fasta file
                # Write anchors
                try:
                    write_ref_genome_files_from_dataframe(self.genomic_dataframe_dict['J'], fln_genomicJs,
                                                          fln_J_gene_CDR3_anchors)
                except Exception as e:
                    print("ERROR: write_ref_genome ", self.genomic_dataframe_dict)
                    raise e

            if 'D' in self.genomic_dataframe_dict:
                # Write fasta file
                # Write anchors
                try:
                    write_ref_genome_files_from_dataframe(self.genomic_dataframe_dict['D'], fln_genomicDs)
                except Exception as e:
                    print("WARNING: write_ref_genome D genes not found.") #, self.genomic_dataframe_dict)
                    pass

        except Exception as e:
            raise e

    def write_ref_genome_dir(self, ref_genome_dir=None):
        """
        Write genome references files in directory ref_genome_dir

        """
        # TODO: MAKE DIRECTORY
        import os
        os.makedirs(ref_genome_dir)
        ref_genome_dir = ref_genome_dir + "/"
        fln_genomicVs = ref_genome_dir + "genomicVs.fasta"
        fln_genomicDs = ref_genome_dir + "genomicDs.fasta"
        fln_genomicJs = ref_genome_dir + "genomicJs.fasta"
        fln_V_gene_CDR3_anchors = ref_genome_dir + "V_gene_CDR3_anchors.csv"
        fln_J_gene_CDR3_anchors = ref_genome_dir + "J_gene_CDR3_anchors.csv"
        self.write_ref_genome(fln_genomicVs=fln_genomicVs, fln_genomicDs=fln_genomicDs,
                              fln_genomicJs=fln_genomicJs,
                              fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                              fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors)

    def write_mdldata_dir(self, model_dir_path, sep=';', b_igor_directory=True):
        """
        Export IgorModel and IgorRefGenome
        :param model_dir_path: Directory to save model in IGoR directory structure
        :param sep: Field separator, default ';'
        :param b_igor_directory: If True default IGoR directory structure, default False.
        """
        try:
            import pathlib
            pathlib.Path(model_dir_path).mkdir(parents=True, exist_ok=True)
            if b_igor_directory:
                # os.makedirs(model_dir_path, exist_ok=True)
                os.makedirs(model_dir_path + "/models", exist_ok=True)
                fln_dict = get_default_fln_names_for_model_dir(model_dir_path)
                self.write_model(fln_dict['fln_model_parms'], fln_dict['fln_model_marginals'])

                os.makedirs(model_dir_path + "/ref_genome", exist_ok=True)
                fln_dict.pop('fln_model_parms', None)
                fln_dict.pop('fln_model_marginals', None)
                ref_genome = self.parms.get_IgorRefGenome()
                ref_genome.write_ref_genome_dir(model_dir_path+"/ref_genome", sep=sep)
            else:
                # os.makedirs(model_dir_path, exist_ok=True)
                os.makedirs(model_dir_path, exist_ok=True)
                fln_model_parms = model_dir_path + "/" + "model_params.txt"
                fln_model_marginals = model_dir_path + "/" + "model_marginals.txt"
                self.write_model(fln_model_parms, fln_model_marginals)
                ref_genome = self.parms.get_IgorRefGenome()
                ref_genome.write_ref_genome_dir(model_dir_path, sep=sep)
        except Exception as e:
            raise e
        # self.write_ref_genome(**fln_dict)

    def get_nicknames_for_gene_segment(self, strGene):
        """
        Return tuple ('GeneChoice_nickname', 'Five_prime_nickname', 'Three_prime_nickname')
        to be used to construct a scenario sequence with the function get_gene_segment
        :param strGene: V, D or J are only accepted.
        """
        if strGene in ['V', 'D', 'J']:
            str_seq_type = str(strGene).upper() + "_gene"
            list_seq_type_events = [event for event in self.parms.Event_list if event.seq_type == str_seq_type]
            # separate GeneChoice and Deletion
            str_event_type = 'GeneChoice'
            list_event_type_GeneChoice = [event for event in list_seq_type_events if event.event_type == str_event_type]
            event_nickname_GeneChoice = list_event_type_GeneChoice[
                0].nickname  # in principle is only 1, but for future extensions could be 2 (D gene)

            str_event_type = 'Deletion'
            list_event_type_Deletion = [event for event in list_seq_type_events if event.event_type == str_event_type]

            # 5 prime
            str_seq_side = 'Five_prime'
            list_seq_side_Five_prime = [event for event in list_event_type_Deletion if event.seq_side == str_seq_side]
            if len(list_seq_side_Five_prime) > 0:
                event_nickname_Five_prime = list_seq_side_Five_prime[0].nickname
            else:
                event_nickname_Five_prime = None
            # 3 prime
            str_seq_side = 'Three_prime'
            list_seq_side_Three_prime = [event for event in list_event_type_Deletion if event.seq_side == str_seq_side]
            if len(list_seq_side_Three_prime) > 0:
                event_nickname_Three_prime = list_seq_side_Three_prime[0].nickname
            else:
                event_nickname_Three_prime = None

            return event_nickname_GeneChoice, event_nickname_Five_prime, event_nickname_Three_prime

        elif strGene in ['VD']:
            # FIXME: THIS SHOULD BE FIX IN IGOR FOR VD, but is too late to change it!!!!
            str_seq_type = str(strGene).upper() + "_genes"
            list_seq_type_events = [event for event in self.parms.Event_list if event.seq_type == str_seq_type]
            # separate GeneChoice and Deletion
            str_event_type = 'DinucMarkov'
            list_event_type_DinucMarkov = [event for event in list_seq_type_events if event.event_type == str_event_type]
            event_nickname_DinucMarkov = list_event_type_DinucMarkov[0].nickname
            return event_nickname_DinucMarkov
        elif strGene in ['DJ', 'VJ']:
            str_seq_type = str(strGene).upper() + "_gene"
            list_seq_type_events = [event for event in self.parms.Event_list if event.seq_type == str_seq_type]
            # separate GeneChoice and Deletion
            str_event_type = 'DinucMarkov'
            list_event_type_DinucMarkov = [event for event in list_seq_type_events if
                                           event.event_type == str_event_type]
            event_nickname_DinucMarkov = list_event_type_DinucMarkov[0].nickname
            return event_nickname_DinucMarkov

    def get_gene_segment_dict(self, strGene:str, ps_scenario:pd.Series):
        """Return cuted gene or expanded with palidromic insertions for a scenario
        :param strGene: 'V', 'D', 'J', 'VD', 'DJ', or 'VJ'
        :param ps_scenario: scenario as a pandas Series.
        """
        if strGene in ['V', 'D', 'J']:
            tuple_nickname_gene_segment = self.get_nicknames_for_gene_segment(strGene)
            list_realization_value = list()
            str_description = ""
            for ev_nickname in tuple_nickname_gene_segment:
                if ev_nickname is None:
                    list_realization_value.append(None)
                else:
                    ev_realiz = self.realization(ps_scenario, ev_nickname)
                    list_realization_value.append(ev_realiz.value)
                    str_description = str_description + ev_realiz.name + " ("+ ev_nickname+": " + str(ev_realiz.id) + ")"

            str_gene_template = list_realization_value[0]
            int_gene_5_del = list_realization_value[1]
            int_gene_3_del = list_realization_value[2]
            gene_segment_dict = get_gene_segment(str_gene_template, int_gene_5_del=int_gene_5_del,
                                                 int_gene_3_del=int_gene_3_del)

            if not int_gene_5_del is None:
                str_description = str_description + ", 5'del : " + str(int_gene_5_del)

            if not int_gene_3_del is None:
                str_description = str_description + ", 3'del : " + str(int_gene_3_del)
            gene_segment_dict['gene_description'] = str_description

            return gene_segment_dict
        elif strGene in ['VD', 'VJ']:
            nickname_dinucl = self.get_nicknames_for_gene_segment(strGene)
            ev_realiz_dinucl = self.realization(ps_scenario, nickname_dinucl)
            dinucl_segment_dict = collections.OrderedDict()

            dinucl_segment_dict['gene_segment'] = "".join(ev_realiz_dinucl.value)
            dinucl_segment_dict['gene_description'] = strGene +", ins: "+ str(len(ev_realiz_dinucl.value))

            return dinucl_segment_dict
        elif strGene in ['DJ']:
            nickname_dinucl = self.get_nicknames_for_gene_segment(strGene)
            ev_realiz_dinucl = self.realization(ps_scenario, nickname_dinucl)
            dinucl_segment_dict = collections.OrderedDict()

            dinucl_segment_dict['gene_segment'] = "".join(ev_realiz_dinucl.value[::-1])
            dinucl_segment_dict['gene_description'] = strGene +", ins: "+ str(len(ev_realiz_dinucl.value))

            return dinucl_segment_dict
        else:
            return None

    def get_df_scenario_aln_from_scenario(self, ps_scenario):
        """
        Return a Dataframe with the informations need it to show an alignment from a scenario.
        :param ps_scenario: Pandas Series (row from df_scenarios)
        """
        offset = 0
        # TODO: MAKE THIS LIST A GLOBAL VARIBLE IN UTILS WITH A GOOD NAME.
        list_cols_4_alignment = ["segment_description", "gene_description", "gene_template",
                                 "int_gene_5_del", "int_gene_3_del",
                                 "offset", "palindrome_5_end", "gene_ini",
                                 "gene_end", "gene_cut", "palindrome_3_end", "gene_segment"]
        df_scenario_aln = pd.DataFrame(columns=list_cols_4_alignment)
        # FIXME: This is for VDJ, VJ is missing
        VDJ_aln_list = ['V', 'VD', 'D', 'DJ', 'J']
        VJ_aln_list = ['V', 'VJ', 'J']
        if self.parms.event_GeneChoice_D is None:
            choose_aln_list = VJ_aln_list
        else:
            choose_aln_list = VDJ_aln_list
        for ii, strGene in enumerate(choose_aln_list):
            ordered_dicto = self.get_gene_segment_dict(strGene, ps_scenario)
            dicto = dict(ordered_dicto)
            dicto['segment_description'] = strGene
            # dicto['segment_sequence'] = strGene
            dicto['offset'] = offset
            # print(dicto)
            df_scenario_aln.loc[ii] = dicto  # , ignore_index=True)
            offset = offset + len(ordered_dicto['gene_segment'])

        df_scenario_aln.aln_scenario_len = offset

        V_offset = df_scenario_aln.loc[df_scenario_aln['segment_description'] == 'V'].offset.values[0]
        J_offset = df_scenario_aln.loc[df_scenario_aln['segment_description'] == 'J'].offset.values[0]

        try:
            V_anchor = self.V_anchor(ps_scenario[self.parms.event_GeneChoice_V.nickname])
            df_scenario_aln.aln_pos_V_anchor = V_offset + V_anchor
        except Exception as e:
            df_scenario_aln.aln_pos_V_anchor = None
            print(e, "V anchor not found")
            pass
        try:
            J_anchor = self.J_anchor(ps_scenario[self.parms.event_GeneChoice_J.nickname])
            df_scenario_aln.aln_pos_J_anchor = J_offset \
                                               + J_anchor \
                                               - int(df_scenario_aln.loc[ df_scenario_aln['segment_description'] == 'J' ][ 'int_gene_5_del' ]) \
                                               + 3 # to take the  W or F and get the JUNCTION for IMGT notation and CDR3 for us
        except Exception as e:
            df_scenario_aln.aln_pos_J_anchor = None
            print(e, "J anchor not found")
            pass

        # print(df_scenario_aln)
        # print(df_scenario_aln.aln_pos_V_anchor)
        # print(df_scenario_aln.aln_pos_J_anchor)
        # print(df_scenario_aln.aln_scenario_len)
        return df_scenario_aln

    def write_df_scenario_aln_FASTA(self, fln_scenario_fasta, ps_scenario):
        """
        Write a scenario alignment in a fasta file.
        :param fln_scenario_fasta: Filename to write align scenario
        :param ps_scenario: Pandas Series scenario
        """
        with open(fln_scenario_fasta, 'w') as ofile:
            df_scenario_aln = self.get_df_scenario_aln_from_scenario(ps_scenario)
            for index, row in df_scenario_aln.iterrows():
                ngaps_right = df_scenario_aln.aln_scenario_len - (row["offset"]  + len(row["gene_segment"]) )
                fasta_desription = row["gene_description"]
                fasta_sequence = '-'*row["offset"] + row["gene_segment"] + '-'*ngaps_right
                ofile.write(">"+ fasta_desription+"\n")
                ofile.write(fasta_sequence + "\n")

    def get_str_seq_from_ps_scenario(self, ps_scenario:pd.Series):
        """
        Return string sequence from an evaluated scenario without errors
        :param ps_scenario: Pandas Series (row from df_scenarios)
        """
        df_scenario_aln = self.get_df_scenario_aln_from_scenario(ps_scenario)
        return df_scenario_aln['gene_segment'].sum()


    def plot_scenario(self, ps_scenario, nt_lim:Union[None,tuple,list]=None, show_CDR3=True, seq_aligned:Union[None,list,tuple]=None, ax=None):
        """
        Return matplotlib fig, ax
        :param ps_scenario: Pandas Series scenario
        :param nt_lim:Union[None,tuple,list] region limits to show the scenario alignment
        default give boundaries around CDR3, if no anchors in model, show the whole scenario
        :param show_CDR3: Show CDR3 lines default(=True)
        :param seq_aligned: List, tuple or np.array of size 3 (seq_index, str_seq, offset)
        """
        # FIXME: A BETTER WAY TO MANAGE THE SCENARIO SHOULD BE IMPLEMENTED
        #  1. From a list of seq_genes (V, VD, D, VJ, J) associate the events of deletion and insertions
        #  2. For an event like D get the realization of 'd_gene' and its corresponding deletions
        #  3. Now get the array with the markers of the convinient positions using mdl.realization
        #  and observables.
        df_scenario_aln = self.get_df_scenario_aln_from_scenario(ps_scenario)
        if seq_aligned is not None:
            try:
                seq_index = seq_aligned[0]
                str_seq = seq_aligned[1]
                seq_V_offset = seq_aligned[2]
                # V_nickname = self.parms.event_GeneChoice_V.nickname
                # seq_aligned = [seq_index, df_secuencias.loc[seq_index].values[0], 0 - df_offsets.loc[seq_index, ps_scenario[V_nickname]]]
                seq_dict = {'segment_description': 'nt_sequence',
                            'gene_description': "seq_index: " + str(seq_index),
                            'gene_template': np.NaN,
                            'int_gene_5_del': np.NaN,
                            'int_gene_3_del': np.NaN,
                            'offset': int(seq_V_offset),
                            'palindrome_5_end': np.NaN,
                            'gene_ini': np.NaN,
                            'gene_end': np.NaN,
                            'gene_cut': np.NaN,
                            'palindrome_3_end': np.NaN,
                            'gene_segment': str_seq
                            }
                # FIXME: ValueError: could not broadcast input array from shape (281) into shape (0)
                df_scenario_aln.loc[len(df_scenario_aln)] = seq_dict
            except Exception as e:
                raise e
        da_scenario_aln = from_df_scenario_aln_to_da_scenario_aln(df_scenario_aln)

        if ax is None:
            # import matplotlib.pyplot as plt
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(20, 10))
        ax = plot_scenario_from_da_scenario_aln(da_scenario_aln, nt_lim=nt_lim, show_CDR3=show_CDR3, ax=ax)
        return ax

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
                                  key=lambda event: 100 * event.priority - len(
                                      self.xdata[event.nickname].attrs['parents']),
                                  reverse=True)
            sequence_arrengement_dict['VJ_gene'] = VJ_gene_list
            arrengement_list = V_gene_list + VJ_gene_list + J_gene_list
        else:
            D_gene_list = sorted(D_gene_list,
                                 key=lambda event: 100 * event.priority - len(
                                     self.xdata[event.nickname].attrs['parents']),
                                 reverse=True)
            sequence_arrengement_dict['D_gene'] = D_gene_list
            VD_gene_list = [event for event in self.parms.Event_list if event.seq_type == 'VD_genes']
            VD_gene_list = sorted(VD_gene_list,
                                  key=lambda event: 100 * event.priority - len(
                                      self.xdata[event.nickname].attrs['parents']),
                                  reverse=True)
            sequence_arrengement_dict['VD_gene'] = VD_gene_list

            DJ_gene_list = [event for event in self.parms.Event_list if event.seq_type == 'DJ_gene']
            DJ_gene_list = sorted(DJ_gene_list,
                                  key=lambda event: 100 * event.priority - len(
                                      self.xdata[event.nickname].attrs['parents']),
                                  reverse=True)
            sequence_arrengement_dict['DJ_gene'] = VD_gene_list

            arrengement_list = V_gene_list + VD_gene_list + D_gene_list + DJ_gene_list + J_gene_list

        self.sequence_construction_event_list = arrengement_list
        return sequence_arrengement_dict  # arrengement_list

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

    @classmethod
    def make_default_VDJ(cls, df_V_ref_genome, df_D_ref_genome, df_J_ref_genome, lims_deletions=None, lims_insertions=None):
        """Create a default VJ model from V and J genes dataframes
        :param df_V_ref_genome: Pandas Dataframe of Genome reference for V gene with CDR3 anchors
        :param df_D_ref_genome: Pandas Dataframe of Genome reference for D gene
        :param df_J_ref_genome: Pandas Dataframe of Genome reference for J gene with CDR3 anchors
        :param lims_deletions: Tuple with min and maximum value for deletions, e.g. (-4,20). Negative numbers are palidromic insertions
        :param lims_insertions: Tuple with min and maximum value for deletions, e.g. (0,30)
        """
        cls = IgorModel()
        cls.genomic_dataframe_dict = dict()
        cls.genomic_dataframe_dict['V'] = df_V_ref_genome
        cls.genomic_dataframe_dict['D'] = df_D_ref_genome
        cls.genomic_dataframe_dict['J'] = df_J_ref_genome
        cls.parms = IgorModel_Parms.make_default_VDJ(df_V_ref_genome, df_D_ref_genome, df_J_ref_genome,
                                                     lims_deletions=lims_deletions, lims_insertions=lims_insertions)
        cls.marginals = IgorModel_Marginals.make_uniform_from_parms(cls.parms)
        cls.generate_xdata()
        return cls

    @classmethod
    def make_default_VJ(cls, df_V_ref_genome, df_J_ref_genome, lims_deletions=None,
                         lims_insertions=None):
        """Create a default VJ model from V and J genes dataframes
        lims_deletions tuple with min and maximum value for deletions, e.g. (-4,20)
        lims_insertions tuple with min and maximum value for deletions, e.g. (0,30)
        """
        cls = IgorModel()
        cls.genomic_dataframe_dict = dict()
        cls.genomic_dataframe_dict['V'] = df_V_ref_genome
        cls.genomic_dataframe_dict['J'] = df_J_ref_genome
        cls.parms = IgorModel_Parms.make_default_VJ(df_V_ref_genome, df_J_ref_genome,
                                                     lims_deletions=lims_deletions, lims_insertions=lims_insertions)
        cls.marginals = IgorModel_Marginals.make_uniform_from_parms(cls.parms)
        cls.generate_xdata()
        return cls

    @classmethod
    def make_default_from_Dataframe_dict(cls, genomic_dataframe_dict, lims_deletions=None, lims_insertions=None):
        cls = IgorModel()
        try:
            cls.genomic_dataframe_dict = genomic_dataframe_dict
            if 'D' in genomic_dataframe_dict:
                cls = IgorModel.make_default_VDJ(genomic_dataframe_dict['V'],
                                                       genomic_dataframe_dict['D'],
                                                       genomic_dataframe_dict['J'],
                                                       lims_deletions=lims_deletions,
                                                       lims_insertions=lims_insertions)
                # cls.parms = IgorModel_Parms.make_default_VDJ(cls.genomic_dataframe_dict['V'],
                #                                              cls.genomic_dataframe_dict['D'],
                #                                              cls.genomic_dataframe_dict['J'],
                #                                              lims_deletions=lims_deletions,
                #                                              lims_insertions=lims_insertions)
                # cls.parms.attach_anchors()
            else:
                cls = IgorModel.make_default_VJ(cls.genomic_dataframe_dict['V'],
                                                            cls.genomic_dataframe_dict['J'],
                                                            lims_deletions=lims_deletions,
                                                            lims_insertions=lims_insertions)
            return cls
        except Exception as e:
            raise e

    # TODO: MAKE A METHOD TO EXPORT A LINE FROM AN SCENARIO
    #  IN DEV:
    def get_AIRR_from_ps_scenario(self, ps_scenario, v_offset=0):
        """
        Return airr format description from IGoR scenario
        """

        df_events_rearrangement = self.get_df_events_rearrangement()

        # 1. From model get events nickname
        # 2. Make the rearrengement list using pygor AIRR schema.
        # 2.1 GeneChoice Events
        for nickname_GeneChoice in self.event_GeneChoice_nickname_list:
            print(nickname_GeneChoice)
        # 2.2 Deletions Events
        for nickname_Deletion in self.event_Deletion_nickname_list:
            print(nickname_Deletion)
        # 2.3 Insertions Events
        for nickname_Insertion in self.event_Insertion_nickname_list:
            print(nickname_Insertion)

        self.event_DinucMarkov_nickname_list
        # 3. return row for pandas dataframe or dict.

        pass



    def get_AIRR_VDJ_rearragement_dict_from_scenario(self, scenario, str_sequence, v_offset=0, pgen=None, junction=None,
                                                     junction_aa=None):
        # get_AIRR_VDJ_rearragement_dict_from_scenario(scenario, indexed_seq.seq_index, indexed_seq.sequence)
        # airr_dict = dict()

        from .AIRR import AIRR_VDJ_rearrangement

        realizations_ids_dict = scenario.realizations_ids_dict
        realization_dict = self.get_realizations_dict_from_scenario_dict(realizations_ids_dict)

        v_segment, vd_segment, d_segment, dj_segment, j_segment = self.construct_sequence_VDJ_from_realization_dict(
            realization_dict)

        airr_vdj = AIRR_VDJ_rearrangement(sequence_id=scenario.seq_index, sequence=str_sequence)

        airr_vdj.v_data.call = realization_dict['v_choice'].name
        airr_vdj.d_data.call = realization_dict['d_gene'].name
        airr_vdj.j_data.call = realization_dict['j_choice'].name

        airr_vdj.sequence_alignment = v_segment['gene_segment'] + vd_segment['gene_segment'] + d_segment[
            'gene_segment'] + dj_segment['gene_segment'] + j_segment['gene_segment']

        airr_vdj.np1 = v_segment['palindrome_3_end'] + vd_segment['gene_segment']
        airr_vdj.np2 = dj_segment['gene_segment'] + j_segment['palindrome_5_end']

        airr_vdj.pgen = pgen

        airr_vdj.junction = junction
        airr_vdj.junction_aa = junction_aa

        airr_vdj.rev_comp = False

        # FIXME: CORRECT CIGAR FORMAT TEMPORARY SOLUTION
        airr_vdj.v_data.cigar = str(len(v_segment['gene_cut'])) + "M"
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
        airr_vdj.np1 = vd_segment['gene_segment']  # This include the palindromic insertions

        # D
        airr_vdj.p5d_length = len(d_segment['palindrome_5_end'])
        airr_vdj.d_data.germline_start = d_segment['gene_ini'] + 1
        airr_vdj.d_data.germline_end = d_segment['gene_end'] + 1
        airr_vdj.d_data.sequence_start = airr_vdj.np1_length + (
                airr_vdj.v_data.sequence_end - airr_vdj.v_data.sequence_start - 1)
        airr_vdj.d_data.sequence_end = airr_vdj.d_data.sequence_start + len(d_segment['gene_cut']) - 1
        airr_vdj.p3d_length = len(d_segment['palindrome_3_end'])

        # DJ
        airr_vdj.n2_length = realization_dict['dj_ins'].value
        airr_vdj.np2_length = airr_vdj.p5d_length + airr_vdj.n2_length + airr_vdj.p3d_length
        airr_vdj.np2 = dj_segment['gene_segment']  # This include the palindromic insertions

        # J
        airr_vdj.p5j_length = len(j_segment['palindrome_5_end'])
        airr_vdj.j_data.germline_start = j_segment['gene_ini'] + 1
        airr_vdj.j_data.germline_end = j_segment['gene_end'] + 1
        airr_vdj.j_data.sequence_start = airr_vdj.np2_length + (
                airr_vdj.d_data.sequence_end - airr_vdj.d_data.sequence_start - 1)
        airr_vdj.j_data.sequence_end = airr_vdj.j_data.sequence_start + len(j_segment['gene_cut']) - 1

        return airr_vdj.to_dict()

    def get_AIRR_VJ_rearragement_dict_from_scenario(self, scenario, str_sequence, v_offset=0, pgen=None, junction=None,
                                                    junction_aa=None):
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
        airr_vj.d_data.call = None  # realization_dict['d_gene'].name
        airr_vj.j_data.call = realization_dict['j_choice'].name

        airr_vj.sequence_alignment = v_segment['gene_segment'] + vj_segment['gene_segment'] + j_segment['gene_segment']

        airr_vj.np1 = v_segment['palindrome_3_end'] + vj_segment['gene_segment'] + j_segment['palindrome_5_end']
        airr_vj.np2 = None

        airr_vj.pgen = pgen

        airr_vj.junction = junction
        airr_vj.junction_aa = junction_aa

        airr_vj.rev_comp = False

        # FIXME: CORRECT CIGAR FORMAT TEMPORARY SOLUTION
        airr_vj.v_data.cigar = str(len(v_segment['gene_cut'])) + "M"
        airr_vj.d_data.cigar = None  # str(len(d_segment['gene_cut'])) + "M"
        airr_vj.j_data.cigar = str(len(j_segment['gene_cut'])) + "M"

        airr_vj.v_data.score = 5 * len(v_segment['gene_cut'])
        airr_vj.d_data.score = None  # 5 * len(d_segment['gene_cut'])
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

    def get_dataframe_from_fln_generated_realizations_werr(self, igor_fln_generated_realizations_werr, sep=';'):
        try:
            print("igor_fln_generated_realizations_werr: ", igor_fln_generated_realizations_werr)
            df2 = pd.read_csv(igor_fln_generated_realizations_werr, sep=';').set_index('seq_index')
            # df2.sort_values(by=['seq_index', 'scenario_rank'], ascending=[True, True])
            events_name__nickname_dict = self.parms.get_event_dict('name', 'nickname')
            events_nickname__event_type_dict = self.parms.get_event_dict('nickname', 'event_type')
            events_nickname__seq_type_dict = self.parms.get_event_dict('nickname', 'seq_type')
            df2.rename(columns=events_name__nickname_dict, inplace=True)

            for column_name in df2.columns:
                try:
                    if column_name in events_nickname__event_type_dict.keys():
                        if events_nickname__event_type_dict[column_name] == 'GeneChoice':
                            # remove parenthesis and make it an int column
                            df2[column_name] = df2[column_name].apply(lambda x: int(eval(x)))
                            seq_type = events_nickname__seq_type_dict[column_name]
                            str_gene_type = seq_type[0].lower()
                            # gene_call_column_name = (str_gene_type+"_call")
                            # df2[gene_call_column_name] = df2[column_name].apply(lambda x: self.parms[column_name].name.loc[x])
                            df2[column_name].apply(lambda x: self.parms[column_name].name.loc[x])

                        elif events_nickname__event_type_dict[column_name] == 'Insertion':
                            # Change to insertions values
                            df2[column_name] = df2[column_name].apply(lambda x: int(eval(x)))
                        elif events_nickname__event_type_dict[column_name] == 'Deletion':
                            # Change to deletions values
                            df2[column_name] = df2[column_name].apply(lambda x: int(eval(x)))
                        elif events_nickname__event_type_dict[column_name] == 'DinucMarkov':
                            # Change to deletions values
                            df2[column_name] = df2[column_name].apply(lambda x: eval(x.replace("(", "[").replace(")", "]")))
                        else:
                            print("IgorModel.get_dataframe_from_fln_generated_realizations_werr: column not found!")
                    else:
                        if (column_name == 'Errors') or (column_name == 'Mismatches'):
                            df2[column_name] = df2[column_name].apply(lambda x: eval(x.replace("(", "[").replace(")", "]")))

                except Exception as e:
                    pass

            return df2
        except Exception as e:
            raise e

    def get_dataframe_scenarios(self, fln_scenarios):
        """
        Return dataframe scenarios from fln_scenarios file.
        :param fln_scenarios: filename of IGoR scenarios file.
        """
        df_scenarios = self.get_dataframe_from_fln_generated_realizations_werr(fln_scenarios)
        self.get_df_normalize_prob(df_scenarios)
        return df_scenarios



    @staticmethod
    def get_df_normalize_prob(df_scenarios):
        df_scenarios['norm_scenario_proba_cond_seq'] = get_df_normalize_prob(df_scenarios)

    def get_probability_matrix_from_event_list_and_scenarios_dataframe(self, event_list:list, df_scenarios:pd.DataFrame):
        """Get probability xarray tensor from IGoR's scenarios dataframe for a given list of events
        :param event_list: Nickname's events list.
        :param df_scenarios: Dataframe with nicknames as headers.
        :return: xarray of joint probability for event_list calculated from
        the weighted ocurrencies in df_scenarios.
        """
        # Initialize xarray tensor
        da_events_prob = self.get_zero_xarray_from_list(event_list)
        df_scenarios['norm_scenario_proba_cond_seq'] = get_df_normalize_prob(df_scenarios)
        aaa = df_scenarios.groupby(event_list)['norm_scenario_proba_cond_seq'].apply(lambda x: x.sum())
        for iii, value in aaa.iteritems():
            coordenadas = dict(zip(aaa.index.names, iii))
            # print(coordenadas, value)  # da_vj_zero[coordenadas])
            da_events_prob[coordenadas] = value

        return da_events_prob

    def get_IgorEvent_realization(self, ps_scenario, event_nickname: Union[None, str, list]=None):
        if event_nickname is None:
            try:
                realization_dict = dict()
                for nickname in self.get_events_nicknames_list():
                    realization_dict[nickname] = self.parms.get_Event_realization(event_nickname, ps_scenario[event_nickname])
                return realization_dict
            except Exception as e:
                raise e
        else:
            try:
                return self.parms.get_Event_realization(event_nickname, ps_scenario[event_nickname])
            except Exception as e:
                raise e


    def get_IgorEvent_realization_for_nickname(self, ps_scenario, event_nickname:str):
        try:
            return self.parms.get_Event_realization(event_nickname, ps_scenario[event_nickname])
        except Exception as e:
            raise e

    def realization(self, ps_scenario, event_nickname:str)->IgorEvent_realization:
        """
        Return realization of scenario.
        """
        try:
            id = ps_scenario[event_nickname]
            # FIXME: isinstance(id, pd.Series)
            realiz = self.parms.Event_dict[event_nickname].loc[id]

            if isinstance(id, list):
                return IgorEvent_realization.from_tuple(id, realiz.value.values, realiz.name.values)

            if isinstance(realiz, pd.DataFrame):
                return IgorEvent_realization.from_tuple(id.values, realiz.value.values, realiz.name.values)
            else:
                return IgorEvent_realization.from_tuple(id, realiz['value'], realiz['name'])
            # return IgorEvent_realization.from_tuple(np.array(id), np.array(aver.value), np.array(aver.name))
        except Exception as e:
            raise e
        # return self.parms.Event_dict['v_choice'].loc[id]
        # return self.get_IgorEvent_realization_for_nickname(ps_scenario, event_nickname)


    def get_df_realizations_dinucl(self, df_scenarios, event_nickname):
        """Return a new dataframe of the Dinucl Markov event with columns id, value and name"""
        id = df_scenarios[event_nickname]
        if id.dtype == 'object':
            # iterate over
            df_realizations = pd.DataFrame(id)
            df_realizations = df_realizations.rename(columns={event_nickname: 'id'})
            df_realizations['value'] = df_realizations['id'].apply(
                lambda x: self.parms.Event_dict[event_nickname].loc[x].value.values)
            df_realizations['name'] = df_realizations['id'].apply(
                lambda x: self.parms.Event_dict[event_nickname].loc[x].name.values)

        return df_realizations

    def realizations_dict(self, ps_scenario):
        """
        Return Ordered dictionary of realization of scenario.
        """
        try:
            # realizations_dict = collections.OrderedDict()
            # for event_nickname in self.parms.Event_dict.keys():
            #     realizations_dict[event_nickname] = self.realization(ps_scenario, event_nickname)
            # return realizations_dict

            realizations_dict = collections.OrderedDict()
            for event_nickname in [x.nickname for x in self.parms.get_Event_list_sorted()]:
                realizations_dict[event_nickname] = self.realization(ps_scenario, event_nickname)
            return realizations_dict

            # return IgorEvent_realization.from_tuple(np.array(id), np.array(aver.value), np.array(aver.name))
        except Exception as e:
            raise e


    def get_df_realizations(self, df_scenarios, event_nickname:Union[None, str]=None):
        """
        Return a dataframe column with the id, value and name column of the realization
        """

        event_realization = self.realization(df_scenarios, event_nickname)
        # TODO: TO A DATAFRAME
        data = {'id': event_realization.id,
                'value': event_realization.value,
                'name': event_realization.name}

        return pd.DataFrame(data, index=df_scenarios.index)


    def w_average_function_df_scenarios(self, observable_func, df_scenarios:pd.DataFrame):
        """Return average of function weigthed with the probability scenarios"""
        average_value = (df_scenarios['norm_scenario_proba_cond_seq'] * self.get_observable_from_df_scenarios(observable_func, df_scenarios)).sum()
        return average_value

    def w_mean_df_scenarios(self, column_name:str, df_scenarios:pd.DataFrame):
        """Return weighted mean with the normalized probabilities for each
        scenario (norm_scenario_proba_cond_seq)
        :param column_name: column name of df_scenario to calculate the average
        :param df_scenarios: Scenarios with normalize probability. Loaded with self.get_dataframe_scenarios()
        """
        # group_column_name = df_scenarios.groupby(column_name)['norm_scenario_proba_cond_seq'].apply(
        #     lambda x: x.sum())
        return (df_scenarios['norm_scenario_proba_cond_seq'] * df_scenarios[column_name]).sum()

    def w_variance_df_scenarios(self, colname_1:str, df_scenarios:pd.DataFrame):
        """Return weighted covariance with the normalized probabilities of the column names given for each
        scenario (norm_scenario_proba_cond_seq)
        :param colname_1: column name of df_scenario to calculate the weighted covariance
        :param colname_2: column name of df_scenario to calculate the weighted covariance
        :param df_scenarios: Scenarios with normalize probability. Loaded with self.get_dataframe_scenarios()
        """
        w_mean_colname_1 = self.w_mean_df_scenarios(colname_1, df_scenarios)

        # FIXME: FINISH THIS COV = w_mean_df_scenarios(
        #  self.w_average_function_df_scenarios( lambda product_x_y, df_scenarios)
        # self.w_average_function_df_scenarios()
        return (df_scenarios['norm_scenario_proba_cond_seq'] * ((df_scenarios[colname_1] - w_mean_colname_1)**2)).sum()

    def w_covariance_df_scenarios(self, colname_1:str, colname_2:str, df_scenarios:pd.DataFrame):
        """Return weighted covariance with the normalized probabilities of the column names given for each
        scenario (norm_scenario_proba_cond_seq)
        :param colname_1: column name of df_scenario to calculate the weighted covariance
        :param colname_2: column name of df_scenario to calculate the weighted covariance
        :param df_scenarios: Scenarios with normalize probability. Loaded with self.get_dataframe_scenarios()
        """
        w_mean_colname_1 = self.w_mean_df_scenarios(colname_1, df_scenarios)
        w_mean_colname_2 = self.w_mean_df_scenarios(colname_2, df_scenarios)
        # FIXME: FINISH THIS COV = w_mean_df_scenarios(
        #  self.w_average_function_df_scenarios( lambda product_x_y, df_scenarios)
        # self.w_average_function_df_scenarios()
        return (df_scenarios['norm_scenario_proba_cond_seq'] * ((df_scenarios[colname_1] - w_mean_colname_1)*(df_scenarios[colname_2] - w_mean_colname_2))).sum()


    def get_P_marginal_from_df_scenarios_cols(self, df_scenarios, colname_list):
        """Get marginalize probabilities of df_scenarios"""
        return df_scenarios.groupby(colname_list)['norm_scenario_proba_cond_seq'].apply(
            lambda x: x.sum())

    def get_P_from_scenarios_cols(self, df_scenarios, colname_list):
        """
        Return xarray with marginalize probabilities of listed columns in dataframe scenarios df_scenarios
        :param df_scenarios: Scenarios with normalize probability. Loaded with self.get_dataframe_scenarios()
        :param colname_list: List of variables preserve for marginalization
        """
        df_marginal = self.get_P_marginal_from_df_scenarios_cols(df_scenarios, colname_list)
        # groupby_colname = df_scenarios.groupby(colname_list)['norm_scenario_proba_cond_seq'].apply(
        #     lambda x: x.sum())
        da = df_marginal.to_xarray()
        da.values = np.nan_to_num(da.values, 0)
        return da


    def get_CDR3_nt_pos_from_ps_scenario(self, ps_scenario, V_offset=None):
        """
        Return tuple of anchors positions V_anchor_in_seq, J_anchor_in_seq
        :param ps_scenario: Pandas Series of one scenario.
        :param V_offset: Offset from IGoR's V alignments.
        """

        try:
            if V_offset is None:
                V_offset = 0
            V_offset = int(V_offset)
            len_GeneChoice = 0
            len_Insertion = 0
            len_Deletion = 0

            # print("ps_scenario: ", ps_scenario)
            for event_nickname in self.event_GeneChoice_nickname_list:
                len_GeneChoice += len(self.realization(ps_scenario, event_nickname).value)
                # print(event_nickname, len_GeneChoice)
            for event_nickname in self.event_Insertion_nickname_list:
                len_Insertion += self.realization(ps_scenario, event_nickname).value
            for event_nickname in self.event_Deletion_nickname_list:
                len_Deletion += -self.realization(ps_scenario, event_nickname).value

            len_scenario = len_GeneChoice + len_Insertion + len_Deletion
            # print('*'*80)
            # print("len_scenario: ", len_scenario)

            V_nickname = self.event_GeneChoice_V_nickname
            J_nickname = self.event_GeneChoice_J_nickname

            V_choice_realization = self.realization(ps_scenario, V_nickname)
            J_choice_realization = self.realization(ps_scenario, J_nickname)

            V_anchor_in_seq = V_offset + int( self.V_anchor(V_choice_realization.id) )
            J_anchor_in_seq = V_offset + int( self.J_anchor(J_choice_realization.id) ) - len(J_choice_realization.value) + len_scenario + 3

            return V_anchor_in_seq, J_anchor_in_seq
        except Exception as e:
            raise e

    def get_VDJ_CDR3_from_scenario(self, ps_scenario):
        """
        Return the numbers of amino acids in vd insertions
        """
        try:
            v_choice = self.realization(ps_scenario, 'v_choice')
            j_choice = self.realization(ps_scenario, 'j_choice')
            d_gene = self.realization(ps_scenario, 'd_gene')

            v_3_del = self.realization(ps_scenario, 'v_3_del')
            d_5_del = self.realization(ps_scenario, 'd_5_del')
            d_3_del = self.realization(ps_scenario, 'd_3_del')
            j_5_del = self.realization(ps_scenario, 'j_5_del')

            # vd_ins = mdl.realization(ps_scenario, 'vd_ins')
            vd_dinucl = self.realization(ps_scenario, 'vd_dinucl')

            # dj_ins = mdl.realization(ps_scenario, 'dj_ins')
            dj_dinucl = self.realization(ps_scenario, 'dj_dinucl')

            v_anchor = self.V_anchor(v_choice.id)
            j_anchor = self.J_anchor(j_choice.id)

            # TODO: mdl.get_CDR3_seq(ps_scenario)
            ##### V_Gene
            v_gene_len = len(v_choice.value)
            # mdl.get_CDR3_seq(ps_scenario)
            v_ini = 0
            v_end = v_gene_len
            str_v_3_palidrome = ""
            if v_3_del.value < 0:
                str_v_3_palidrome = dna_complementary((v_choice.value[v_3_del.value:])[::-1])
            else:
                v_end = v_end - v_3_del.value

            str_V_segment = v_choice.value[v_ini:v_end] + str_v_3_palidrome

            ##### D_gene
            d_gene_len = len(d_gene.value)
            d_ini = 0
            d_end = d_gene_len
            str_d_5_palidrome = ""
            if d_5_del.value < 0:
                int_ini = 0
                str_d_5_palidrome = dna_complementary((d_gene.value[:-d_5_del.value])[::-1])
            else:
                d_ini = d_5_del.value

            str_d_3_palidrome = ""
            if d_3_del.value < 0:
                str_d_3_palidrome = dna_complementary((d_gene.value[d_3_del.value:])[::-1])
            else:
                d_end = d_end - d_3_del.value

            str_D_segment = str_d_5_palidrome + d_gene.value[d_ini:d_end] + str_d_3_palidrome

            ##### J_gene
            j_gene_len = len(j_choice.value)
            j_ini = 0
            j_end = j_gene_len
            str_j_5_palindrome = ""
            if j_5_del.value < 0:
                j_ini = 0
                str_j_5_palindrome = dna_complementary((j_choice.value[:-j_5_del.value])[::-1])
            else:
                j_ini = j_5_del.value

            str_J_segment = str_j_5_palindrome + j_choice.value[j_ini:j_end]

            str_VD_segment = "".join(vd_dinucl.value)
            str_DJ_segment = "".join(dj_dinucl.value[::-1])

            if (v_anchor > v_end) or (j_anchor < j_ini):
                return np.NaN
            else:
                str_sequence = str_V_segment[
                               v_anchor:] + str_VD_segment + str_D_segment + str_DJ_segment + str_J_segment[:j_anchor]
                if len(str_sequence) % 3 == 0:
                    return dna_translate(str_sequence)
                else:
                    return np.NaN
        except Exception as e:
            return None

    def get_VDJ_CDR3_from_df_scenario(self, df_scenario):
        return self.get_observable_from_df_scenarios(self.get_VDJ_CDR3_from_scenario, df_scenario)

    def get_mutual_information_events_from_df_scenarios(self, df_scenarios, event_nickname_x, event_nickname_y):
        """
        Return mutual information in log10 of the desired events
        """
        P_x_y = self.get_P_from_scenarios_cols(df_scenarios, [event_nickname_x, event_nickname_y])
        P_x = self.get_P_from_scenarios_cols(df_scenarios, [event_nickname_x])
        P_y = self.get_P_from_scenarios_cols(df_scenarios, [event_nickname_y])
        I_X_Y = get_D_KL_from_xarray(P_x_y, P_x, P_y)
        return I_X_Y

    def get_mutual_information_from_df_scenarios(self, df_scenarios):
        """
        Return an xarray with the information the mutual information calculated from scenarios dataframe
        :param df_scenarios: Scenarios with normalize probability. Loaded with self.get_dataframe_scenarios()
        """
        try:
            #self.get_sorted_events_nicknames_list()
            dict_nickname_event_type = self.parms.get_event_dict('nickname', 'event_type')
            dict_events = {key: val for key, val in dict_nickname_event_type.items() if val != 'DinucMarkov'}
            event_lista_nicknames = list(dict_events.keys())
            data_0 = np.zeros((len(event_lista_nicknames), len(event_lista_nicknames)))
            da_mi = xr.DataArray(data_0, dims=('x', 'y'), coords={'x': event_lista_nicknames, 'y': event_lista_nicknames})
            da_mi.name = 'mutual_information'

            import itertools
            # for event_nickname_x, event_nickname_y in itertools.product(event_lista_nicknames, event_lista_nicknames):
            # mutual information I(X, Y) = I(Y, X)
            for event_nickname_x, event_nickname_y in itertools.combinations_with_replacement(event_lista_nicknames, 2):
                if event_nickname_x != event_nickname_y:
                    mi = self.get_mutual_information_events_from_df_scenarios(df_scenarios, event_nickname_x,
                                                                             event_nickname_y)
                    da_mi.loc[{"x": event_nickname_x, "y": event_nickname_y}] = mi
                    da_mi.loc[{"x": event_nickname_y, "y": event_nickname_x}] = mi
                else:
                    da_mi.loc[{"x": event_nickname_x, "y": event_nickname_y}] = 0.0
                    # print(event_nickname_x, event_nickname_y, mi)
            return da_mi
        except Exception as e:
            raise e
        # P_x_y = self.get_P_from_scenarios_cols(df_scenarios, [event_nickname_x, event_nickname_y])
        # P_x = self.get_P_from_scenarios_cols(df_scenarios, [event_nickname_x])
        # P_y = self.get_P_from_scenarios_cols(df_scenarios, [event_nickname_y])

    @staticmethod
    def plot_mutual_information(da_mi, ax=None, **kwargs):
        try:

            da_mi_xticks = np.arange(len(da_mi.x))
            da_mi_yticks = np.arange(len(da_mi.y))
            if ax is None:
                import matplotlib.pyplot as plt
                fig, ax = plt.subplots(figsize=(10, 10))

            da_mi.assign_coords(x=da_mi_xticks, y=da_mi_yticks).plot(ax=ax, cmap='gnuplot2_r', **kwargs)

            ax.set_xticks(da_mi_xticks)
            ax.set_yticks(da_mi_yticks)
            ax.set_xticklabels(da_mi.x.values)
            ax.set_yticklabels(da_mi.y.values)
            ax.set_xlabel(None)
            ax.set_ylabel(None)
            ax.set_aspect('equal')
            return ax
        except Exception as e:
            raise e

    def add_Edge(self, parent_nickname:str, child_nickname:str):
        """
        Add an Edge on Bayes network from parent to child, a new dimension is added
        to the conditional probabilities.
        :param parent_nickname: Nickname of parent event.
        :param child_nickname: Nickname of child event.
        """
        self.parms.add_Edge(parent_nickname, child_nickname)
        # MODIFY THE MARGINALS
        import copy
        original_marginal = copy.deepcopy(self.marginals[child_nickname])

        import xarray as xr
        da_parent_ones = xr.ones_like(self.xdata[parent_nickname][parent_nickname]) # FIXME don't use this! SHOULD BE ONE DIMENSIONAL IN GENERAL!!!!!!
        da_child = self.xdata[child_nickname]
        da_new_child = da_parent_ones * da_child
        self.marginals.marginals_dict[child_nickname] = da_new_child.values
        # TODO: update marginals.network_dict

        #for key, value in parms.Edges_dict.items():
        self.marginals.network_dict[parent_nickname] = self.parms.Edges_dict[parent_nickname] + [parent_nickname]
        self.marginals.network_dict[child_nickname] = self.parms.Edges_dict[child_nickname] + [child_nickname]
        self.generate_xdata()

    def remove_Edge(self, parent_nickname:str, child_nickname:str):
        self.parms.remove_Edge(parent_nickname, child_nickname)
        da_child = self.xdata[child_nickname]
        da_child.sum()







class IgorScenario:
    def __init__(self, seq_index:Union[None, int]=None,
                 scenario_rank:Union[None, int]=None,
                 scenario_proba_cond_seq:Union[None, int]=None,
                 realizations_ids_dict:Union[None, dict]=None):
        if seq_index is None:
            self.seq_index = -1
        else:
            self.seq_index = seq_index

        if scenario_rank is None:
            self.scenario_rank = -1
        else:
            self.scenario_rank = scenario_rank

        if scenario_proba_cond_seq is None:
            self.scenario_proba_cond_seq = -1
        else:
            self.scenario_proba_cond_seq = scenario_proba_cond_seq

        if realizations_ids_dict is None:
            self.realizations_ids_dict = dict()
        else:
            self.realizations_ids_dict = realizations_ids_dict

        # self.events_ordered_list = list()

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
    def get_scenario_fasta(self, mdl: IgorModel):
        str_fasta = ""
        # sort events to construct fasta sequence:
        mdl.parms.Event_list
        for key in mdl.xdata.keys():
            self.realizations_ids_dict[key]

        return str_fasta

    def set_model(self, mdl: IgorModel):
        """ Initiate scenario dictionary with a IgorModel """
        for key in mdl.xdata.keys():
            self.realizations_ids_dict[key] = -1

    # TODO: in DEV - FINISH THIS METHOD
    def set_model_from_headers(self, header_line: str):
        # seq_index;scenario_rank;scenario_proba_cond_seq;GeneChoice_V_gene_Undefined_side_prio7_size35;GeneChoice_J_gene_Undefined_side_prio7_size14;GeneChoice_D_gene_Undefined_side_prio6_size2;Deletion_V_gene_Three_prime_prio5_size21;Deletion_D_gene_Five_prime_prio5_size21;Deletion_D_gene_Three_prime_prio5_size21;Deletion_J_gene_Five_prime_prio5_size23;Insertion_VD_genes_Undefined_side_prio4_size31;DinucMarkov_VD_genes_Undefined_side_prio3_size16;Insertion_DJ_gene_Undefined_side_prio2_size31;DinucMarkov_DJ_gene_Undefined_side_prio1_size16;Mismatches
        header_line = "seq_index;scenario_rank;scenario_proba_cond_seq;GeneChoice_V_gene_Undefined_side_prio7_size35;GeneChoice_J_gene_Undefined_side_prio7_size14;GeneChoice_D_gene_Undefined_side_prio6_size2;Deletion_V_gene_Three_prime_prio5_size21;Deletion_D_gene_Five_prime_prio5_size21;Deletion_D_gene_Three_prime_prio5_size21;Deletion_J_gene_Five_prime_prio5_size23;Insertion_VD_genes_Undefined_side_prio4_size31;DinucMarkov_VD_genes_Undefined_side_prio3_size16;Insertion_DJ_gene_Undefined_side_prio2_size31;DinucMarkov_DJ_gene_Undefined_side_prio1_size16;Mismatches"
        header_fields = header_line.split(";")
        events_list = header_fields[3:]

    # FIXME:
    @classmethod
    def load_FromLineBestScenario(cls, line, delimiter=";"):
        # seq_index;scenario_rank;scenario_proba_cond_seq;GeneChoice_V_gene_Undefined_side_prio7_size35;GeneChoice_J_gene_Undefined_side_prio7_size14;GeneChoice_D_gene_Undefined_side_prio6_size2;Deletion_V_gene_Three_prime_prio5_size21;Deletion_D_gene_Five_prime_prio5_size21;Deletion_D_gene_Three_prime_prio5_size21;Deletion_J_gene_Five_prime_prio5_size23;Insertion_VD_genes_Undefined_side_prio4_size31;DinucMarkov_VD_genes_Undefined_side_prio3_size16;Insertion_DJ_gene_Undefined_side_prio2_size31;DinucMarkov_DJ_gene_Undefined_side_prio1_size16;Mismatches
        cls = IgorScenario()
        linesplit = line.split(delimiter)
        for ii in range(len(linesplit)):
            # TODO: find a better way to do this, if is a list keep it as list
            if (ii in [11, 13, 14]):
                linesplit[ii] = linesplit[ii]
            else:
                linesplit[ii] = linesplit[ii].replace("(", "").replace(")", "")

    @classmethod
    def load_FromSQLRecord(cls, sqlRecordScenario: list, sql_scenario_name_type_list: list):
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

    @classmethod
    def load_from_dict(self, dicto):
        dicto_copy = dicto.copy()
        cls = IgorScenario()
        if 'seq_index' in dicto_copy:
            cls.seq_index = dicto_copy['seq_index']
            dicto_copy.pop('seq_index')
        if 'scenario_rank' in dicto_copy:
            cls.scenario_rank = dicto_copy['scenario_rank']
            dicto_copy.pop('scenario_rank')
        if 'scenario_proba_cond_seq' in dicto_copy:
            cls.scenario_proba_cond_seq = dicto_copy['scenario_proba_cond_seq']
            dicto_copy.pop('scenario_proba_cond_seq')
        cls.realizations_ids_dict = dicto_copy
        return cls


    def export_to_AIRR_line(self, scenario_col_list: list, sep='\t'):
        str_line = ""
        self.seq_index = -1
        self.scenario_rank = -1
        self.scenario_proba_cond_seq = -1

        # n_d_5_del = self.mdlParms.Event_dict[strEv].loc[self.id_d_5_del]['value']
        # name_D = self.mdlParms.Event_dict[strEv].loc[self.id_d_gene]['name']
        # header_list=['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'v_score', 'd_score', 'j_score'])
        # sequence_id	sequence	rev_comp	productive	v_call	d_call	j_call	c_call	sequence_alignment	germline_alignment	junction	junction_aa	v_score	v_cigar	d_score	d_cigar	j_score	j_cigar	c_score	c_cigar	vj_in_frame	stop_codon	v_identity	v_evalue	d_identity	d_evalue	j_identity	j_evalue	v_sequence_start	v_sequence_end	v_germline_start	v_germline_end	d_sequence_start	d_sequence_end	d_germline_start	d_germline_end	j_sequence_start	j_sequence_end	j_germline_start	j_germline_end	junction_length	np1_length	np2_length	duplicate_count	consensus_count
        airr_header_list = ["sequence_id", "sequence", "rev_comp", "productive", "v_call", "d_call", "j_call", "c_call",
                            "sequence_alignment", "germline_alignment", "junction", "junction_aa", "v_score", "v_cigar",
                            "d_score", "d_cigar", "j_score", "j_cigar", "c_score", "c_cigar", "vj_in_frame",
                            "stop_codon", "v_identity", "v_evalue", "d_identity", "d_evalue", "j_identity", "j_evalue",
                            "v_sequence_start", "v_sequence_end", "v_germline_start", "v_germline_end",
                            "d_sequence_start", "d_sequence_end", "d_germline_start", "d_germline_end",
                            "j_sequence_start", "j_sequence_end", "j_germline_start", "j_germline_end",
                            "junction_length", "np1_length", "np2_length", "duplicate_count", "consensus_count"]

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


class IgorTask:
    """
    This class should encapsulate all
    the input parameters and output files when IGoR run.
    """

    def __init__(self, igor_exec_path=None, igor_datadir=None,
                 igor_models_root_path=None, igor_species=None, igor_chain=None,
                 igor_model_dir_path=None,
                 igor_path_ref_genome=None, fln_genomicVs=None, fln_genomicDs=None, fln_genomicJs=None,
                 fln_V_gene_CDR3_anchors=None, fln_J_gene_CDR3_anchors=None,
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
                 igor_fln_db=None,
                 mdl:Union[None,IgorModel] = None,
                 genomes:Union[None,IgorRefGenome] = None
                 ):
        # To execute IGoR externally
        self.igor_exec_path = igor_exec_path
        self.igor_datadir = igor_datadir

        # To load default models and genomic templates
        self.igor_models_root_path = igor_models_root_path  # igor models paths where all species and chains are stored.
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
        self.igor_fln_indexed_functionality = None # FIXME: ADD to batch update etc, etc, etc
        # inference
        self.igor_fln_infer_final_marginals = igor_fln_infer_final_marginals
        self.igor_fln_infer_final_parms = igor_fln_infer_final_parms
        self.igor_fln_infer_likelihoods = None
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

        # Temporary files used to get models in a local directory to generate sequences, for instance
        self.igor_mdldata_dir = None
        self.igor_fln_mdldata_parms = None
        self.igor_fln_mdldata_marginals = None
        self.igor_fln_mdldata_genomicVs = None
        self.igor_fln_mdldata_genomicDs = None
        self.igor_fln_mdldata_genomicJs = None
        self.igor_fln_mdldata_V_gene_CDR3_anchors = None
        self.igor_fln_mdldata_J_gene_CDR3_anchors = None


        self.igor_fln_db = igor_fln_db

        # TODO: experimental dictionary to check status of igor batch associated files
        # almost each of these files correspond to a sql table
        self.batch_data = igor_batch_dict

        self.igor_db = IgorSqliteDB(igor_fln_db=self.igor_fln_db)
        self.igor_db_bs = None

        self.b_read_seqs = False
        self.b_align = False
        self.b_infer = False
        self.b_evaluate = False
        self.b_generate = False

        self.mdl = mdl  # IgorModel()
        self.genomes = genomes  # IgorRefGenome() #{ 'V' : IgorRefGenome(), 'D' : IgorRefGenome(), 'J' : IgorRefGenome() }

        self.df_infer_likelihoods = None

        self.igor_align_dict_options = copy.deepcopy(igor_align_dict_options)

        self.igor_infer_dict_options = copy.deepcopy(igor_infer_dict_options)

        self.igor_evaluate_dict_options = copy.deepcopy(igor_evaluate_dict_options)

        self.igor_output_dict_options = copy.deepcopy(igor_output_dict_options)

        self.igor_generate_dict_options = copy.deepcopy(igor_generate_dict_options)

        try:
            if self.igor_batchname is None:
                self.gen_random_batchname()
        except Exception as e:
            e_message = ""
            import sys
            raise type(e)(str(e) + e_message).with_traceback(sys.exc_info()[2])

        try:
            if self.igor_wd is None:
                self.gen_igor_wd()
            else:
                # if not None and path doesnt exist create it.
                import pathlib
                pathlib.Path(self.igor_wd).mkdir(parents=True, exist_ok=True)
        except Exception as e:
            print(e)
            raise e

        try:
            if self.igor_exec_path is None:
                self.igor_exec_path = run_get_igor_exec_path()
        except Exception as e:
            print(e)
            raise e

        try:
            if self.igor_datadir is None:
                self.run_datadir()
        except Exception as e:
            print(e)
            raise e

        try:
            Q_species_chain = (self.igor_chain is not None) and (self.igor_species is not None)
            Q_model_parms = (self.igor_model_parms_file is not None)
            Q_fln_db = (self.igor_fln_db is not None)
            if True in [Q_species_chain, Q_model_parms, Q_fln_db]:
                self.update_model_filenames()
        except Exception as e:
            pass

        if self.mdl is None:
            try:
                self.load_mdl_from_db()
            except Exception as e:
                try:
                    self.load_IgorModel()
                except Exception as e:
                    # print("Model was not loaded!")
                    pass

        if self.genomes is None:
            try:
                # FIXME ADD DATABASE TRY
                self.load_IgorRefGenome()
            except Exception as e:
                pass



        # try:
        #     if self.mdl is None:
        #         self.load_mdl_from_db()
        #         self.load_IgorModel()
        # except Exception as e:
        #     print(e)
        #     raise e

    def __repr__(self):
        str_repr = ""
        for key, value in self.to_dict().items():
            str_repr = str_repr + str(key) + " = " + str(value) + "\n"
        return str_repr

    def to_dict(self):
        dicto = {
            "igor_exec_path": self.igor_exec_path,
            "igor_datadir": self.igor_datadir,
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
            "b_align": self.b_align,
            "b_infer": self.b_infer,
            "b_evaluate": self.b_evaluate,
            "b_generate": self.b_generate
        }
        return dicto

    def load_IgorRefGenome(self, igor_path_ref_genome=None):
        try:
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
        except Exception as e:
            raise e

    def make_model_default_VJ_from_genomes_dir(self, igor_path_ref_genome=None):
        try:
            self.load_IgorRefGenome(igor_path_ref_genome=igor_path_ref_genome)
            mdl_parms = IgorModel_Parms.make_default_VJ(self.genomes.df_genomicVs, self.genomes.df_genomicJs)
            mdl_marginals = IgorModel_Marginals.make_uniform_from_parms(mdl_parms)
            self.mdl = IgorModel.load_from_parms_marginals_object(mdl_parms, mdl_marginals)
        except Exception as e:
            print("ERROR: ", e)

    def make_model_default_VDJ_from_genomes_dir(self, igor_path_ref_genome=None):
        try:
            self.load_IgorRefGenome(igor_path_ref_genome=igor_path_ref_genome)
            mdl_parms = IgorModel_Parms.make_default_VDJ(self.genomes.df_genomicVs, self.genomes.df_genomicDs,
                                                         self.genomes.df_genomicJs)
            mdl_marginals = IgorModel_Marginals.make_uniform_from_parms(mdl_parms)
            self.mdl = IgorModel.load_from_parms_marginals_object(mdl_parms, mdl_marginals)
        except Exception as e:
            print("ERROR: ", e)

    def make_model_default_VDJ_from_fasta_files(self, fln_genomicVs: Union[None, str, Path] = None,
                                                fln_genomicJs: Union[None, str, Path] = None,
                                                fln_genomicDs: Union[None, str, Path] = None):
        """
        Make a default VDJ model from files
        """
        try:
            return 0
        except Exception as e:
            raise e

    def load_IgorModel(self, igor_model_parms_file: Union[None, str] = None,
                       igor_model_marginals_file: Union[None, str] = None,
                       fln_V_gene_CDR3_anchors: Union[None, str] = None,
                       fln_J_gene_CDR3_anchors: Union[None, str] = None):
        try:
            if igor_model_parms_file is not None:
                self.igor_model_parms_file = igor_model_parms_file
            if igor_model_marginals_file is not None:
                self.igor_model_marginals_file = igor_model_marginals_file

            if fln_V_gene_CDR3_anchors is not None:
                self.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors

            if fln_J_gene_CDR3_anchors is not None:
                self.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

            if ((self.igor_species is None) or (self.igor_chain is None)):
                # self.mdl = IgorModel.load_from_txt(self.igor_model_parms_file, self.igor_model_marginals_file)
                self.mdl = IgorModel(model_parms_file = self.igor_model_parms_file,
                                     model_marginals_file=self.igor_model_marginals_file,
                                     fln_V_gene_CDR3_anchors= self.fln_V_gene_CDR3_anchors,
                                     fln_J_gene_CDR3_anchors= self.fln_J_gene_CDR3_anchors)
            else:
                try:
                    # self.mdl = IgorModel.load_default(self.igor_species, igor_option_path_dict[self.igor_chain])
                    self.mdl = IgorModel.load_default(self.igor_species, igor_option_path_dict[self.igor_chain])
                except KeyError as ke:
                    self.mdl = IgorModel.load_default(self.igor_species, self.igor_chain)
                except Exception as e:
                    raise e

            # print("Model loaded")
        except Exception as e:
            e_message = "WARNING: IgorTask.load_IgorModel" + str(self)
            import sys
            raise type(e)(str(e) + e_message).with_traceback(sys.exc_info()[2])

    def load_IgorModel_from_infer_files(self, igor_fln_infer_final_parms: Union[None, str] = None,
                                        igor_fln_infer_final_marginals: Union[None, str] = None,
                                        fln_V_gene_CDR3_anchors: Union[None, str] = None,
                                        fln_J_gene_CDR3_anchors: Union[None, str] = None):
        """
        Load IgorModel from inferred model files.
        :param igor_fln_infer_final_parms: Path of inferred model parms file.
        :param igor_fln_infer_final_marginals: Path of inferred model marginals file
        """
        try:
            if igor_fln_infer_final_parms is not None:
                self.igor_fln_infer_final_parms = igor_fln_infer_final_parms

            if igor_fln_infer_final_marginals is not None:
                self.igor_fln_infer_final_marginals = igor_fln_infer_final_marginals

            if fln_V_gene_CDR3_anchors is not None:
                self.fln_V_gene_CDR3_anchors = fln_V_gene_CDR3_anchors

            if fln_J_gene_CDR3_anchors is not None:
                self.fln_J_gene_CDR3_anchors = fln_J_gene_CDR3_anchors

            import pathlib
            # import os
            # if os.path.isfile(self.fln_V_gene_CDR3_anchors):
            if self.fln_V_gene_CDR3_anchors is not None:
                if not pathlib.Path(self.fln_V_gene_CDR3_anchors).is_file():
                    self.fln_V_gene_CDR3_anchors = None

            # if os.path.isfile(self.fln_J_gene_CDR3_anchors):
            if self.fln_J_gene_CDR3_anchors is not None:
                if not pathlib.Path(self.fln_J_gene_CDR3_anchors).is_file():
                    self.fln_J_gene_CDR3_anchors = None


            self.mdl = IgorModel(model_parms_file=self.igor_fln_infer_final_parms,
                                 model_marginals_file=self.igor_fln_infer_final_marginals,
                                 fln_V_gene_CDR3_anchors=self.fln_V_gene_CDR3_anchors,
                                 fln_J_gene_CDR3_anchors=self.fln_J_gene_CDR3_anchors)

            try:
                self.df_infer_likelihoods = pd.read_csv(self.igor_fln_infer_likelihoods, sep=';')
            except Exception as e:
                print("Likelihoods files not found: ", self.igor_fln_infer_likelihoods)
                raise e

        except Exception as e:
            e_message = "ERROR: IgorTask.load_IgorModel_inferred"
            # print(e)
            raise e
        else:
            return self.mdl

    @classmethod
    def default_model(cls, specie, chain, igor_wd=None, model_parms_file=None, model_marginals_file=None, **kwargs):
        """Return an IgorTask object"""
        try:
            cls = IgorTask()
            cls.igor_species = specie
            cls.igor_chain = igor_option_path_dict[chain]
            if igor_wd is not None:
                cls.igor_wd = igor_wd
            # cls.igor_modeldirpath =  model_parms_file
            cls.run_datadir()
            cls.igor_model_dir_path = cls.igor_models_root_path + "/" + cls.igor_species + "/" + cls.igor_chain
            cls.update_model_filenames(igor_model_dir_path=cls.igor_model_dir_path)
            cls.igor_path_ref_genome = cls.igor_model_dir_path + "/" + "ref_genome"
            cls.update_ref_genome()

            if model_parms_file is None:
                cls.igor_model_parms_file = cls.igor_model_dir_path + "/models/model_parms.txt"
                cls.igor_model_marginals_file = cls.igor_model_dir_path + "/models/model_marginals.txt"
                cls.mdl = IgorModel(model_parms_file=cls.igor_model_parms_file,
                                    model_marginals_file=cls.igor_model_marginals_file,
                                    fln_V_gene_CDR3_anchors=cls.fln_V_gene_CDR3_anchors,
                                    fln_J_gene_CDR3_anchors=cls.fln_J_gene_CDR3_anchors)
                cls.load_IgorRefGenome()
            return cls
        except Exception as e:
            raise e

    def gen_igor_wd(self):
        # p = subprocess.run("pwd", shell=True, capture_output=True, text=True)
        # line = p.stdout.readline()
        # self.igor_wd = line.decode("utf-8").replace('\n', '')
        self.igor_wd = run_get_igor_wd()

    def gen_random_batchname(self):
        try:
            # p = subprocess.Popen("head /dev/urandom | tr -dc A-Za-z0-9 | head -c10", shell=True, stdout=subprocess.PIPE)
            # line = p.stdout.readline()
            # self.igor_batchname = "dataIGoR" + line.decode("utf-8").replace('\n', '')
            str_random = run_get_random_string()
            self.igor_batchname = "dataIGoR" + str_random
        except Exception as e:
            raise e

    def update_model_filenames(self, igor_model_dir_path: Union[None, str] = None,
                               olga_model_dir_path: Union[None, str] = None,
                               igor_models_root_path: Union[None, str] = None):
        """Update model filenames
            :param igor_model_dir_path: Directory path for genome templates.
            If None default is igor_model_dir_path = igor_models_root_path + "/models"
            :param igor_models_root_path: Directory path where different species and chain models.
            If None don't change default value is get it from run_datadir()
            "$(igor -getdatadir)/models/"

        """

        try:
            if igor_models_root_path is not None:
                self.igor_models_root_path = igor_models_root_path

            # if model_path is None use the self.igor_model_dir_path
            if igor_model_dir_path is None:
                # use previously defined igor_model_dir_path
                if self.igor_model_dir_path is None:
                    # if wasn't defined use the current directory
                    igor_model_dir_path = "."
                    if (not (self.igor_species is None)) and (not (self.igor_chain is None)):
                        self.run_datadir()
                        self.igor_model_dir_path = self.igor_models_root_path + "/" + self.igor_species + "/" + \
                                                   igor_option_path_dict[self.igor_chain]
                    else:
                        self.igor_model_dir_path = igor_model_dir_path
            else:
                # if a model_path is provided then override it
                self.igor_model_dir_path = igor_model_dir_path

            if olga_model_dir_path is not None:
                self.igor_model_parms_file = olga_model_dir_path + "/model_params.txt"
                self.igor_model_marginals_file = olga_model_dir_path + "/model_marginals.txt"
                self.igor_path_ref_genome = olga_model_dir_path
            else:

                self.igor_model_parms_file = self.igor_model_dir_path + "/models/model_parms.txt"
                self.igor_model_marginals_file = self.igor_model_dir_path + "/models/model_marginals.txt"
                self.igor_path_ref_genome = self.igor_model_dir_path + "/ref_genome/"
        except Exception as e:
            e_message = "WARNING: IgorTask.update_model_filenames: " + str(self.igor_model_dir_path)
            import sys
            raise type(e)(str(e) + e_message).with_traceback(sys.exc_info()[2])

    def update_ref_genome(self, igor_path_ref_genome: Union[None, str] = None,
                          igor_model_dir_path: Union[None, str] = None,
                          genomes: Union[None, IgorRefGenome] = None,
                          fln_genomicVs: Union[None, str] = None,
                          fln_genomicDs: Union[None, str] = None,
                          fln_genomicJs: Union[None, str] = None,
                          fln_V_gene_CDR3_anchors: Union[None, str] = None,
                          fln_J_gene_CDR3_anchors: Union[None, str] = None):
        """Assign names to ref_genome files gene templates and CDR3 anchors
        :param igor_path_ref_genome: Directory path for genome templates.
        If None igor_path_ref_genome = igor_model_dir_path + "/ref_genome"
        :param igor_model_dir_path: Character to delimitate csv file.
        :param genomes:Union[None,IgorRefGenome]
        """
        try:
            if igor_model_dir_path is not None:
                self.igor_model_dir_path = igor_model_dir_path
            # else:
                # self.igor_model_dir_path = self.igor_wd + "/" + self.igor_batchname + "_mdldata"

            if igor_path_ref_genome is not None:
                self.igor_path_ref_genome = igor_path_ref_genome
            else:
                self.igor_path_ref_genome = self.igor_model_dir_path + "/ref_genome"  # default path

            if genomes is not None:
                self.genomes = genomes

            if self.genomes is None:
                self.genomes = IgorRefGenome()

            self.genomes.update_fln_names(path_ref_genome=self.igor_path_ref_genome,
                                          fln_genomicVs=fln_genomicVs, fln_genomicDs=fln_genomicDs,
                                          fln_genomicJs=fln_genomicJs,
                                          fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                                          fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors)

            self.fln_genomicVs = self.genomes.fln_genomicVs
            self.fln_genomicJs = self.genomes.fln_genomicJs
            self.fln_genomicDs = self.genomes.fln_genomicDs
            self.fln_V_gene_CDR3_anchors = self.genomes.fln_V_gene_CDR3_anchors
            self.fln_J_gene_CDR3_anchors = self.genomes.fln_J_gene_CDR3_anchors
        except Exception as e:
            e_message = "ERROR: IgorTask.update_ref_genome \n" + str(self.to_dict())
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def update_batch_filenames(self, igor_batchname=None, igor_wd=None):
        # reads
        try:
            if igor_batchname is not None:
                self.igor_batchname = igor_batchname

            if igor_wd is not None:
                self.igor_wd = igor_wd

            if self.igor_wd is None:
                self.gen_igor_wd()

            if self.igor_batchname is None:
                self.gen_random_batchname()

            self._update_align_batch_filenames(igor_batchname=igor_batchname, igor_wd=igor_wd)

            self._update_infer_batch_filenames(igor_batchname=igor_batchname, igor_wd=igor_wd)

            self._update_evaluate_batch_filenames(igor_batchname=igor_batchname, igor_wd=igor_wd)

            self._update_output_batch_filenames(igor_batchname=igor_batchname, igor_wd=igor_wd)

            # Set all files as not existing by default
            import os.path
            for file_id in igor_file_id_list:
                self.batch_data[file_id]['status'] = os.path.isfile(self.batch_data[file_id]['filename'])
            # database
            self.igor_fln_db = self.igor_wd + "/" + self.igor_batchname + ".db"

            tmp_prefix_aligns = self.igor_wd + "/aligns/" + self.igor_batchname
            self.batch_data['indexed_sequences']['filename'] = tmp_prefix_aligns + "_indexed_sequences.csv"
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
        except Exception as e:
            e_message = "ERROR: IgorTask.update_batch_filenames "
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

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

    def _update_align_batch_filenames(self, igor_batchname=None, igor_wd=None):
        """Update align filenames using batchname and igor_wd"""
        self.igor_fln_indexed_sequences = self.igor_wd + "/aligns/" + self.igor_batchname + "_indexed_sequences.csv"
        # aligns
        self.igor_fln_indexed_CDR3 = self.igor_wd + "/aligns/" + self.igor_batchname + "_indexed_CDR3s.csv"

        self.igor_fln_align_V_alignments = self.igor_wd + "/aligns/" + self.igor_batchname + "_V_alignments.csv"
        self.igor_fln_align_J_alignments = self.igor_wd + "/aligns/" + self.igor_batchname + "_J_alignments.csv"
        self.igor_fln_align_D_alignments = self.igor_wd + "/aligns/" + self.igor_batchname + "_D_alignments.csv"

        self.igor_fln_indexed_functionality = self.igor_wd + "/aligns/" + self.igor_batchname + "_indexed_functionality.csv"


    def _update_infer_batch_filenames(self, igor_batchname=None, igor_wd=None):
        """Update inference filenames using batchname and igor_wd"""
        # inference
        tmpstr = self.igor_wd + "/" + self.igor_batchname + "_inference/"
        self.igor_fln_infer_final_parms = tmpstr + "final_parms.txt"
        self.igor_fln_infer_final_marginals = tmpstr + "final_marginals.txt"
        self.igor_fln_infer_likelihoods = tmpstr + "likelihoods.out"

    def _update_evaluate_batch_filenames(self, igor_batchname=None, igor_wd=None):
        """Update evaluate filenames using batchname and igor_wd"""
        # evaluate
        tmpstr = self.igor_wd + "/" + self.igor_batchname + "_evaluate/"
        self.igor_fln_evaluate_final_parms = tmpstr + "final_parms.txt"
        self.igor_fln_evaluate_final_marginals = tmpstr + "final_marginals.txt"

    def _update_output_batch_filenames(self, igor_batchname=None, igor_wd=None):
        """Update output filenames using batchname and igor_wd"""
        # output
        tmpstr = self.igor_wd + "/" + self.igor_batchname + "_output/"
        self.igor_fln_output_pgen = tmpstr + "Pgen_counts.csv"
        self.igor_fln_output_scenarios = tmpstr + "best_scenarios_counts.csv"
        self.igor_fln_output_coverage = tmpstr + "coverage.csv"

    def _update_generate_batch_filenames(self, igor_batchname=None, igor_wd=None):
        """Update generate filenames using batchname and igor_wd"""
        # generate
        tmpstr = self.igor_wd + "/" + self.igor_batchname + "_generated/"
        self.igor_fln_generated_realizations_werr = tmpstr + "generated_realizations_werr.csv"
        self.igor_fln_generated_seqs_werr = tmpstr + "generated_seqs_werr.csv"
        self.igor_fln_generation_info = tmpstr + "generated_seqs_werr.out"

    def _update_mdldata_batch_filenames(self):
        self.igor_mdldata_dir = self.igor_wd + "/" + self.igor_batchname + "_mdldata/"
        fln_dict = get_default_fln_names_for_model_dir(self.igor_mdldata_dir)

        self.igor_fln_mdldata_parms = fln_dict['fln_model_parms']
        self.igor_fln_mdldata_marginals = fln_dict['fln_model_marginals']
        self.igor_fln_mdldata_genomicVs = fln_dict['fln_genomicVs']
        self.igor_fln_mdldata_genomicDs = fln_dict['fln_genomicDs']
        self.igor_fln_mdldata_genomicJs = fln_dict['fln_genomicJs']

        self.igor_fln_mdldata_V_gene_CDR3_anchors = fln_dict['fln_V_gene_CDR3_anchors']
        self.igor_fln_mdldata_J_gene_CDR3_anchors = fln_dict['fln_J_gene_CDR3_anchors']


    @classmethod
    def load_from_batchname(cls, batchname, wd=None, ):
        cls = IgorTask()
        cls.igor_path_ref_genome
        cls.igor_model_dir_path
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
        cmd = self.igor_exec_path + " -run_demo"
        return run_command_print(cmd)
        # return run_command(cmd)

    def run_datadir(self):
        # cmd = self.igor_exec_path + " -getdatadir"
        # self.igor_datadir = run_command(cmd).replace('\n', '')
        # self.igor_models_root_path = self.igor_datadir + "/models/"
        self.igor_datadir = run_get_igor_datadir()
        self.igor_models_root_path = self.igor_datadir + "/models/"

    def run_read_seqs(self, igor_read_seqs=None):
        """
        Run IGoR's -read_seqs options
        """
        try:
            #TODO: FIXME igor_read_seqs is different that the input sequences
            if igor_read_seqs is not None:
                self.igor_read_seqs = igor_read_seqs

            import pathlib
            pathlib.Path(self.igor_wd).mkdir(parents=True, exist_ok=True)

            "igor -set_wd $WDPATH -batch foo -read_seqs ../demo/murugan_naive1_noncoding_demo_seqs.txt"
            cmd = self.igor_exec_path
            cmd = cmd + " -set_wd " + self.igor_wd
            cmd = cmd + " -batch " + self.igor_batchname
            cmd = cmd + " -read_seqs " + self.igor_read_seqs
            # TODO: if self.igor_read_seqs extension fastq then convert to csv and copy and create the file in aligns. Overwrite if necesserasy
            print(cmd)
            cmd_stdout = run_command(cmd)
            # subprocess.run(cmd, shell=True, capture_output=True, text=True)
            # cmd_stdout = run_command_print(cmd)
            self.igor_fln_indexed_sequences = self.igor_wd + "/aligns/" + self.igor_batchname + "_indexed_sequences.csv"
            self.b_read_seqs = True  # FIXME: If run_command success then True
            # return cmd_stdout
        except Exception as e:
            raise e

    def run_align(self, igor_read_seqs=None, igor_threads=None, igor_align_dict_opts:Union[dict, None]=None):
        # "igor -set_wd ${tmp_dir} -batch ${randomBatch} -species
        # ${species} -chain ${chain} -align --all"
        try:
            import os.path

            if igor_align_dict_opts is not None:
                self.igor_align_dict_options = igor_align_dict_opts

            if igor_read_seqs is not None:
                self.igor_read_seqs = igor_read_seqs

            if self.b_read_seqs is False:
                self.run_read_seqs(igor_read_seqs=igor_read_seqs)

            import pathlib
            pathlib.Path(self.igor_wd).mkdir(parents=True, exist_ok=True)

            if self.igor_mdldata_dir is not None:

                cmd = self.igor_exec_path
                if self.igor_threads is not None:
                    cmd = cmd + " -threads " + str(int(self.igor_threads))
                cmd = cmd + " -set_wd " + self.igor_wd
                cmd = cmd + " -batch " + self.igor_batchname
                # I think that the safests is to use the
                cmd = cmd + " -set_genomic "

                if os.path.isfile(self.igor_fln_mdldata_genomicVs):
                    cmd = cmd + " --V " + self.igor_fln_mdldata_genomicVs

                if self.igor_fln_mdldata_genomicDs is not None:
                    if os.path.isfile(self.igor_fln_mdldata_genomicDs):
                        cmd = cmd + " --D " + self.igor_fln_mdldata_genomicDs
                if os.path.isfile(self.igor_fln_mdldata_genomicJs):
                    cmd = cmd + " --J " + self.igor_fln_mdldata_genomicJs

                if os.path.isfile(self.igor_fln_mdldata_V_gene_CDR3_anchors) or \
                        os.path.isfile(self.igor_fln_mdldata_V_gene_CDR3_anchors):
                    cmd = cmd + " -set_CDR3_anchors "

                    if os.path.isfile(self.igor_fln_mdldata_V_gene_CDR3_anchors):
                        cmd = cmd + " --V " + self.igor_fln_mdldata_V_gene_CDR3_anchors

                    if os.path.isfile(self.igor_fln_mdldata_J_gene_CDR3_anchors):
                        cmd = cmd + " --J " + self.igor_fln_mdldata_J_gene_CDR3_anchors

                cmd = cmd + " -align " + command_from_dict_options(self.igor_align_dict_options)

            else:
                cmd = self.igor_exec_path
                if self.igor_threads is not None:
                    cmd = cmd + " -threads " + str(int(self.igor_threads))
                cmd = cmd + " -set_wd " + self.igor_wd
                cmd = cmd + " -batch " + self.igor_batchname
                # I think that the safests is to use the
                cmd = cmd + " -set_genomic "

                if os.path.isfile(self.genomes.fln_genomicVs):
                    cmd = cmd + " --V " + self.genomes.fln_genomicVs
                if self.genomes.fln_genomicDs is not None:
                    if os.path.isfile(self.genomes.fln_genomicDs):
                        cmd = cmd + " --D " + self.genomes.fln_genomicDs
                if os.path.isfile(self.genomes.fln_genomicJs):
                    cmd = cmd + " --J " + self.genomes.fln_genomicJs

                if os.path.isfile(self.genomes.fln_V_gene_CDR3_anchors) or \
                        os.path.isfile(self.genomes.fln_J_gene_CDR3_anchors):
                    cmd = cmd + " -set_CDR3_anchors "

                    if os.path.isfile(self.genomes.fln_V_gene_CDR3_anchors):
                        cmd = cmd + " --V " + self.genomes.fln_V_gene_CDR3_anchors

                    if os.path.isfile(self.genomes.fln_J_gene_CDR3_anchors):
                        cmd = cmd + " --J " + self.genomes.fln_J_gene_CDR3_anchors

                cmd = cmd + " -align " + command_from_dict_options(self.igor_align_dict_options)


            # cmd = self.igor_exec_path
            # cmd = cmd + " -set_wd " + self.igor_wd
            # cmd = cmd + " -batch " + self.igor_batchname
            # # TODO: USE COSTUM MODEL OR USE SPECIFIED SPECIES?
            # # I think that the safests is to use the
            # # FIXME: CHANGE TO CUSTOM GENOMICS
            # cmd = cmd + " -set_genomic "
            #
            # if os.path.isfile(self.genomes.fln_genomicVs):
            #     cmd = cmd + " --V " + self.genomes.fln_genomicVs
            # if os.path.isfile(self.genomes.fln_genomicDs):
            #     cmd = cmd + " --D " + self.genomes.fln_genomicDs
            # if os.path.isfile(self.genomes.fln_genomicJs):
            #     cmd = cmd + " --J " + self.genomes.fln_genomicJs
            #
            # cmd = cmd + " -set_CDR3_anchors "
            #
            # if os.path.isfile(self.genomes.fln_V_gene_CDR3_anchors):
            #     cmd = cmd + " --V " + self.genomes.fln_V_gene_CDR3_anchors
            # if os.path.isfile(self.genomes.fln_J_gene_CDR3_anchors):
            #     cmd = cmd + " --J " + self.genomes.fln_J_gene_CDR3_anchors
            #
            # cmd = cmd + " -align " + command_from_dict_options(self.igor_align_dict_options)
            # # return cmd

            print(cmd)
            cmd_stdout = run_command_print(cmd)
            # run_command_no_output(cmd)
            self._update_align_batch_filenames()
            self.b_align = True  # FIXME: If run_command success then True
            return cmd_stdout
        except Exception as e:
            raise e

    # FIXME: FINISH METHOD
    def align(self, input_sequences: Union[None, str, Path, pd.DataFrame, np.array, list] = None,
              mdl: Union[None, str, Path, IgorModel, IgorModel_Parms] = None,
              igor_wd=None, batch_clean=True, b_best_align_only=False,
              igor_align_dict_opts:Union[None, dict]=None):
        df_alignments_dict = {'V': None, 'D': None, 'J': None, 'CDR3': None}
        try:
            # by default the igor_wd is the current directory unless something else is
            tmp_evaluate_dir = tempfile.TemporaryDirectory(prefix='igor_aligning_', dir='.')
            if igor_wd is None:
                igor_wd = tmp_evaluate_dir.name
            else:
                igor_wd = self.igor_wd

            # 3. Write Sequences in file if file not exist
            # FIXME: -IN DEV- CHANGE THIS IN A WAY TO NOT DELETE igor_read_seqs.

            # if a file is created during the process delete it

            # 4. Export model and ref_genome to model_dir
            self.write_mdldata_dir()
            # path_mdl_data = self.igor_wd + "/" + self.igor_batchname + "_mdldata"
            self.update_model_filenames(igor_model_dir_path=self.igor_mdldata_dir)
            self.update_ref_genome(igor_model_dir_path=self.igor_mdldata_dir)
            self.update_batch_filenames()
            # self.mdl.write_mdldata_dir(path_mdl_data)

            # 5. Run infer model

            import copy
            if mdl is not None:
                self.mdl = copy.deepcopy(mdl)

            # ---best_align_only

            if input_sequences is None:
                # use the self.igor_read_seqs
                # Do not delete self.igor_read_seqs
                self.run_align(igor_read_seqs=self.igor_read_seqs, igor_align_dict_opts=igor_align_dict_opts)
            else:
                # write a temporary file
                tmp_igor_read_seqs = self.igor_read_seqs
                fln_input_sequences = self.igor_wd + "/" + self.igor_batchname + "_input_sequences.csv"
                write_sequences_to_file(input_sequences, fln_input_sequences)

                # 4. Export model and ref_genome to model_dir
                path_mdl_data = igor_wd + "/" + self.igor_batchname + "_mdldata"
                self.update_model_filenames(igor_model_dir_path=path_mdl_data)
                self.update_ref_genome(igor_model_dir_path=path_mdl_data)
                self.update_batch_filenames()
                self._update_mdldata_batch_filenames()
                if self.mdl.parms.event_GeneChoice_D is None:
                    self.igor_fln_mdldata_genomicDs = None
                self.mdl.write_mdldata_dir(path_mdl_data)

                self.run_align(igor_read_seqs=fln_input_sequences, igor_align_dict_opts=igor_align_dict_opts)
                import os
                os.unlink(fln_input_sequences)
                self.igor_read_seqs = tmp_igor_read_seqs

                # self.igor_fln_indexed_sequences
                # self.igor_fln_align_V_alignments
                # self.igor_fln_align_D_alignments
                # self.igor_fln_align_J_alignments
                # self.igor_fln_indexed_CDR3
                if os.path.isfile(self.igor_fln_align_V_alignments):
                    df_V_aligments = pd.read_csv(self.igor_fln_align_V_alignments, sep=';')
                    df_alignments_dict['V'] = df_V_aligments
                if os.path.isfile(self.igor_fln_align_D_alignments):
                    df_D_aligments = pd.read_csv(self.igor_fln_align_D_alignments, sep=';')
                    df_alignments_dict['D'] = df_D_aligments
                if os.path.isfile(self.igor_fln_align_J_alignments):
                    df_J_aligments = pd.read_csv(self.igor_fln_align_J_alignments, sep=';')
                    df_alignments_dict['J'] = df_J_aligments

        except Exception as e:
            raise e
        else:
            # TODO: RETURN ALIGNMENT? OR CDR3 ALIGN ONLY?
            return df_alignments_dict
        finally:
            self._run_clean_batch_aligns()

    def _run_evaluate(self, igor_read_seqs=None, N_scenarios=None, Pgen=True,
                      igor_model_parms_file=None,
                      igor_model_marginals_file=None,
                      fln_V_gene_CDR3_anchors=None,
                      fln_J_gene_CDR3_anchors=None,
                      igor_fln_db=None,
                      igor_db=None,
                      mdl: Union[IgorModel, None] = None, igor_threads=None,
                      igor_evaluate_dict_opts:Union[dict, None]=None,
                      igor_align_dict_opts:Union[dict, None]=None):
        # "igor -set_wd $WDPATH -batch foo -species human -chain beta
        # -evaluate -output --scenarios 10"
        try:
            if igor_evaluate_dict_opts is not None:
                self.igor_evaluate_dict_options = igor_evaluate_dict_opts

            if igor_align_dict_opts is not None:
                self.igor_align_dict_options = igor_align_dict_opts

            # print(self.to_dict())
            import os.path
            if igor_read_seqs is not None:
                self.igor_read_seqs = igor_read_seqs

            if mdl is not None:
                self.mdl = mdl

            if self.mdl is None:
                try:
                    self.load_IgorModel(igor_model_parms_file=igor_model_parms_file,
                                        igor_model_marginals_file=igor_model_marginals_file,
                                        fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                                        fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors)
                except:
                    try:
                        self.load_mdl_from_db(igor_fln_db=igor_fln_db, igor_db=igor_db)
                    except:
                        pass
                    pass
            else:
                self._update_mdldata_batch_filenames()
                if self.mdl.parms.event_GeneChoice_D is None:
                    self.igor_fln_mdldata_genomicDs = None
                self.write_mdldata_dir(self.igor_mdldata_dir)

            import pathlib
            pathlib.Path(self.igor_wd).mkdir(parents=True, exist_ok=True)

            if self.b_align is False:
                try:
                    self.run_align(igor_read_seqs=self.igor_read_seqs, igor_align_dict_opts=igor_align_dict_opts)
                except Exception as e:
                    raise e


            if self.mdl is not None:
                self._update_mdldata_batch_filenames()
                self.write_mdldata_dir(self.igor_mdldata_dir)

                cmd = self.igor_exec_path
                if self.igor_threads is not None:
                    cmd = cmd + " -threads " + str(int(self.igor_threads))
                cmd = cmd + " -set_wd " + self.igor_wd
                cmd = cmd + " -batch " + self.igor_batchname
                cmd = cmd + " -set_custom_model " + self.igor_fln_mdldata_parms + " " + self.igor_fln_mdldata_marginals

                if os.path.isfile(self.igor_fln_mdldata_V_gene_CDR3_anchors) or \
                        os.path.isfile(self.igor_fln_mdldata_J_gene_CDR3_anchors):
                    cmd = cmd + " -set_CDR3_anchors "
                    if os.path.isfile(self.igor_fln_mdldata_V_gene_CDR3_anchors):
                        cmd = cmd + " --V " + self.igor_fln_mdldata_V_gene_CDR3_anchors
                    if os.path.isfile(self.igor_fln_mdldata_J_gene_CDR3_anchors):
                        cmd = cmd + " --J " + self.igor_fln_mdldata_J_gene_CDR3_anchors
                # here the evaluation
                self.igor_output_dict_options["--scenarios"]['active'] = True
                if N_scenarios is not None:
                    self.igor_output_dict_options["--scenarios"]['value'] = str(N_scenarios)
                self.igor_output_dict_options["--Pgen"]['active'] = Pgen
                cmd = cmd + " -evaluate " + command_from_dict_options(self.igor_evaluate_dict_options)
                cmd = cmd + " -output " + command_from_dict_options(self.igor_output_dict_options)

            else:
                cmd = self.igor_exec_path
                cmd = cmd + " -set_wd " + self.igor_wd
                cmd = cmd + " -batch " + self.igor_batchname
                cmd = cmd + " -set_custom_model " + self.igor_model_parms_file + " " + self.igor_model_marginals_file

                # here the evaluation
                self.igor_output_dict_options["--scenarios"]['active'] = True
                if N_scenarios is not None:
                    self.igor_output_dict_options["--scenarios"]['value'] = str(N_scenarios)
                self.igor_output_dict_options["--Pgen"]['active'] = Pgen
                cmd = cmd + " -evaluate " + command_from_dict_options(self.igor_evaluate_dict_options)
                cmd = cmd + " -output " + command_from_dict_options(self.igor_output_dict_options)


            print(cmd)
            # FIXME: REALLY BIG FLAW USE DICTIONARY FOR THE SPECIE AND CHAIN
            # self.mdl = IgorModel.load_default(self.igor_species, igor_option_path_dict[self.igor_chain], modelpath=self.igor_models_root_path)
            run_command(cmd)
            self._update_evaluate_batch_filenames()
            # run_command_no_output(cmd)
            # self.b_evaluate = True # FIXME: If run_command success then Truerun_infer
        except Exception as e:
            raise e

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
        self.igor_output_dict_options["--scenarios"]['active'] = True
        if N_scenarios is not None:
            self.igor_output_dict_options["--scenarios"]['value'] = str(N_scenarios)
        self.igor_output_dict_options["--Pgen"]['active'] = False
        cmd = cmd + " -evaluate " + command_from_dict_options(self.igor_evaluate_dict_options)
        cmd = cmd + " -output " + command_from_dict_options(self.igor_output_dict_options)
        # return cmd
        print(cmd)
        # FIXME: REALLY BIG FLAW USE DICTIONARY FOR THE SPECIE AND CHAIN
        # self.mdl = IgorModel.load_default(self.igor_species, igor_option_path_dict[self.igor_chain], modelpath=self.igor_models_root_path)
        # run_command(cmd)
        run_command_print(cmd)

    def generate(self, N_seqs=None, mdl=None,
                 igor_wd=None, igor_batchname=None,
                 igor_model_parms_file=None, igor_model_marginals_file=None,
                 igor_db=None, igor_fln_db=None,
                 igor_species=None, igor_chain=None, clean_batch=True, return_scenarios=False, return_df=True):
        """
        Generate Sequences using IgorTask
        """
        try:
            if N_seqs is None:
                N_seqs = 1

            if igor_wd is not None:
                self.igor_wd = igor_wd # tmp_generate_dir.name

            # with tempfile.TemporaryDirectory(prefix='igor_generating_', dir='.') as tmp_generate_dirname:
            #     self.igor_wd = tmp_generate_dirname
            if (igor_species is not None) and (igor_chain is not None):
                try:
                    self.mdl = IgorModel.load_default(igor_species, igor_chain)
                except Exception as e:
                    raise e

            if mdl is not None:
                self.mdl = mdl

            self.update_batch_filenames()

            self.igor_mdldata_dir = self.igor_wd + "/" + self.igor_batchname + "_mdldata"
            self.write_mdldata_dir(self.igor_mdldata_dir)
            # TODO: SHOULD I UPDATE HERE THE VARIABLES igor_fln_mdldata_genomicVs,
            #  igor_fln_mdldata_V_gene_CDR3_anchors, igor_fln_mdl_parms
            self.update_model_filenames(self.igor_mdldata_dir)
            self.update_ref_genome(self.igor_mdldata_dir)

            self._run_generate(N_seqs)
            # TODO: ADD COLUMNS OF EVENTS
            pd_sequences = get_dataframe_from_fln_generated_seqs_werr(self.igor_fln_generated_seqs_werr)

        except Exception as e:
            raise e
        else:
            if return_scenarios:
                df_scenarios = self.mdl.get_dataframe_from_fln_generated_realizations_werr(
                    self.igor_fln_generated_realizations_werr)
                return pd_sequences.join(df_scenarios)
            else:
                return pd_sequences  # mdl_inferrred

        finally:
            if clean_batch:
                self._run_clean_batch_generate()
                self._run_clean_batch_mdldata()
                # self.run_clean_batch()

    def get_dataframe_from_fln_generated_seqs_werr(self, igor_fln_generated_seqs_werr=None):
        if igor_fln_generated_seqs_werr is not None:
            self.igor_fln_generated_seqs_werr = igor_fln_generated_seqs_werr

        return get_dataframe_from_fln_generated_seqs_werr(self.igor_fln_generated_seqs_werr)

    def get_dataframe_from_fln_generated_realizations_werr(self, igor_fln_generated_realizations_werr=None,
                                                           mdl:Union[None, IgorModel]=None):
        if igor_fln_generated_realizations_werr is not None:
            self.igor_fln_generated_realizations_werr = igor_fln_generated_realizations_werr

        if mdl is not None:
            self.mdl = mdl

        return self.mdl.get_dataframe_from_fln_generated_realizations_werr(self.igor_fln_generated_realizations_werr)

    def evaluate(self, input_sequences: Union[str, pd.DataFrame, np.ndarray, Path],
                 N_scenarios = None, mdl:IgorModel = None, igor_wd=None, igor_batchname=None, clean_batch=True, airr_format=True, db=True,
                 igor_evaluate_dict_options:Union[dict,None]=None):
        """
        Return evaluation of sequences
        """
        try:
            # by default the igor_wd is the current directory unless something else is
            tmp_evaluate_dir = tempfile.TemporaryDirectory(prefix='igor_evaluating_', dir='.')
            if igor_wd is None:
                igor_wd = tmp_evaluate_dir.name
            else:
                igor_wd = self.igor_wd

            if igor_batchname is not None:
                self.igor_batchname = igor_batchname

            # FIXME: WHAT HAPPEN IF I WANT TO PRESERVE THE WORKING DIRECTORY?????

            # self.igor_wd = igor_wd

            # 2. Copy model to IgorTask
            import copy
            if mdl is not None:
                self.mdl = copy.deepcopy(mdl)

            # 3. Write Sequences in file if file not exist
            fln_input_sequences = igor_wd + "/" + self.igor_batchname + "input_sequences.csv"
            write_sequences_to_file(input_sequences, fln_input_sequences)

            # 4. Export model and ref_genome to model_dir
            path_mdl_data = igor_wd + "/" + self.igor_batchname + "_mdldata"
            self.update_model_filenames(igor_model_dir_path=path_mdl_data)
            self.update_ref_genome(igor_model_dir_path=path_mdl_data)
            self.update_batch_filenames()
            self.mdl.write_mdldata_dir(path_mdl_data)

            # 5. Run evaluate model
            self._run_evaluate(igor_read_seqs=fln_input_sequences, N_scenarios=N_scenarios)

            # Save evaluations in database
            try:
                if db:
                    self.create_db()
                    self.load_db()

                    base_fln_output = self.igor_fln_db.split(".db")[0]
                    output_fln_prefix = base_fln_output + "_rearrangement"
                    output_fln_airr = output_fln_prefix + ".tsv"
                    if airr_format:
                        self.igor_db.export_IgorBestScenarios_to_AIRR(output_fln_airr)
                        pd_airr_rearrangement = pd.read_csv(output_fln_airr, sep='\t')
                else:
                    pd_airr_rearrangement = self.mdl.get_dataframe_from_fln_generated_realizations_werr(
                        self.igor_fln_output_scenarios)
            except Exception as e:
                raise e

        except Exception as e:
            raise e
        else:
            return pd_airr_rearrangement
        finally:
            tmp_evaluate_dir.cleanup()
            if clean_batch:
                self.run_clean_batch()


    def infer(self, input_sequences: Union[None, str, Path, pd.DataFrame, np.array, list] = None,
              model: Union[None, str, Path, IgorModel, IgorModel_Parms] = None,
              igor_wd=None, batch_clean=True):
        """Run igor infer with new data and model
        :param input_sequences: Union[None, str, Path, pd.DataFrame, np.array, list] = None
        :param model: Union[None, str, Path, IgorModel, IgorModel_Parms] = None
        """
        try:
            if igor_wd is not None:
                self.igor_wd = igor_wd
            # 3. Write Sequences in file if file not exist
            # FIXME: -IN DEV- CHANGE THIS IN A WAY TO NOT DELETE igor_read_seqs.

            # if a file is created during the process delete it

            # 4. Export model and ref_genome to model_dir
            self.write_mdldata_dir()
            # path_mdl_data = self.igor_wd + "/" + self.igor_batchname + "_mdldata"
            self.update_model_filenames(igor_model_dir_path=self.igor_mdldata_dir)
            self.update_ref_genome(igor_model_dir_path=self.igor_mdldata_dir)
            self.update_batch_filenames()
            # self.mdl.write_mdldata_dir(path_mdl_data)

            # 5. Run infer model
            if input_sequences is None:
                # use the self.igor_read_seqs
                # Do not delete self.igor_read_seqs
                self._run_infer(igor_read_seqs=self.igor_read_seqs)
            else:
                # write a temporary file
                tmp_igor_read_seqs = self.igor_read_seqs
                fln_input_sequences = self.igor_wd + "/" + self.igor_batchname + "_input_sequences.csv"
                write_sequences_to_file(input_sequences, fln_input_sequences)
                self._run_infer(igor_read_seqs=fln_input_sequences)
                import os
                os.unlink(fln_input_sequences)
                self.igor_read_seqs = tmp_igor_read_seqs

            self.load_IgorModel_from_infer_files()
        except Exception as e:
            raise e
        else:
            return self.mdl
        finally:
            self._run_clean_batch_infer()



    def _run_infer(self, igor_read_seqs=None,
                   igor_model_parms_file=None,
                   igor_model_marginals_file=None,
                   fln_V_gene_CDR3_anchors=None,
                   fln_J_gene_CDR3_anchors=None,
                   igor_fln_db=None,
                   igor_db=None,
                   mdl:Union[IgorModel, None]=None, N_iter=None,
                   igor_infer_dict_opts=None,
                   igor_align_dict_opts=None):
        """Run inference and return IgorModel object"""
        # "igor -set_wd $WDPATH -batch foo -species human -chain beta
        # -evaluate -output --scenarios 10"

        try:
            if igor_align_dict_opts is not None:
                self.igor_align_dict_options = igor_align_dict_opts

            if igor_infer_dict_opts is not None:
                self.igor_infer_dict_options = igor_infer_dict_opts

            if N_iter is not None:
                if self.igor_infer_dict_options is None:
                    self.igor_infer_dict_options = copy(igor_infer_dict_options)
                self.igor_infer_dict_options['--N_iter']['active'] = True
                self.igor_infer_dict_options['--N_iter']['value'] = str(N_iter)

            if igor_read_seqs is not None:
                self.igor_read_seqs = igor_read_seqs

            if mdl is not None:
                self.mdl = mdl

            if self.mdl is None:
                try:
                    # 1. Load from default from igor_species and igor_chain
                    self.load_IgorModel(igor_model_parms_file=igor_model_parms_file,
                       igor_model_marginals_file=igor_model_marginals_file,
                       fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                       fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors)
                except:
                    try:
                        self.load_mdl_from_db(igor_fln_db=igor_fln_db, igor_db=igor_db)
                    except:
                        pass
                    pass

            import pathlib
            pathlib.Path(self.igor_wd).mkdir(parents=True, exist_ok=True)

            if self.b_align is False:
                self.run_align(igor_read_seqs=igor_read_seqs, igor_align_dict_opts=igor_align_dict_opts)
                print("== Alignment finished! ==")

            if self.mdl is not None:
                self._update_mdldata_batch_filenames()
                self.write_mdldata_dir(self.igor_mdldata_dir)
                # if self.igor_mdldata_dir is not None:
                #     self.update_model_filenames(self.igor_mdldata_dir)
                #     self.update_ref_genome()

                cmd = self.igor_exec_path
                cmd = cmd + " -set_wd " + self.igor_wd
                cmd = cmd + " -batch " + self.igor_batchname
                cmd = cmd + " -set_custom_model " + self.igor_fln_mdldata_parms + " " + self.igor_fln_mdldata_marginals

                if os.path.isfile(self.igor_fln_mdldata_V_gene_CDR3_anchors) or \
                    os.path.isfile(self.igor_fln_mdldata_J_gene_CDR3_anchors):
                    cmd = cmd + " -set_CDR3_anchors "
                    if os.path.isfile(self.igor_fln_mdldata_V_gene_CDR3_anchors):
                        cmd = cmd + " --V " + self.igor_fln_mdldata_V_gene_CDR3_anchors
                    if os.path.isfile(self.igor_fln_mdldata_J_gene_CDR3_anchors):
                        cmd = cmd + " --J " + self.igor_fln_mdldata_J_gene_CDR3_anchors
                # here the evaluation
                cmd = cmd + " -infer "
                cmd = cmd + " " + command_from_dict_options(self.igor_infer_dict_options)

            else:
                cmd = self.igor_exec_path
                cmd = cmd + " -set_wd " + self.igor_wd
                cmd = cmd + " -batch " + self.igor_batchname
                cmd = cmd + " -set_custom_model " + self.igor_model_parms_file + " " + self.igor_model_marginals_file
                if (self.fln_V_gene_CDR3_anchors is not None) or \
                        (self.fln_V_gene_CDR3_anchors is not None):
                    cmd = cmd + " -set_CDR3_anchors "
                    if self.fln_V_gene_CDR3_anchors is not None:
                        cmd = cmd + " --V " + self.fln_V_gene_CDR3_anchors
                    if self.fln_J_gene_CDR3_anchors is not None:
                        cmd = cmd + " --J " + self.fln_J_gene_CDR3_anchors
                # here the evaluation
                cmd = cmd + " -infer "
                cmd = cmd + " " + command_from_dict_options(self.igor_infer_dict_options)


            print(cmd)
            output = run_command_print(cmd)

            # run_command_no_output(cmd)
            self.b_infer = True  # FIXME: If run_command success then True
            self._update_infer_batch_filenames()
            self.load_IgorModel_from_infer_files()
        except Exception as e:
            raise e
        else:
            return self.mdl
            # return output

        # finally:
        #     self.run_clean_batch()

    def _run_generate(self, N_seqs:int=1, mdl:Union[None,IgorModel]=None, seed=None, igor_wd=None, igor_batchname=None,
                      igor_model_parms_file=None, igor_model_marginals_file=None,
                      fln_V_gene_CDR3_anchors=None, fln_J_gene_CDR3_anchors=None,
                      igor_db=None, igor_fln_db=None,
                      igor_species=None, igor_chain=None, #return_df=False,
                      fln_output_prefix:Union[None, str]=None, clean_batch=False):
        """
        Run IGoR generate command line.
        :param N_seqs: Integer number of sequences to generate (default 1)
        :param mdl: IgorModel object to generate sequences if None is provide it uses the self.mdl (default None)
        :param seed: Seed to generate random sequences (default None). If None then IGoR's chooses a random seed.
        :param igor_wd: Working directory to execuate IGoR, if None it uses self.igor_wd which default value is the current diretory(default None).
        :param igor_batchname: IGoR batch option to identify the IGoR's execution(default None).
        :param igor_model_parms_file: If no model is specified in self.mdl or mdl option, with this option a model_parms path can be used (default None).
        :param igor_model_marginals_file: Same as igor_model_parms_file (default None).
        :param fln_V_gene_CDR3_anchors: CDR3 IGoR's anchors V path to file (default None).
        :param fln_J_gene_CDR3_anchors: CDR3 IGoR's anchors J path to file (default None).
        :param igor_db: A database within a IGoR model can be used to get the model (default None).
        :param igor_fln_db: A database file can be used to generate sequences (default None).
        :param igor_species: IGoR's name of species (default None).
        :param igor_chain: IGoR's name of chain (default None).
        """

        try:
            if seed is not None:
                self.igor_generate_dict_options['--seed']['active'] = True
                self.igor_generate_dict_options['--seed']['value'] = str(seed)

            if mdl is not None:
                self.mdl = mdl

            if igor_species is not None:
                self.igor_species = igor_species

            if igor_chain is not None:
                self.igor_chain = igor_chain

            if igor_wd is not None:
                self.igor_wd = igor_wd

            if igor_batchname is not None:
                self.igor_batchname = igor_batchname

            if igor_model_parms_file is not None:
                self.igor_model_parms_file = igor_model_parms_file

            if igor_model_marginals_file is not None:
                self.igor_model_marginals_file = igor_model_marginals_file

            if igor_fln_db is not None:
                self.igor_fln_db = igor_fln_db

            if igor_db is not None:
                self.igor_db = igor_db

            if self.mdl is None:
                try:
                    # 1. Load from default from igor_species and igor_chain
                    self.load_IgorModel(igor_model_parms_file=igor_model_parms_file,
                       igor_model_marginals_file=igor_model_marginals_file,
                       fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                       fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors)
                except:
                    try:
                        self.load_mdl_from_db(igor_fln_db=igor_fln_db, igor_db=igor_db)
                    except:
                        pass
                    pass
            import pathlib
            pathlib.Path(self.igor_wd).mkdir(parents=True, exist_ok=True)
            self._update_mdldata_batch_filenames()
            self.write_mdldata_dir(self.igor_mdldata_dir)

            self.update_model_filenames(self.igor_mdldata_dir)
            self.update_ref_genome()


            cmd = self.igor_exec_path
            cmd = cmd + " -set_wd " + self.igor_wd
            cmd = cmd + " -batch " + self.igor_batchname
            cmd = cmd + " -set_custom_model " + self.igor_model_parms_file + " " + self.igor_model_marginals_file
            if (self.fln_V_gene_CDR3_anchors is not None) or \
                    (self.fln_V_gene_CDR3_anchors is not None):
                cmd = cmd + " -set_CDR3_anchors "
                if self.fln_V_gene_CDR3_anchors is not None:
                    cmd = cmd + " --V " + self.fln_V_gene_CDR3_anchors
                if self.fln_J_gene_CDR3_anchors is not None:
                    cmd = cmd + " --J " + self.fln_J_gene_CDR3_anchors
            # if N_seqs is not None:
            cmd = cmd + " -generate " + str(N_seqs)
            cmd = cmd + " " + command_from_dict_options(self.igor_generate_dict_options)
            # else:
            #     cmd = cmd + " -generate "
            print(cmd)

            run_command(cmd)
            # run_command_print(cmd)

            self._update_generate_batch_filenames()
            # path_generated = self.igor_wd + "/" + self.igor_batchname + "_generated/"
            # self.igor_fln_generated_realizations_werr = path_generated + "generated_realizations_werr.csv"
            # self.igor_fln_generated_seqs_werr = path_generated + "generated_seqs_werr.csv"
            # self.igor_fln_generation_info = path_generated + "generated_seqs_werr.out"
            self.b_generate = True

            # FIXME: LOAD TO DATABASE CREATE PROPER TABLES FOR THIS
            # import pandas as pd
            if fln_output_prefix is not None:
                import os
                self.igor_fln_generated_seqs_werr = self.igor_wd + "/" + self.igor_batchname + "_generated/generated_seqs_werr.csv"
                self.igor_fln_generated_realizations_werr = self.igor_wd + "/" + self.igor_batchname + "_generated/generated_realizations_werr.csv"
                self.igor_fln_generation_info = self.igor_wd + "/" + self.igor_batchname + "_generated/generation_info.out"
                output_generated_sequences = fln_output_prefix + "_sequences.csv"
                output_generated_realizations = fln_output_prefix + "_realizations.csv"
                output_generated_info = fln_output_prefix + "_info.out"
                # output_generated_sequences_airr = fln_output_prefix + "_sequences.tsv"

                os.rename(self.igor_fln_generated_seqs_werr, output_generated_sequences)
                os.rename(self.igor_fln_generated_realizations_werr, output_generated_realizations)
                os.rename(self.igor_fln_generation_info, output_generated_info)
                # TODO: IF DIRECTORY EMPTY DELETE IT.

                self.igor_fln_generated_seqs_werr = output_generated_sequences
                self.igor_fln_generated_realizations_werr = output_generated_realizations
                self.igor_fln_generation_info = output_generated_info


        except Exception as e:
            raise e
        else:
            return self.igor_fln_generated_seqs_werr, self.igor_fln_generated_realizations_werr, self.igor_fln_generation_info
        finally:
            if clean_batch:
                self._run_clean_batch_generate()
                self._run_clean_batch_mdldata()

    def run_generate_to_dataframe(self, N):
        self._run_generate(self, N)

        # FIXME: LOAD TO DATABASE CREATE PROPER TABLES FOR THIS
        # import pandas as pd
        df = pd.read_csv(self.igor_fln_generated_seqs_werr, delimiter=';').set_index('seq_index')
        return df

    def run_clean_batch(self):
        """Clean all files defined with batchname and igor_wd"""
        try:
            self._run_clean_batch_evaluate()
            self._run_clean_batch_aligns()
            self._run_clean_batch_infer()
            self._run_clean_batch_generate()
        except Exception as e:
            print(e)
            pass

        # cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_evaluate"
        # run_command_no_output(cmd)
        # cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_output"
        # run_command_no_output(cmd)
        # cmd = "rm " + self.igor_wd + "/aligns/" + self.igor_batchname + "*.csv"
        # run_command_no_output(cmd)
        # cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_generated"
        # run_command_no_output(cmd)

    def _run_clean_batch_aligns(self):
        try:
            cmd = "rm " + self.igor_wd + "/aligns/" + self.igor_batchname + "*.csv"
            run_command_no_output(cmd)
            cmd = "rm " + self.igor_wd + "/aligns/aligns_info.out"
            run_command_no_output(cmd)
            cmd = "rmdir --ignore-fail-on-non-empty " + self.igor_wd + "/aligns"
            run_command_no_output(cmd)
        except Exception as e:
            raise e

    def _run_clean_batch_mdldata(self):
        try:
            cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_mdldata"
            run_command_no_output(cmd)
        except Exception as e:
            raise e

    def _run_clean_batch_infer(self):
        try:
            cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_inference"
            run_command_no_output(cmd)
            cmd = "rmdir --ignore-fail-on-non-empty " + self.igor_wd + "/" + self.igor_batchname + "_output"
            run_command_no_output(cmd)
        except Exception as e:
            raise e

    def _run_clean_batch_generate(self):
        try:
            cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_generated"
            run_command_no_output(cmd)
        except Exception as e:
            raise e

    def _run_clean_batch_output(self):
        try:
            cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_output"
            run_command_no_output(cmd)
        except Exception as e:
            raise e

    def _run_clean_batch_evaluate(self):
        try:
            cmd = "rm -r " + self.igor_wd + "/" + self.igor_batchname + "_evaluate"
            run_command_no_output(cmd)
            try:
                self._run_clean_batch_output()
            except:
                pass
        except Exception as e:
            raise e


    def create_db(self, igor_fln_db=None):
        if igor_fln_db is not None:
            self.igor_fln_db = igor_fln_db
        if self.igor_fln_db is None:
            # Generate it using batchname
            self.igor_fln_db = self.igor_batchname + ".db"
        self.igor_db = IgorSqliteDB.create_db(self.igor_fln_db)

    def load_db(self):
        """
        Save all batch and models files in a single sqlite database.
        return: None
        """
        try:
            self.load_db_from_indexed_sequences()
        except Exception as e:
            print(e)
            pass
        try:
            self.load_db_from_indexed_cdr3()
        except Exception as e:
            print(e)
            pass
        try:
            self.load_db_from_genomes()
        except Exception as e:
            print(e)
            pass
        try:
            self.load_db_from_alignments()
        except Exception as e:
            print(e)
            pass
        try:
            self.load_IgorModel()
        except Exception as e:
            print(e)
            pass
        try:
            self.load_db_from_models()
        except Exception as e:
            print(e)
            pass
        try:
            self.load_db_from_bestscenarios()
        except Exception as e:
            print(e)
            pass
        try:
            self.load_db_from_pgen()
        except Exception as e:
            print(e)
            pass


    def load_db_from_indexed_sequences(self, igor_fln_indexed_sequences=None):
        """
        Load indexed_sequences from database
        :param igor_fln_indexed_sequences: csv file of indexed_sequences, generated with -read_seqs
        """
        if igor_fln_indexed_sequences is not None:
            self.igor_fln_indexed_sequences = igor_fln_indexed_sequences
        self.igor_db.load_IgorIndexedSeq_FromCSV(self.igor_fln_indexed_sequences)

    # load genome templates from fasta and csv files.
    def load_db_from_genomes(self):
        print("Loading Gene templates ...")
        try:
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
            self.load_db_from_anchors()
            # self.igor_db.load_IgorGeneAnchors_FromCSV("V", self.genomes.fln_V_gene_CDR3_anchors)
            # self.igor_db.load_IgorGeneAnchors_FromCSV("J", self.genomes.fln_J_gene_CDR3_anchors)
            # except Exception as e:
            #     print("ERROR : ", e)
        except Exception as e:
            raise e

    def load_db_from_anchors(self):
        """Load anchors from database"""
        self.igor_db.load_IgorGeneAnchors_FromCSV("V", self.genomes.fln_V_gene_CDR3_anchors)
        self.igor_db.load_IgorGeneAnchors_FromCSV("J", self.genomes.fln_J_gene_CDR3_anchors)


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
        print("Alignments loaded in database in " + str(self.igor_fln_db))

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

    def load_db_from_bestscenarios(self, igor_fln_output_scenarios:Union[None, str]=None,
                                   mdl:Union[None, IgorModel]=None):
        if igor_fln_output_scenarios is not None:
            self.igor_fln_output_scenarios = igor_fln_output_scenarios
        if mdl is not None:
            self.mdl = mdl
        print(self.igor_fln_output_scenarios)
        self.igor_db.load_IgorBestScenarios_FromCSV(self.igor_fln_output_scenarios, self.mdl)



    def load_db_from_pgen(self):
        print(self.igor_fln_output_pgen)
        self.igor_db.load_IgorPgen_FromCSV(self.igor_fln_output_pgen)

    def load_mdl_from_db(self, igor_fln_db: Union[str, None] = None, igor_db: Union[IgorSqliteDB, None] = None):
        """
        Return a IgorModel object in self.mdl from igor_fln_db or igor_db.
        """

        if igor_db is not None:
            self.igor_db = igor_db

        if igor_fln_db is not None:
            self.igor_fln_db = igor_fln_db

        try:
            if self.igor_db is None:
                if self.igor_fln_db is not None:
                    self.create_db(self.igor_fln_db)
            self.mdl = self.igor_db.get_IgorModel()
            print("Model loaded from database")
        except Exception as e:
            e_message = "WARNING: Igor Model was not found in " + str(self.igor_fln_db)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

            # pass
        # return self.mdl

    # def get_IgorModel_from_db(self):
    #     self.mdl = self.igor_db.get_IgorModel()
    #     return self.mdl

    # FIXME: this method should be deprecated!!!
    def load_VDJ_database(self, flnIgorSQL):
        self.flnIgorSQL = flnIgorSQL
        self.igor_db = IgorSqliteDB(flnIgorSQL)
        # FIXME :EVERYTHING
        flnIgorIndexedSeq = self.igor_wd + "/aligns/" + self.igor_batchname + "_indexed_sequences.csv"
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

        flnIgorDB = self.igor_batchname + ".db"
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
        flnIgorBestScenarios = self.igor_wd + "/" + self.igor_batchname + "_output/best_scenarios_counts.csv"
        self.igor_db_bs = IgorSqliteDBBestScenariosVDJ(flnIgorBSSQL)  # IgorDBBestScenariosVDJ.sql
        self.igor_db_bs.createSqliteDB(self.igor_batchname + "_bs.db")
        self.igor_db_bs.load_IgorBestScenariosVDJ_FromCSV(flnIgorBestScenarios)

    def get_pgen_pd(self):
        # load pgen file
        import pandas as pd
        # df = pd.read_csv(self.igor_fln_output_pgen, sep=';')
        # df = df.set_index('seq_index')
        # df = df.sort_index()
        df = get_df_pgen(self.igor_fln_output_pgen)
        df_seq = pd.read_csv(self.igor_fln_indexed_sequences, sep=';')
        df_seq = df_seq.set_index('seq_index').sort_index()
        df_cdr3 = pd.read_csv(self.igor_fln_indexed_CDR3, sep=';')
        df_cdr3 = df_cdr3.set_index('seq_index').sort_index()
        df = df.merge(df_seq, left_index=True, right_index=True)
        df = df.merge(df_cdr3, left_index=True, right_index=True)
        return df

    #### DATABASE METHODS
    def db_ls(self):
        """List igor_db tables"""
        try:
            self.igor_db.list_from_db()
        except Exception as e:
            raise e

    def db_get_naive_align_dict_by_seq_index(self, seq_index):
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

    def db_str_fasta_naive_align_by_seq_index(self, seq_index):
        """ Given an Sequence index and the corresponding alignments vj/ vdj
            return a string with considering only offset"""

        fasta_list = list()
        indexed_sequence, vdj_alignments_dict = self.db_get_naive_align_dict_by_seq_index(seq_index)
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
            str_fasta_description2 = "> " + vdj_alignments_dict[key].strGene_name + ", score : " + str(
                vdj_alignments_dict[key].score)
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
        str_fasta = '\n'.join([fasta_rec[0] + "\n" + fasta_rec[1] for fasta_rec in fasta_list])
        return str_fasta  # , fasta_list

    def db_plot_naive_align_by_seq_index(self, seq_index):
        import Bio.AlignIO
        import io
        aaa = self.db_str_fasta_naive_align_by_seq_index(seq_index)
        aln = Bio.AlignIO.read(io.StringIO(aaa), 'fasta')
        view_alignment(aln)

    def db_export_to_igorfiles(self):
        print("Export: ")
        # --- 1. Indexed Sequences
        if self.igor_db.Q_sequences_in_db() and not (self.igor_fln_indexed_sequences is None):
            try:
                self.igor_db.write_IgorIndexedSeq_to_CSV(self.igor_fln_indexed_sequences)
            except Exception as e:
                print("ERROR: write_IgorIndexedSeq_to_CSV", e)
        else:
            print("No IgorIndexedSeq Table not exported")

        # --- 2. Gene Templates
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

        # --- 3. Alignments
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
                    print("ERROR: write_IgorModel_to_TXT ", e)
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

    def db_export_IgorIndexedSeq(self,
                                 igor_fln_indexed_sequences: Union[None, str] = None):
        """
        Export from database IGoR's indexed_seq files
        :param igor_fln_indexed_sequences: Path of csv file to save IgorIndexedSeq
        """
        try:
            if igor_fln_indexed_sequences is not None:
                self.igor_fln_indexed_sequences = igor_fln_indexed_sequences
            self.igor_db.write_IgorIndexedSeq_to_CSV(self.igor_fln_indexed_sequences)
        except Exception as e:
            e_message = "IgorTask.export_from_db_IgorIndexedSeq : igor_fln_db " + str(self.igor_fln_db)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    # FIXME: in dev
    def db_export_IgorGenomes(self, igor_path_ref_genome: Union[None, str] = None,
                              igor_model_dir_path: Union[None, str] = None,
                              fln_genomicVs: Union[None, str] = None,
                              fln_genomicDs: Union[None, str] = None,
                              fln_genomicJs: Union[None, str] = None,
                              fln_V_gene_CDR3_anchors: Union[None, str] = None,
                              fln_J_gene_CDR3_anchors: Union[None, str] = None):
        """
        Export from database IGoR's indexed_seq files
        :param igor_fln_indexed_sequences: Path of csv file to save IgorIndexedSeq
        :param fln_genomicVs: Path of fasta file with genomic V templates, (default None),
        :param fln_genomicDs: Path of fasta file with genomic D templates, (default None),
        :param fln_genomicJs: Path of fasta file with genomic J templates, (default None),
        :param fln_V_gene_CDR3_anchors: Path of csv file with genomic V templates, (default None),
        :param fln_J_gene_CDR3_anchors: Path of csv file with genomic J templates, (default None)
        """
        try:
            # Assign filenames to export data
            self.update_ref_genome(igor_path_ref_genome=igor_path_ref_genome,
                                   igor_model_dir_path=igor_model_dir_path,
                                   fln_genomicVs=fln_genomicVs, fln_genomicDs=fln_genomicDs,
                                   fln_genomicJs=fln_genomicJs,
                                   fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors,
                                   fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors)

            self.igor_db.write_IgorGeneTemplate_to_fasta("V", self.fln_genomicVs)
            self.igor_db.write_IgorGeneTemplate_to_fasta("J", self.fln_genomicJs)
            try:
                self.igor_db.write_IgorGeneTemplate_to_fasta("D", self.fln_genomicDs)
            except Exception as e:
                print("D genes could not be exported")
                pass

            # TODO: ADD SUPPORT FOR sep=',' (OLGA)
            try:
                self.igor_db.write_IgorGeneAnchors_to_CSV("V", self.fln_V_gene_CDR3_anchors)
            except Exception as e:
                print(e)
                pass
            try:
                self.igor_db.write_IgorGeneAnchors_to_CSV("J", self.fln_J_gene_CDR3_anchors)
            except Exception as e:
                print(e)
                pass

        except Exception as e:
            e_message = "IgorTask.export_from_db_IgorGenomes : igor_fln_db " + str(self.igor_fln_db)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def write_mdldata_dir(self, igor_mdldata_dir:Union[str, None, Path] = None,
                          mdl:Union[IgorModel, None] = None):
        if igor_mdldata_dir is not None:
            self.igor_mdldata_dir = igor_mdldata_dir

        if self.igor_mdldata_dir is None:
            self.igor_mdldata_dir = self.igor_wd + "/" + self.igor_batchname + "_mdldata"

        if self.igor_path_ref_genome is None:
            self.igor_path_ref_genome = self.igor_mdldata_dir + "/ref_genome/"
        if mdl is not None:
            self.mdl = copy.deepcopy(mdl)

        self.mdl.write_mdldata_dir(self.igor_mdldata_dir)

    #### AIRR methods ###
    def parse_scenarios_to_airr(self, igor_fln_output_scenarios, airr_fln_output_scenarios):
        # 1. Read header of and make a list
        # open(igor_fln_output_scenarios)
        # 2.
        pass

    def get_dataframe_scenarios(self, igor_fln_output_scenarios:Union[None, str]=None,
                                                           mdl:Union[None, IgorModel]=None):
        try:
            if igor_fln_output_scenarios is None:
                igor_fln_output_scenarios = self.igor_fln_output_scenarios

            if mdl is not None:
                self.mdl = mdl

            return self.mdl.get_dataframe_from_fln_generated_realizations_werr(
                self.igor_fln_output_scenarios)

        except Exception as e:
            raise e


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




#####################################################################################

def naive_align(input_sequences: Union[None, str, Path, pd.DataFrame, np.array, list] = None,
              mdl: Union[None, str, Path, IgorModel, IgorModel_Parms] = None,
              igor_wd=None, igor_threads=None, igor_batchname=None, batch_clean=True,
              igor_align_dict_opts:Union[None, dict]=None):
    """Return preprocessing with IGoR
    :param input_sequences: Sequences to align as dataframe, numpy array, list or string.
    :param mdl: IgorModel
    :param igor_wd: IGoR working directory
    :param igor_threads: IGoR number of threads option
    :param igor_align_dict_opts: Additional options for naive alignment.
    """
    tmp_dir = tempfile.TemporaryDirectory(prefix='igor_aligning_', dir='.')

    df_alignments_dict = dict()
    if igor_wd is None:
        igor_wd = tmp_dir.name
        batch_clean = True
    else:
        os.system("mkdir -p " + igor_wd)
        batch_clean = False

    # 2. Create an IgorTask
    import copy
    mdl_copy = copy.deepcopy(mdl)
    task = IgorTask(igor_wd=igor_wd, mdl=mdl_copy)
    task.igor_fln_indexed_functionality = igor_wd + "/aligns/" + task.igor_batchname + "_indexed_functionality.csv"

    if igor_align_dict_opts is None:
        task.igor_align_dict_options = copy.deepcopy(igor_align_dict_options)
        task.igor_align_dict_options['--all']['dict_options']['---best_align_only']['active'] = True
        task.igor_align_dict_options['--all']['dict_options']['---best_align_only']['value'] = "true"


    fln_input_sequences = igor_wd + "/" + task.igor_batchname + "input_sequences.csv"

    # 3. Write Sequences in file if file not exist
    write_sequences_to_file(input_sequences, fln_input_sequences)

    # 4. Export model and ref_genome to model_dir
    path_mdl_data = task.igor_wd + "/" + task.igor_batchname + "_mdldata"
    task.update_model_filenames(igor_model_dir_path=path_mdl_data)
    task.update_ref_genome(igor_model_dir_path=path_mdl_data)
    task.update_batch_filenames()
    task.mdl.write_mdldata_dir(path_mdl_data)

    df_indexed_functionality = None
    df_indexed_CDR3 = None
    # 5. Run infer model
    try:
        task.run_align(igor_read_seqs=fln_input_sequences, igor_threads=igor_threads, igor_align_dict_opts=igor_align_dict_opts)
        # if os.path.isfile(self.igor_fln_align_V_alignments):
        #     df_V_aligments = pd.read_csv(self.igor_fln_align_V_alignments, sep=';')
        #     df_alignments_dict['V'] = df_V_aligments
        # if os.path.isfile(self.igor_fln_align_D_alignments):
        #     df_D_aligments = pd.read_csv(self.igor_fln_align_D_alignments, sep=';')
        #     df_alignments_dict['D'] = df_D_aligments
        # if os.path.isfile(self.igor_fln_align_J_alignments):
        #     df_J_aligments = pd.read_csv(self.igor_fln_align_J_alignments, sep=';')
        #     df_alignments_dict['J'] = df_J_aligments
        if os.path.isfile(task.igor_fln_indexed_functionality):
            df_indexed_functionality = pd.read_csv(task.igor_fln_indexed_functionality, sep=';', index_col='seq_index')
            df_indexed_functionality['functionality'] = df_indexed_functionality['functionality'].apply(
                lambda x: False if (x == 0) else True)
            df_indexed_functionality['in_frame'] = df_indexed_functionality['in_frame'].apply(
                lambda x: False if (x == 0) else True)
            df_indexed_functionality['stop_codon'] = df_indexed_functionality['stop_codon'].apply(
                lambda x: False if (x == 0) else True)
            #return df_indexed_functionality
        if os.path.isfile(task.igor_fln_indexed_CDR3):
            df_indexed_CDR3 = pd.read_csv(task.igor_fln_indexed_CDR3, sep=';', index_col='seq_index')

        return df_indexed_functionality, df_indexed_CDR3
    except Exception as e:
        raise e
    finally:
        tmp_dir.cleanup()

def generate(Nseqs, mdl:IgorModel, igor_wd=None, igor_batchname=None,
             seed=None,
             clean_batch=True, return_scenarios=False):
    """Return pandas dataframe with generated sequences Only sequences, not scenarios
    :param Nseqs: Number of sequences to be generated.
    :param mdl: IgorModel instance to generate sequences.
    :param igor_wd: Working directory use to run IGoR to generate sequences.
    :param seed: Seed to generate sequences.
    :param batch_clean: If True clean all files used with IGoR, if false keep all files, default True.
    :param return_scenarios: Flag if True return scenarios dataframe.
    """
    try:
        tmp_generate_dir = tempfile.TemporaryDirectory(prefix='igor_generating_', dir='.')
        if igor_wd is None:
            igor_wd = tmp_generate_dir.name

        task = IgorTask(mdl=mdl, igor_wd=igor_wd, igor_batchname=igor_batchname)

        if seed is not None:
            task.igor_generate_dict_options['--seed']['active'] = True
            task.igor_generate_dict_options['--seed']['value'] = str(seed)

        pd_sequences = task.generate(N_seqs=Nseqs, return_scenarios=return_scenarios)

    except Exception as e:
        raise e
    else:
        return pd_sequences
    finally:
        if clean_batch:
            tmp_generate_dir.cleanup()

    # try:
    #     # TODO:
    #     # 1. Use a IgorModel to create an IgorTask
    #     # 2. Create a temporary directory
    #     # 3. generate sequences return_df = True
    #
    #     with tempfile.TemporaryDirectory(prefix='igor_generating_', dir='.') as tmp_generate_dirname:
    #         task = IgorTask(mdl=mdl, igor_wd=tmp_generate_dirname)
    #
    #         task.gen_random_batchname()
    #
    #         task.update_batch_filenames()
    #         path_mdl_data = tmp_generate_dirname + "/" + task.igor_batchname + "_mdldata"
    #         task.update_model_filenames(igor_model_dir_path=path_mdl_data)
    #         task.update_ref_genome()
    #         task.mdl.write_model(task.igor_model_parms_file, task.igor_model_marginals_file)
    #         mdl_inferrred = task._run_generate(Nseqs, mdl=mdl, return_df=True)
    #
    #         pd_sequences = task._run_generate(Nseqs, return_df=True)
    #
    # except Exception as e:
    #     raise e
    # else:
    #     return pd_sequences #mdl_inferrred
    # finally:
    #     task.run_clean_batch()


def infer(input_sequences:Union[str, list, pd.DataFrame, np.ndarray, Path],
          mdl:IgorModel, igor_wd=None, batch_clean=True, return_likelihoods=True,
          N_iter=None,
          igor_infer_dict_opts:Union[None, dict]=None,
          igor_align_dict_opts:Union[None, dict]=None)->IgorModel:
    """
    Returns inferred IgorModel from input_sequences with initial model mdl.
    :param input_sequences: String sequences in list, numpy array, pandas dataframe or file.
    :param mdl: IgorModel object.
    :param igor_wd: Working directory use to run IGoR to infer.
    :param batch_clean: If True clean all files used with IGoR, if false keep all files, default True.
    :param return_likelihoods: If True return likelihoods pandas dataframe with likelihoods
    :param N_iter: Number of iterations to infer model.
    :param igor_infer_dict_opts: Additional IGoR options to pass for inference.
    :param igor_align_dict_opts: Additional IGoR options to pass for alignment.
    """
    try:
        import tempfile
        # batch_clean = False
        # 1. Create a temporary directory igor_wd=tmp_dir.name
        tmp_dir = tempfile.TemporaryDirectory(prefix='igor_inferring_', dir='.')
        if igor_wd is None:
            igor_wd = tmp_dir.name
            batch_clean = True
        else:
            os.system("mkdir -p " + igor_wd)
            batch_clean = False

        # 2. Create an IgorTask
        import copy
        mdl_copy = copy.deepcopy(mdl)
        task = IgorTask(igor_wd=igor_wd, mdl=mdl_copy)
        fln_input_sequences = igor_wd + "/" + task.igor_batchname + "input_sequences.csv"

        # 3. Write Sequences in file if file not exist
        write_sequences_to_file(input_sequences, fln_input_sequences)

        # 4. Export model and ref_genome to model_dir
        path_mdl_data = task.igor_wd + "/" + task.igor_batchname + "_mdldata"
        task.update_model_filenames(igor_model_dir_path=path_mdl_data)
        task.update_ref_genome(igor_model_dir_path=path_mdl_data)
        task.update_batch_filenames()
        task.mdl.write_mdldata_dir(path_mdl_data)
        # print(task)

        # 5. Run infer model
        task._run_infer(igor_read_seqs=fln_input_sequences, N_iter=N_iter, igor_infer_dict_opts=igor_infer_dict_opts, igor_align_dict_opts=igor_align_dict_opts)
        task.load_IgorModel_from_infer_files()

        """
        path_mdl_data = task.igor_wd + "/" + task.igor_batchname + "_mdldata"
        task.update_model_filenames(igor_model_dir_path=path_mdl_data)
        task.update_ref_genome()
        task.mdl.write_model(task.igor_model_parms_file, task.igor_model_marginals_file)
        mdl_inferrred = task.run_generate(Nseqs, mdl=mdl, return_df=True)
        """

    except Exception as e:
        raise e
    else:
        if return_likelihoods == True:
            return task.mdl, task.df_infer_likelihoods
        else:
            return task.mdl
    finally:
        if batch_clean:
            task.run_clean_batch()
        tmp_dir.cleanup()


def evaluate(input_sequences:Union[str, pd.DataFrame, np.ndarray, list, tuple, Path],
             mdl:IgorModel, N_scenarios=None, igor_wd:Union[None, str, Path]=None,
             airr_format=False, batch_clean=True, use_db=False, fln_output:Union[None, str, Path]=None, b_V_offset=False,
             igor_evaluate_dict_opts:Union[None, dict]=None,
             igor_align_dict_opts:Union[None, dict]=None):
    """
    Evaluate input sequences with provided model
    :param input_sequences:Union[str, pd.DataFrame, np.ndarray, Path]
    :param mdl:IgorModel
    :param N_scenarios: Number of scenarios requested, default returns 1.
    :param airr_format: True return pandas dataframe as AIRR format
    :param fln_output: Output filename to save scenarios
    :param b_V_offset: Flag to return V_offset alignments as pandas dataframe.
    :param igor_wd: working directory
    :param batch_clean: Remove all temporary files True by default.
    :param use_db: Flag to save data in sqlite database.
    :param igor_evaluate_dict_opts: Additional IGoR options for evaluation.
    :param igor_align_dict_opts: Additional IGoR options for aligns.
    """

    # Run evaluate

    import tempfile
    pd_airr_rearrangement = None
    df_V_offsets = None
    try:
        # batch_clean = False
        # 1. Create a temporary directory igor_wd=tmp_dir.name
        tmp_dir = tempfile.TemporaryDirectory(prefix='igor_evaluating_', dir='.')
        # if igor_wd is set then use that directory, create it if doesn't exist, but
        # if igor_wd is None then use the temporary directory.
        if igor_wd is None:
            igor_wd = tmp_dir.name
            batch_clean = True
        else:
            os.system("mkdir -p " + igor_wd)
            batch_clean = False

        # 2. Create an IgorTask
        import copy
        mdl_copy = copy.deepcopy(mdl)
        task = IgorTask(igor_wd=igor_wd, mdl=mdl_copy)
        fln_input_sequences = igor_wd + "/" + task.igor_batchname + "input_sequences.csv"

        # 3. Write Sequences in file if file not exist
        write_sequences_to_file(input_sequences, fln_input_sequences)

        # 4. Export model and ref_genome to model_dir
        path_mdl_data = task.igor_wd + "/" + task.igor_batchname + "_mdldata"
        task.update_model_filenames(igor_model_dir_path=path_mdl_data)
        task.update_ref_genome(igor_model_dir_path=path_mdl_data)
        task.update_batch_filenames()
        task.mdl.write_mdldata_dir(path_mdl_data)
        # print(task)

        # 5. Run evaluate model
        task._run_evaluate(igor_read_seqs=fln_input_sequences, N_scenarios=N_scenarios,
                           igor_evaluate_dict_opts=igor_evaluate_dict_opts,
                           igor_align_dict_opts=igor_align_dict_opts)

        if airr_format:
            # Save evaluations in database
            try:
                if use_db:
                    task.create_db()
                    task.load_db()
                    # task.load_db_from_indexed_sequences()
                    # task.load_db_from_indexed_cdr3()
                    # task.load_db_from_genomes()
                    # task.load_db_from_alignments()
                    # task.load_IgorModel()
                    # task.load_db_from_models()
                    # task.load_db_from_bestscenarios()
                    # task.load_db_from_pgen()

                    base_fln_output = task.igor_fln_db.split(".db")[0]
                    output_fln_prefix = base_fln_output
                    output_fln_airr = output_fln_prefix + ".tsv"
                    task.igor_db.export_IgorBestScenarios_to_AIRR(output_fln_airr)
                    pd_airr_rearrangement = pd.read_csv(output_fln_airr, sep='\t')
                else:
                    pd_airr_rearrangement = task.mdl.get_dataframe_from_fln_generated_realizations_werr(
                        task.igor_fln_output_scenarios)
            except Exception as e:
                raise e

        else:
            try:
                pd_airr_rearrangement = task.mdl.get_dataframe_from_fln_generated_realizations_werr(
                    task.igor_fln_output_scenarios)
                pd_igor_pgen = pd.read_csv(task.igor_fln_output_pgen, sep=';', index_col='seq_index')
                pd_airr_rearrangement['Pgen_estimate'] = pd_igor_pgen
                if fln_output is not None:
                    pd_airr_rearrangement.to_csv(fln_output, sep=';')
            except Exception as e:
                raise e

        # TODO: IN DEV RETURN OFFSET FOR V GENES ALIGNMENTS
        try:
            if b_V_offset:
                if os.path.isfile(task.igor_fln_align_V_alignments):
                    df_V_aligments = pd.read_csv(task.igor_fln_align_V_alignments, sep=';', usecols=['seq_index', 'gene_name', 'score', 'offset'])
                    V_nickname = task.mdl.parms.event_GeneChoice_V.nickname
                    gene_name_2_gene_id_dict = {v: k for k, v in task.mdl.parms[V_nickname]['name'].to_dict().items()}
                    gene_name_2_gene_id_dict # aver_dict
                    df_V_aligments[V_nickname] = df_V_aligments['gene_name'].apply(lambda x: gene_name_2_gene_id_dict[x])
                    del df_V_aligments['gene_name']
                    df_V_aligments.dropna(inplace=True)
                    df_V_offsets = df_V_aligments.groupby(by=['seq_index', V_nickname]).apply(lambda g: g.sort_values('score', ascending=False)) ['offset']
                    df_V_offsets = df_V_offsets.droplevel(None)


        except Exception as e:
            print("V alignments not "+task.igor_fln_align_V_alignments+"found.", e)
            pass

    except Exception as e:
        raise e
    else:
        if b_V_offset :
            return pd_airr_rearrangement, df_V_offsets
        else:
            return pd_airr_rearrangement
    finally:
        if batch_clean:
            task.run_clean_batch()
        tmp_dir.cleanup()


def evaluate_pgen(input_sequences:Union[str, pd.DataFrame, np.ndarray, Path],
             mdl:IgorModel, igor_wd=None, batch_clean=True, airr_format=False, pgen_columns:Union[None, list]=None,
                  igor_evaluate_dict_options:Union[None,dict]=None):
    """
    Evaluate input sequences with provided model
    :param input_sequences:Union[str, pd.DataFrame, np.ndarray, Path]
    :param mdl:IgorModel
    :param batch_clean: Remove all temporary files True by default.
    """
    # columns = ['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'pgen', 'scenario_rank', 'scenario_proba_cond_seq']
    if pgen_columns is None:
        if airr_format:
            pgen_columns = ['sequence_id', 'sequence', 'v_call', 'd_call', 'j_call', 'pgen']
        else:
            pgen_columns = ['Pgen_estimate']
    return evaluate(input_sequences, mdl, N_scenarios=1, igor_wd=igor_wd, batch_clean=batch_clean, airr_format=airr_format, igor_evaluate_dict_opts=igor_evaluate_dict_options)[pgen_columns]

def get_CDR3_VJ_anchors(mdl_ii:IgorModel, df_sequences, df_scenarios, df_V_offsets=None):
    try:
        if df_V_offsets is None:
            print("V offsets dataframe is required")
            return None
        # print(df_sequences)
        df_scenarios_VJ = df_scenarios.groupby(
            by=['seq_index', mdl_ii.event_GeneChoice_V_nickname, mdl_ii.event_GeneChoice_J_nickname], as_index=False).apply(
            lambda g: g[g['scenario_proba_cond_seq'] == g['scenario_proba_cond_seq'].max()])
        df_scenarios_VJ.index = df_scenarios_VJ.index.droplevel()
        # print(df_scenarios_VJ)

        df_scenarios_VJ['nt_sequence'] = df_sequences
        df_scenarios_VJ['V_offset'] = df_scenarios_VJ.apply(lambda x:
            df_V_offsets.loc[x.name].loc[x[mdl_ii.event_GeneChoice_V_nickname]],
            axis=1)
        df_scenarios_VJ['VJ_anchors_in_seq'] = df_scenarios_VJ.apply(
            lambda x: mdl_ii.get_CDR3_nt_pos_from_ps_scenario(x, x['V_offset']),
            axis=1)
        # print(df_scenarios_VJ)
        # df_scenario_complete['VJ_anchors_in_seq'] = df_scenario_complete.apply(lambda x: x['nt_sequence'][slice(*mdl_ii.get_CDR3_nt_pos_from_ps_scenario(x, x['V_offset']))], axis=1)
        df_scenarios_VJ['CDR3_nt'] = df_scenarios_VJ.apply(
            lambda x: x['nt_sequence'][slice(*x['VJ_anchors_in_seq'])],
            axis=1)
        df_scenarios_VJ['CDR3_aa'] = df_scenarios_VJ.apply(lambda x: dna_translate(x['CDR3_nt']), axis=1)

        df_scenarios_VJ['v_call'] = df_scenarios_VJ[mdl_ii.event_GeneChoice_V_nickname].apply(
            lambda x: mdl_ii.parms[mdl_ii.event_GeneChoice_V_nickname]['name'].loc[x])
        df_scenarios_VJ['j_call'] = df_scenarios_VJ[mdl_ii.event_GeneChoice_J_nickname].apply(
            lambda x: mdl_ii.parms[mdl_ii.event_GeneChoice_J_nickname]['name'].loc[x])

        return df_scenarios_VJ[ ['nt_sequence', 'V_offset', 'VJ_anchors_in_seq', 'CDR3_nt', 'CDR3_aa', 'v_call', 'j_call'] ]
    except Exception as e:
        raise e


def get_sonia_input(mdl_ii:IgorModel, df_sequences, df_scenarios, df_V_offsets=None):
    """
    Return dataframe to use it as input to sonia
    """
    try:
        df_CDR3_VJ = get_CDR3_VJ_anchors(mdl_ii, df_sequences, df_scenarios, df_V_offsets)
        df_input_sequences_sonia = df_CDR3_VJ[['CDR3_aa', 'v_call', 'j_call']]
        df_input_sequences_sonia['v_call'] = df_input_sequences_sonia['v_call'].apply(lambda x: x.split('*')[0])
        df_input_sequences_sonia['j_call'] = df_input_sequences_sonia['j_call'].apply(lambda x: x.split('*')[0])
        return list(df_input_sequences_sonia.to_numpy() )
    except Exception as e:
        raise e


#############################################
# Alias and Functions to get direct objects
def get_default_IgorModel(species, chain):
    """Return a default IGoR's model"""
    return IgorModel.load_default(species, chain)

def get_IgorModel_from_IgorRefGenome(ref_genome:IgorRefGenome):
    """Return a IgorModel from a IgorRefGenome"""
    return IgorModel.make_default_model_from_IgorRefGenome(ref_genome)

def get_imgt_list_species():
    """Return list of available species in IMGT website"""
    return IgorRefGenome.get_imgt_list_species()

def get_IgorRefGenome_VDJ_from_IMGT(imgt_species, imgt_chain):
    return IgorRefGenome.load_VDJ_from_IMGT_website(imgt_species, imgt_chain)

def get_IgorRefGenome_VJ_from_IMGT(imgt_species, imgt_chain):
    return IgorRefGenome.load_VJ_from_IMGT_website(imgt_species, imgt_chain)

RefGenome = IgorRefGenome
Model = IgorModel

################### SONIA METHODS ############ FIXME: CHANGE IT TO A BETTER LOCATION
try:
    import sonia
    from sonia.sequence_generation import SequenceGeneration
    from sonia.evaluate_model import EvaluateModel
    from sonia.sonia_leftpos_rightpos import SoniaLeftposRightpos

    def infer_SONIA_from_directory(mdl_dirname:Union[str, Path], data_seqs,
                            num_gen_seqs=int(2e5), reset_gen_seqs=True, custom_model_folder=None, add_error=False,
                            custom_error=None,
                            epochs=10, batch_size=5000, initialize=True, seed=None, validation_split=0.2, monitor=False,
                            verbose=0):
        """
        Return SONIA model
        :param mdl_dirname: Path, str or IgorModel object use to infer the selection model with SONIA
        :param data_seqs: List of arrays with format [ 'CDR3_aa', 'v_call', 'j_call' ]
        """
        qm_ii = None
        try:
            qm_ii = SoniaLeftposRightpos(data_seqs=data_seqs, custom_pgen_model=mdl_dirname)
            qm_ii.add_generated_seqs(num_gen_seqs=num_gen_seqs)
            qm_ii.infer_selection(epochs=epochs)
            # Save model if flag is active
            # if b_save_model:
            #     import shutil
            #     try:
            #         qm_ii.save_model(mdl_dirname)
            #     except shutil.SameFileError:
            #         pass
            #     except Exception as e:
            #         raise e
        except Exception as e:
            raise e

        return qm_ii

    def save_SONIA_mdl_dirname(qm:SoniaLeftposRightpos, mdl_dirname):
        import shutil
        try:
            qm.save_model(mdl_dirname)
        except shutil.SameFileError:
            pass
        except Exception as e:
            raise e


    aa_str_to_id_dict = {aa: aa_id for aa_id, aa in enumerate(list('ACDEFGHIKLMNPQRSTVWY'))}

    def ds_SONIA_model(model_dir_path):
        fln_features = model_dir_path + "/" + "features.tsv"
        df_features = pd.read_csv(fln_features, sep=',')
        b_vj = df_features['Feature'].apply(lambda x: x[0] == 'v')
        df_vj = df_features.loc[b_vj]

        df_vj.loc[:, 'V'] = df_vj['Feature'].apply(lambda x: x.split(';')[0])
        df_vj.loc[:, 'J'] = df_vj['Feature'].apply(lambda x: x.split(';')[1])
        ds_VJ = df_vj.set_index(['V', 'J']).to_xarray()

        b_len = df_features['Feature'].apply(lambda x: x[0] == 'l')
        df_len = df_features.loc[b_len]
        df_len.loc[:, 'Length'] = df_len['Feature'].apply(lambda x: int(x[1:]))
        df_len.set_index('Length', inplace=True)

        b_aa = df_features['Feature'].apply(lambda x: x[0] == 'a')
        df_aa = df_features.loc[b_aa]
        df_aa.loc[:, 'aa'] = df_aa['Feature'].apply(lambda x: str(x[1:2]))
        df_aa.loc[:, 'aa_id'] = df_aa['aa'].apply(lambda str_aa: aa_str_to_id_dict[str_aa])
        df_aa.loc[:, 'pos'] = df_aa['Feature'].apply(lambda x: int(x[2:]))

        ds_aa = df_aa.set_index(['aa_id', 'pos']).to_xarray()
        ds_aa_Left = ds_aa.loc[{'pos': ds_aa['pos'] >= 0}]
        ds_aa_Right = ds_aa.loc[{'pos': ds_aa['pos'] < 0}]

        return df_len, ds_aa_Left, ds_aa_Right, ds_VJ


except ImportError as error:
    pass
except Exception as exception:
    # Output unexpected Exceptions.
    print(exception, False)
    print(exception.__class__.__name__ + ": " + exception.message)



# if isinstance(mdl_dirname, IgorModel):
#     mdl = mdl_dirname
#     # Create a temporary directory
#     with tempfile.TemporaryDirectory(prefix='SONIA_inferring_', dir='.') as tmp_generate_dirname:
#         mdl.write_mdldata_dir(tmp_generate_dirname, sep=',', b_igor_directory=False)
#         qm_ii = SoniaLeftposRightpos(data_seqs=data_seqs, custom_pgen_model=tmp_generate_dirname)
#         qm_ii.add_generated_seqs(num_gen_seqs=num_gen_seqs)
#         qm_ii.infer_selection(epochs=epochs)

# def infer_SONIA(mdl, data_seqs,
#                 num_gen_seqs = int(2e5), reset_gen_seqs = True, custom_model_folder = None, add_error=False,custom_error=None,
#                 epochs = 10, batch_size=5000, initialize = True, seed = None,validation_split=0.2, monitor=False,verbose=0):
#     import sonia
#     from sonia.sequence_generation import SequenceGeneration
#     from sonia.evaluate_model import EvaluateModel
#     from sonia.sonia_leftpos_rightpos import SoniaLeftposRightpos
#
#     import tempfile
#     with tempfile.TemporaryDirectory(prefix='SONIA_inferring_', dir='.') as tmp_generate_dirname:
#         mdl.write_mdldata_dir(tmp_generate_dirname, sep=',', b_igor_directory=False)
#         qm_ii = SoniaLeftposRightpos(data_seqs=data_seqs, custom_pgen_model=tmp_generate_dirname)
#         qm_ii.add_generated_seqs(num_gen_seqs = num_gen_seqs)
#         qm_ii.infer_selection(epochs=epochs)
#     return qm_ii
