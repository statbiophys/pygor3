#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 11:02:28 2019

@author: alfaceor
"""

_flag_verbose_sqlite=False
import sqlite3
import csv
from Bio import SeqIO
#import IgorAlignment_data
from .IgorSQL import *
import numpy as np
import pandas as pd

###################
# TODO: Generalize to other kinds of database like postgres
###################

class IgorSqliteDB:
    """
    Class to create and load table or database with sequences
    """
    def __init__(self, igor_fln_db=None): #=None):
        # if flnIgorSQL is None:
        #     self.flnIgorSQL = "IgorDB.sql" # FIXME : VERY BAD WRAPP WITH INIT
        # else:
        #     self.flnIgorSQL = flnIgorSQL

        self.flnIgorSQL = "" # FIXME: this shouldn't be need it!

        self.fln_db  = None #flnIgorDB
        
        self.flnIgorIndexedSeq = ""
        self.flnIgorIndexedCDR3 = ""
        
        self.flnVGeneTemplate   = ""
        self.flnDGeneTemplate   = ""
        self.flnJGeneTemplate   = ""
        
        self.flnVAlignments     = ""
        self.flnDAlignments     = ""
        self.flnJAlignments     = ""
        
        self.flnVAnchors        = ""
        self.flnJAnchors        = ""

        self.flnIgorModel_Parms = ""
        self.flnIgorModel_Marginals = "" #

        self.flnIgorBestScenarios = ""
        self.flnIgorPgen = ""

        if igor_fln_db is not None:
            self.fln_db = igor_fln_db

        self.conn               = None

        self.sql_BestScenarios_cols_list = None
        self.sqlcmd_ins_bs = None
        self.sql_IgorBestScenarios_cols_name_type_list = None
        self.mdl = None
    
    #@classmethod
    def createSqliteDB(self, flnIgorDB):
        # FIXME: THIS METHOD SHOULD BE SOMETHING LIKE EXECUTE SQL SCRIPT and in particular
        """
        Create a SQLite database with the flnIgorDB sql script.
        """
        self.fln_db = flnIgorDB
        self.conn = None
        try:
            self.conn = sqlite3.connect(flnIgorDB)
            qry = open(self.flnIgorSQL, 'r').read()
            cur = self.conn.cursor()
            cur.executescript(qry)
            self.conn.commit()
            cur.close()
            #self.conn.close()
        except sqlite3.Error as e:
            print(e)

    @classmethod
    def create_db(cls, flnIgorDB):
        """
        Connect (or create if not exits) with filename
        """
        cls = IgorSqliteDB()
        cls.fln_db = flnIgorDB
        cls.connect_db()
        return cls

    def connect_db(self, flnIgorDB=None):
        """
        Connect (or create if not exits) to database
        """
        if flnIgorDB is not None:
            self.fln_db = flnIgorDB
        self.conn = sqlite3.connect(self.fln_db)

    def close_db(self):
        self.conn.close()

    def executescript(self, cur, str_query):
        cur.executescript(str_query)

    def execute_query(self, str_query):
        """
        Execute sql script in the SQLite database .
        """
        try:
            self.conn = sqlite3.connect(self.fln_db)
            cur = self.conn.cursor()
            cur.executescript(str_query)
            self.conn.commit()
            cur.close()
            # self.conn.close()
        except sqlite3.Error as e:
            print('str_query : ', str_query)
            print("sqlite3.ERROR : ", e)

    def execute_select_query(self, str_query):
        """
        Execute sql script in the SQLite database .
        """
        try:
            self.conn = sqlite3.connect(self.fln_db)
            cur = self.conn.cursor()
            cur.execute(str_query)
            # self.conn.commit()
            record = cur.fetchall()
            # print(record)
            # cur.close()
            return record
            # self.conn.close()
        except sqlite3.Error as e:
            # raise e
            e_message = "sqlite3.ERROR : " + str(self.fln_db)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def execute_select_query_fetchone(self, str_query):
        """
        Execute sql script in the SQLite database .
        """
        try:
            self.conn = sqlite3.connect(self.fln_db)
            cur = self.conn.cursor()
            cur.execute(str_query)
            #self.conn.commit()
            record = cur.fetchone()
            # print(record)
            # cur.close()
            return record
            #self.conn.close()
        except sqlite3.Error as e:
            print("sqlite3.ERROR : ", e)

    def createSqliteDB_tmp(self):
        # TODO: create database in base of IgorSQL scripts
        self.conn = None
        try:
            self.conn = sqlite3.connect(self.fln_db)
            qry = sqlcmd_ct['indexed_sequences']
            cur = self.conn.cursor()
            cur.executescript(qry)
            self.conn.commit()
            cur.close()
            # self.conn.close()
        except sqlite3.Error as e:
            print(e)

    # TODO: load database by receiving a list of identifiers defined in IgorDictionaries
    def load_db(self, **kwargs):
        """
        Return a parameter
        """
        for key in kwargs.keys():
            # Create a database.
            print( sqlcmd_ct[key] )
            #print(key)
            self.insert_IgorIndexedSeq_FromCSVline()

            # cur = self.conn.cursor()
            # try:
            #     cur.execute('BEGIN TRANSACTION')
            #     with open(flnIgorIndexedSeq) as fp:
            #         csvline = fp.readline()
            #         while csvline:
            #             csvline = fp.readline()
            #             # print(csvline)
            #             self.insert_IgorIndexedSeq_FromCSVline(cur, csvline)
            #     cur.execute('COMMIT')
            #     # self.conn.commit()
            # except sqlite3.Error as e:
            #     print(e)

    def load_VDJ_Database(self, flnIgorIndexedSeq, flnVGeneTemplate, flnDGeneTemplate, flnJGeneTemplate, flnVAlignments, flnDAlignments, flnJAlignments):
        
        self.load_IgorIndexedSeq_FromCSV(flnIgorIndexedSeq )
        
        self.load_IgorGeneTemplate_FromFASTA("V",flnVGeneTemplate)
        self.load_IgorGeneTemplate_FromFASTA("D",flnDGeneTemplate)
        self.load_IgorGeneTemplate_FromFASTA("J",flnJGeneTemplate)
        
        self.load_IgorAlignments_FromCSV("V", flnVAlignments)
        self.load_IgorAlignments_FromCSV("D", flnDAlignments)
        self.load_IgorAlignments_FromCSV("J", flnJAlignments)
    
    ###############################################
    ####  IgorIndexedSeq Table Methods
    ###############################################

    def load_IgorIndexedSeq_FromCSV(self, flnIgorIndexedSeq):
        """
        Insert indexed sequence in database from csv igor indexed_seqs file.
        :param conn:
        :param csvline:
        :return:
        """
        self.execute_query(sqlcmd_ct['indexed_sequences'])
        self.flnIgorIndexedSeq = flnIgorIndexedSeq
        cur = self.conn.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            with open(flnIgorIndexedSeq) as fp:
                csvline = fp.readline()
                while csvline:
                    csvline = fp.readline()
                    # print(csvline)
                    self.insert_IgorIndexedSeq_FromCSVline(cur, csvline)
            cur.execute('COMMIT')
            # self.conn.commit()
        except sqlite3.Error as e:
            print(e)

    def insert_IgorIndexedSeq_FromCSVline(self, cur, csvline):
        """
        Insert IGoR indexed_sequences on Database flnIgorDB
        :param csvline:
        """
        sql = ''' INSERT INTO IgorIndexedSeq(seq_index,sequence)
                  VALUES(?,?) '''
        
        csvline = csvline.replace('\n','')
        data = tuple(csvline.split(";"))
        if len(data) == 2:
            try:
                #cur = self.conn.cursor()
                cur.execute(sql, data)
                #self.conn.commit()
            except sqlite3.Error as e:
                print(data)
                print(e)
                pass

    def fetch_IgorIndexedSeq_records(self):
        """
        Fetch seq_index and sequence from Igor database.
        :param seq_index: string to specify the type of gene V, D or J
        :return:
        """

        #print(gene_name)
        sqlSelect = "SELECT seq_index, sequence FROM IgorIndexedSeq;"
        #print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        #if strGene == 'D':
        #    print(sqlSelect)
        records = cur.fetchall()
        return records

    def fetch_IgorIndexedSeq_By_seq_index(self, seq_index):
        """
        Fetch seq_index and sequence in Igor database.
        :param seq_index: string to specify the type of gene V, D or J
        :return: 
        """

        #print(gene_name)
        sqlSelect = "SELECT seq_index, sequence FROM IgorIndexedSeq WHERE seq_index = "+str(seq_index)+";"
        #print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        #if strGene == 'D':
        #    print(sqlSelect)
        record = cur.fetchone()
        return (record)

    def get_IgorIndexedSeq_By_seq_index(self, seq_index):
        record = self.fetch_IgorIndexedSeq_By_seq_index(seq_index)
        from .IgorIO import IgorIndexedSequence
        return IgorIndexedSequence.load_FromSQLRecord(record)

    def fetch_IgorIndexedSeq_By_seq_indexList(self, seq_indexList):
        """
        Fetch seq_index and sequence in Igor database.
        :param seq_index: string to specify the type of gene V, D or J
        :return: 
        """

        #print(gene_name)
        strSeq_indexList = str( tuple( sorted( set( seq_indexList ) ) )  )

        sqlSelect = "SELECT seq_index, sequence FROM IgorIndexedSeq WHERE seq_index IN "+strSeq_indexList+";"
        #print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        #if strGene == 'D':
        #    print(sqlSelect)
        record = cur.fetchall()
        return (record)

    def fetch_IgorIndexedSeq_indexes(self):
        seq_indexes_list = self.execute_select_query("SELECT seq_index FROM IgorIndexedSeq;")
        seq_indexes_list = list( map(lambda x: x[0], seq_indexes_list) )
        return seq_indexes_list

    def write_IgorIndexedSeq_to_CSV(self, flnIgorIndexedSeq, sep=";"):
        print("Saving indexed_sequences to file: ", flnIgorIndexedSeq)
        # print(gene_name)
        sqlSelect = "SELECT seq_index, sequence FROM IgorIndexedSeq;"
        # print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        # if strGene == 'D':
        #    print(sqlSelect)
        records = cur.fetchall()

        import pathlib
        import os
        directory_name = os.path.dirname(flnIgorIndexedSeq)
        pathlib.Path(directory_name).mkdir(parents=True, exist_ok=True)

        with open(flnIgorIndexedSeq, "w") as ofile:
            ofile.write("seq_index" + sep + "sequence" + "\n")
            for record in records:
                ofile.write(str(record[0])+sep+str(record[1])+"\n")

    def delete_IgorIndexedSeq_Tables(self):
        """
        Method to delete IgorIndexedSeq table.
        """
        tables_list = self.get_list_of_tables_with_name('IgorIndexedSeq')
        for tablename in tables_list:
            sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
            self.execute_query(sql_cmd_drop)

    ###############################################
    ####  IgorXGeneTemplate Tables Methods
    ###############################################

    def load_IgorGeneTemplate_FromFASTA(self, strGene, flnGeneTemplate):
        """
        Insert D Gene templates in database from fasta files used by IGoR.
        :param flnIgorGeneTemplate: Fasta file
        """
        # TODO: ADD A PANDAS OBJECT FOR RAPID ACCESS TO GENOMIC TEMPLATES.
        filename = {'V': self.flnVGeneTemplate, 'D' : self.flnDGeneTemplate, 'J': self.flnJGeneTemplate }
        filename[strGene] = flnGeneTemplate
        print(strGene, filename[strGene])
        # Create table if don't exits
        self.execute_query(sqlcmd_ct['genomic'+strGene + 's'])
        with open(filename[strGene], "r") as handle:
            for gene_id, bioRecord in enumerate(SeqIO.parse(handle, "fasta") ):
                self.insert_IgorGeneTemplate_FromBioRecord(strGene, gene_id, bioRecord)

    def insert_IgorGeneTemplate_FromBioRecord(self, strGene, gene_id, bioRecord):
        """
        Insert IGoR Gene template in Database flnIgorDB
        :param strGene: string to specify the type of gene V, D or J
        :param gene_id: id to identify the gene template
        :param bioRecord: Biopython record of the inserted sequence
        """
        sql = "INSERT INTO Igor"+strGene+"GeneTemplate("+strGene.lower()+"gene_id,gene_name,sequence) VALUES(?,?,?) "

        data = tuple([gene_id, str(bioRecord.description).strip(), str(bioRecord.seq)])
        try:
            cur = self.conn.cursor()
            if data is not None:
                cur.execute(sql, data)
            self.conn.commit()
        except sqlite3.Error as e:
            print(e)
            pass

    def fetch_IgorGeneTemplate_By_gene_name(self, strGene, gene_name):
        """
        Fetch Gene templates in database from fasta files used by IGoR.
        :param strGene: string to specify the type of gene V, D or J
        :param flnIgorGeneTemplate: Fasta file
        """
        sqlSelect = "SELECT * FROM Igor"+strGene.upper()+"GeneTemplate WHERE gene_name =\""+gene_name+"\";"
        #print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        #if strGene == 'D':
        #    print(sqlSelect)
        record = cur.fetchone()
        return (record)        
    
    def fetch_IgorGeneTemplate_By_gene_id(self, strGene, gene_id):
        """
        Fetch Gene templates in database from fasta files used by IGoR.
        :param strGene: string to specify the type of gene V, D or J
        :param gene_id: 
        """
        sqlSelect = "SELECT * FROM Igor"+strGene.upper()+"GeneTemplate WHERE "+strGene.lower()+"gene_id ="+str(gene_id)+";"
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        record = cur.fetchone()
        return (record)

    def write_IgorGeneTemplate_to_fasta(self, strGene, flnGeneTemplate, sep=";"):
        try:
            print("Saving gene templates to file: ", flnGeneTemplate)
            sqlSelect = "SELECT * FROM Igor"+strGene.upper()+"GeneTemplate;"
            cur = self.conn.cursor()
            cur.execute(sqlSelect)
            records = cur.fetchall()
            # str_file_header = "seq_index" + sep + "sequence" + "\n"
            import pathlib
            import os
            directory_name = os.path.dirname(flnGeneTemplate)
            pathlib.Path(directory_name).mkdir(parents=True, exist_ok=True)
            with open(flnGeneTemplate, "w") as ofile:
                # ofile.write(str_file_header)
                for record in records:
                    ofile.write(">" + str(record[1]) + "\n")
                    ofile.write(str(record[2]) + "\n")

        except Exception as e:
            raise e

    def delete_IgorGeneTemplate_Tables(self):
        """
        Method to delete IgorIndexedSeq table.
        """
        tables_list = self.get_list_of_tables_with_name('Igor%GeneTemplate')
        print("Tables to delete : ", tables_list)
        for tablename in tables_list:
            sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
            self.execute_query(sql_cmd_drop)


    def load_IgorGeneAnchors_FromCSV(self, strGene, flnGeneAnchors):
        # import pdb; pdb.set_trace()

        self.execute_query(sqlcmd_ct['gene'+strGene+'CDR3Anchors'])
        # FIXME: FILENAMES WERE NOT STORED
        filename = {'V': self.flnVAnchors, 'J': self.flnJAnchors}
        filename[strGene] = flnGeneAnchors

        print("Loading Gene Anchors from ", flnGeneAnchors)
        cur = self.conn.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            with open(filename[strGene], "r") as fp:
                csvline = fp.readline()
                while csvline:
                    csvline = fp.readline()
                    self.insert_IgorGeneAnchors_FromCSVline(strGene, cur, csvline)

            cur.execute('COMMIT')
            # self.conn.commit()
        except sqlite3.Error as e:
            print("sqlite3.ERROR : ", e)

    def insert_IgorGeneAnchors_FromCSVline(self, strGene, cur, csvline):
        """
        Insert IGoR indexed_CDR3_sequences in Database flnIgorDB
        :param csvline:
        """
        sql = "INSERT INTO Igor"+strGene+\
              "GeneCDR3Anchors("+strGene.lower()+"gene_id, anchor_index, function) VALUES(?,?,?) "

        csvline = csvline.replace('\n', '')
        data = csvline.split(";")

        try:
            data[0] = self.fetch_IgorGeneTemplate_By_gene_name(strGene, data[0])[0]
        except TypeError:
            pass
        except Exception as e:
            print("data : ", data)
            print("ERROR : ", e)
            pass

        try:
            if len(data) == 2:
                # NO FUNCTION SPECIFIED
                data = data + [None]
                cur.execute(sql, data)
            elif len(data) == 3:
                cur.execute(sql, data)
            else:
                pass #print(data)

        except sqlite3.Error as e:
            print(len(data), data)
            print("sqlite3.ERROR : ", e)
            pass

    # TODO: fetch records from database
    def fetch_IgorGeneAnchors_By_Gene(self, strGene):
        sqlSelect = "SELECT * FROM Igor" + strGene.upper() + "GeneCDR3Anchors"
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        records = cur.fetchall()
        return (records)

    def write_IgorGeneAnchors_to_CSV(self, strGene:str, flnGeneAnchors, sep=';'):
        """Export Gene anchors to csv IGoR file
        :param strGene: Gene letter
        """
        print("Saving IGoR's "+strGene+" gene anchors to file: ", flnGeneAnchors)
        sqlcmd_select_template = """
                    SELECT gene.gene_name,
                            anch.anchor_index, anch.function 
                    FROM Igor{upper}GeneTemplate gene 
                    LEFT JOIN Igor{upper}GeneCDR3Anchors anch 
                    ON gene.{lower}gene_id = anch.{lower}gene_id;
                    """

        sqlSelect = sqlcmd_select_template.format(upper=strGene.upper(), lower=strGene.lower())

        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        records = cur.fetchall()

        import pathlib
        import os
        directory_name = os.path.dirname(flnGeneAnchors)
        pathlib.Path(directory_name).mkdir(parents=True, exist_ok=True)

        str_file_header = "gene"+sep+"anchor_index"+sep+"function\n"
        with open(flnGeneAnchors, "w") as ofile:
            ofile.write(str_file_header)
            for record in records:
                if not record[1] is None:
                    ofile.write(str(record[0]) + sep)
                    ofile.write(str(record[1]) + sep)
                    ofile.write(str(record[2]) + '\n')

    def delete_IgorGeneAnchors_Tables(self):
        """
        Method to delete IgorIndexedSeq table.
        """
        tables_list = self.get_list_of_tables_with_name('Igor%GeneCDR3Anchors')
        print("Tables to delete : ", tables_list)
        for tablename in tables_list:
            sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
            self.execute_query(sql_cmd_drop)


    def fetch_IgorGenomicData_By_Gene(self, strGene):
        if strGene == 'D':
            sqlcmd_select_template = """
                                SELECT gene.{lower}gene_id, gene.gene_name, gene.sequence 
                                FROM Igor{upper}GeneTemplate gene;
                                """
        else:
            sqlcmd_select_template = """
                    SELECT gene.{lower}gene_id, gene.gene_name, gene.sequence, 
                            anch.anchor_index, anch.function 
                    FROM Igor{upper}GeneTemplate gene 
                    LEFT JOIN Igor{upper}GeneCDR3Anchors anch 
                    ON gene.{lower}gene_id = anch.{lower}gene_id;
                    """
        sqlSelect = sqlcmd_select_template.format(upper=strGene.upper(), lower=strGene.lower())

        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        records = cur.fetchall()
        return (records)

    def get_IgorGenomicDataFrame_dict(self):
        """return dataframes genomic_data """
        import pandas as pd
        genomic_data_dict = dict()
        V_records = self.fetch_IgorGenomicData_By_Gene("V")
        J_records = self.fetch_IgorGenomicData_By_Gene("J")

        columnas = ['id', 'gene_name', 'sequence', 'anchor_index', 'function']
        df_V = pd.DataFrame.from_records(V_records, columns=columnas, index='id')
        genomic_data_dict["V"] = df_V
        df_J = pd.DataFrame.from_records(J_records, columns=columnas, index='id')
        genomic_data_dict["J"] = df_J

        try:
            D_records = self.fetch_IgorGenomicData_By_Gene("D")
            columnas = ['id', 'gene_name', 'sequence']
            df_D = pd.DataFrame.from_records(D_records, columns=columnas, index='id')
            genomic_data_dict["D"] = df_D
        except Exception as e:
            print("No D genes were found in database ", self.fln_db)
            print("WARNING: ", e)
            pass


        return genomic_data_dict


    ###############################################
    ####  IgorXAlignments Tables Methods
    ###############################################
    
    def load_IgorAlignments_FromCSV(self, strGene, flnAlignments):
        """
        Insert Gene templates in database from fasta files used by IGoR.
        :param strGene: string to specify the type of gene V, D or J
        :param flnIgorGeneTemplate: Fasta file
        """

        # Load database from file
        strGene = strGene.upper()
        filename = {'V': self.flnVAlignments, 'D' : self.flnDAlignments, 'J': self.flnJAlignments}
        filename[strGene] = flnAlignments
        # Create table if don't exits
        self.execute_query(sqlcmd_ct[strGene + '_alignments'])
        try:
            cur = self.conn.cursor()
            cur.execute('BEGIN TRANSACTION')
            with open(filename[strGene]) as fp:
                csvline = fp.readline()
                while csvline:
                    csvline = fp.readline()
                    #print(csvline)
                    self.insert_IgorAlignments_FromCSVline(cur, strGene, csvline)
            
            cur.execute('COMMIT')
        except sqlite3.Error as e:
            print(e)

    def insert_IgorAlignments_FromCSVline(self, cur, strGene, csvline):
        """
        Insert IGoR Alignments on Database flnIgorDB
        :param strGene: string to specify the type of gene V, D or J
        :param csvline:
        """
        sqlBase = '''INSERT INTO Igor{}Alignments(
                seq_index, 
                {}gene_id, 
                score,
                offset,
                insertions,
                deletions,
                mismatches,
                length,
                offset_5_p_align,
                offset_3_p_align) 
                VALUES(?,?,?,?,?,?,?,?,?,?)
                '''
        
        strGene =  strGene.upper()
        sql = sqlBase.format(strGene.upper(), strGene.lower())
        # search on IgorGeneTemplate the corresponding Gene id 
        csvline = csvline.replace("{","[").replace("}","]").replace('\n','')
        csvlist = csvline.split(";")
        if len(csvlist) == 10 :
            gene_name = csvlist[1]
            gene_name = gene_name.strip()
            gene_id = self.fetch_IgorGeneTemplate_By_gene_name(strGene, gene_name)[0]
            #print(gene_id, gene_name)
            csvlist[1] = str(gene_id)
            data = tuple(csvlist)
            #print(data)
            try:
                #cur = self.conn.cursor()
                cur.execute(sql, data)
                #self.conn.commit()
            except sqlite3.Error as e:
                print(e)
                pass
        else:
            print(csvlist)

    def fetch_IgorAlignments_By_seq_index(self, strGene, seq_index, limit=None):
        """
        Fetch IgorAlignments from database by seq_index.
        :param strGene: string to specify the type of gene V, D or J
        :param seq_index: IgorIndexedSequences index
        """
        sqlSelect = "SELECT * FROM Igor"+strGene.upper()+"Alignments WHERE seq_index=="+str(seq_index)+" ORDER BY score DESC"
        if limit is not None:
            sqlSelect = sqlSelect + " LIMIT " + str(limit)
        sqlSelect = sqlSelect +";"
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        record = cur.fetchall()
        return (record)

    def fetch_IgorAlignments_By_seq_index_and_gene_name(self, strGene, seq_index, gene_name, limit=None):
        sqlSelect = "SELECT * FROM Igor" + strGene.upper() + "Alignments WHERE seq_index==" + str(
            seq_index) + " ORDER BY score DESC"
        if limit is not None:
            sqlSelect = sqlSelect + " LIMIT " + str(limit)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        record = cur.fetchall()
        return (record)

    def fetch_best_IgorAlignments_By_seq_index(self, strGene, seq_index):
        """
        Fetch IgorAlignments from database by seq_index.
        :param strGene: string to specify the type of gene V, D or J
        :param seq_index: IgorIndexedSequences index
        """
        sqlSelect = "SELECT * FROM Igor" + strGene.upper() + "Alignments WHERE seq_index==" + str(
            seq_index) + " ORDER BY score DESC LIMIT 1"
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        record = cur.fetchone()
        return (record)

    def get_best_IgorAlignment_data_By_seq_index(self, strGene, seq_index):
        from .IgorIO import IgorAlignment_data
        best_align_data_record = self.fetch_best_IgorAlignments_By_seq_index(strGene, seq_index)
        # print("best_align_data_record ", best_align_data_record)
        best_align_data = IgorAlignment_data.load_FromSQLRecord(best_align_data_record)
        best_align_data.strGene_class = strGene

        gene_record = self.fetch_IgorGeneTemplate_By_gene_id(best_align_data.strGene_class, best_align_data.gene_id)
        best_align_data.strGene_name = gene_record[1]
        best_align_data.strGene_seq = gene_record[2]
        return best_align_data

    # TODO: Create an IgorAlignment_data instance with a better sql query (join the necessary tables.)
    # TODO:

    def get_IgorAlignment_data_list_query(self, strGene, seq_index, where=None):
        sqlSelect = "SELECT * FROM Igor" + strGene.upper() + "Alignments WHERE seq_index==" + str(seq_index) + ""
        from .IgorIO import IgorAlignment_data
        align_data_records = self.fetch_IgorAlignments_By_seq_index(strGene, seq_index)
        align_data_list = list()
        for align_record in align_data_records:
            align_data = IgorAlignment_data.load_FromSQLRecord(align_record)
            align_data.strGene_class = strGene

            gene_record = self.fetch_IgorGeneTemplate_By_gene_id(align_data.strGene_class, align_data.gene_id)
            align_data.strGene_name = gene_record[1]
            align_data.strGene_seq = gene_record[2]
            align_data_list.append(align_data)
        return align_data_list

    def get_IgorAlignment_data_list_By_seq_index(self, strGene, seq_index, limit=None):
        from .IgorIO import IgorAlignment_data
        align_data_records = self.fetch_IgorAlignments_By_seq_index(strGene, seq_index, limit=limit)
        # best_align_data_record = self.fetch_best_IgorAlignments_By_seq_index(strGene, seq_index)
        # print("best_align_data_record ", align_data_records)
        align_data_list = list()
        for align_record in align_data_records:
            align_data = IgorAlignment_data.load_FromSQLRecord(align_record)
            align_data.strGene_class = strGene

            gene_record = self.fetch_IgorGeneTemplate_By_gene_id(align_data.strGene_class, align_data.gene_id)
            align_data.strGene_name = gene_record[1]
            align_data.strGene_seq = gene_record[2]
            align_data_list.append(align_data)
        return align_data_list

    def get_DataFrame_IgorAlignment_By_seq_index(self, strGene, seq_index, limit=None):
        records = self.fetch_IgorAlignments_By_seq_index(strGene, seq_index, limit=limit)
        pd.DataFrame.from_records(records)

    def appendList_IgorAlignments_data_By_seq_index(self, strGene_class, seq_index, alnDataList=None):
        """
        Append to a list of IgorAlignment_data objects given gene class ("V", "D", "J"), seq_index 
        append a list to append the objects.
        :param strGene_class: string to specify the type of gene V, D or J.
        :param seq_index: IgorIndexedSequences index.
        :param alnDataList: List of IgorAlignment_data objects.
        """
        if alnDataList == None:
            alnDataList = list()
        
        try:
            alignsSQLrecords = self.fetch_IgorAlignments_By_seq_index(strGene_class, seq_index)
            for alignSQLrecord in alignsSQLrecords:
                #print(alignSQLrecord)
                gene_id = alignSQLrecord[1]
                geneTemplateRecord = self.fetch_IgorGeneTemplate_By_gene_id(strGene_class, gene_id)
                strGene_name = geneTemplateRecord[1]
                strGene_seq  = geneTemplateRecord[2]
                aln_data = IgorAlignment_data.IgorAlignment_data.load_FromSQLRecord(alignSQLrecord, strGene_name=strGene_name)
                aln_data.strGene_class = strGene_class
                aln_data.strGene_seq   = strGene_seq
                alnDataList.append(aln_data)
                #print(aln_data.strGene_name, aln_data.score, aln_data.offset, aln_data.insertions)
            return alnDataList
        except Exception as e:
            print(e)
    
    def write_IgorAlignments2Fasta(self, alnDataList):
        print(alnDataList)

    def get_naive_sequence_from_IgorAligment_data(self, seq_index):
        # TODO: use self.dicts and IgorAlignment_data to reconstruct sequence
        # 1. get Indexed Sequence
        indexed_sequence = self.get_IgorIndexedSeq_By_seq_index(seq_index)
        indexed_sequence.offset = 0
        # 2. get Alignments

        # 3. Choose best V and J alignment
        best_v_align_data = self.get_best_IgorAlignment_data_By_seq_index('V', indexed_sequence.seq_index)
        best_j_align_data = self.get_best_IgorAlignment_data_By_seq_index('J', indexed_sequence.seq_index)

        # 4. if D exists then choose the best D gene where offsets are between CDR3 anchors
        d_align_data_list = self.get_IgorAlignment_data_list_By_seq_index('D', indexed_sequence.seq_index)

        # 5. Once D select the highest score.
        # 6. if there is an overlap in V or J segments then, removing that segment

        return ""

    def write_IgorAlignments_to_CSV(self, strGene, flnGeneTemplate, sep=";"):
        print("Exporting alignments to file: ", flnGeneTemplate)

        sqlcmd_select_template = """
            SELECT aln.seq_index, t_gene.gene_name, aln.score, aln.offset, aln.insertions, aln.deletions, 
            aln.mismatches, aln.length, aln.offset_5_p_align, aln.offset_3_p_align
            FROM Igor{upper}Alignments aln 
            INNER JOIN Igor{upper}GeneTemplate t_gene 
            ON aln.{lower}gene_id = t_gene.{lower}gene_id;
            """

        sqlSelect = sqlcmd_select_template.format(upper=strGene.upper(), lower=strGene.lower())
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        records = cur.fetchall()

        import pathlib
        import os
        directory_name = os.path.dirname(flnGeneTemplate)
        pathlib.Path(directory_name).mkdir(parents=True, exist_ok=True)

        str_file_header = "seq_index"+sep+"gene_name"+sep+"score"+sep+"offset"+sep+"insertions"+sep+"deletions"+sep+"mismatches"+sep+"length"+sep+"5_p_align_offset"+sep+"3_p_align_offset" + "\n"
        with open(flnGeneTemplate, "w") as ofile:
            ofile.write(str_file_header)
            for record in records:
                csvline = sep.join(map(str, record)).replace("[", "{").replace("]", "}")
                ofile.write(csvline + "\n")

    def delete_IgorAlignments_Tables(self):
        """
        Method to delete IgorAlignments table.
        """
        tables_list = self.get_list_of_tables_with_name('Igor%Alignments')
        for tablename in tables_list:
            sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
            self.execute_query(sql_cmd_drop)

    # def get_naive_alignment(self, seq_index):
    #     alnDataListV = db.fetch_IgorAlignments_By_seq_index('V', seq_index)
    #     v_align_data = p3.IgorAlignment_data.load_FromSQLRecord(alnDataListV[0])
    #     print(v_align_data.to_dict())
    #     v_gene_seq = db.fetch_IgorGeneTemplate_By_gene_id('V', v_align_data.gene_id)[2]
    #     print(v_gene_seq)

    def load_IgorIndexedCDR3_FromCSV(self, flnIgorIndexedCDR3):
        """
                Insert indexed CDR3 in database from csv igor indexed_seqs file.
                :param conn:
                :param csvline:
                :return:
                """
        self.execute_query(sqlcmd_ct['indexed_CDR3'])

        self.flnIgorIndexedCDR3 = flnIgorIndexedCDR3
        cur = self.conn.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            with open(self.flnIgorIndexedCDR3) as fp:
                csvline = fp.readline()
                while csvline:
                    csvline = fp.readline()
                    # print(csvline)
                    self.insert_IgorIndexedCDR3_FromCSVline(cur, csvline)
            cur.execute('COMMIT')
            # self.conn.commit()
        except sqlite3.Error as e:
            print(e)

    def insert_IgorIndexedCDR3_FromCSVline(self, cur, csvline):
        """
        Insert IGoR indexed_sequences on Database flnIgorDB
        :param csvline:
        """
        # seq_index;v_anchor;j_anchor;CDR3nt;CDR3aa
        sql = ''' INSERT INTO IgorIndexedCDR3(seq_index,v_anchor,j_anchor,CDR3,CDR3_aa)
                  VALUES(?,?,?,?,?) '''

        csvline = csvline.replace('\n', '')
        # FIXME: CHANGE THIS ONLY SAVE THE ANCHORS AND NOT THE CDR3 and CDR3_aa
        # data = tuple(csvline.split(";")[0:2]) # to pick just seq_index, v_anchor, j_anchor
        # if len(data) == 2:
        data = tuple(csvline.split(";"))
        # if len(data) == 2:
        try:
            # cur = self.conn.cursor()
            cur.execute(sql, data)
            # self.conn.commit()
        except sqlite3.Error as e:
            print(data)
            print(e)
            pass

    def fetch_IgorIndexedCDR3_By_seq_index(self, seq_index):
        # print(gene_name)
        sqlSelect = "SELECT * FROM IgorIndexedCDR3 WHERE seq_index = " + str(seq_index) + ";"
        # print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        # if strGene == 'D':
        #    print(sqlSelect)
        record = cur.fetchone()
        return (record)

    def write_IgorIndexedCDR3_to_CSV(self, flnIgorIndexedCDR3, sep=';'):
        print("Saving IGoR's indexed CDR3 to file : ", flnIgorIndexedCDR3)
        sqlSelect = "SELECT * FROM IgorIndexedCDR3;"
        # print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        # if strGene == 'D':
        #    print(sqlSelect)
        records = cur.fetchall()

        import pathlib
        import os
        directory_name = os.path.dirname(flnIgorIndexedCDR3)
        pathlib.Path(directory_name).mkdir(parents=True, exist_ok=True)

        str_file_header = "seq_index;v_anchor;j_anchor;CDR3;CDR3_aa" + "\n"
        with open(flnIgorIndexedCDR3, "w") as ofile:
            ofile.write(str_file_header)
            for record in records:
                csvline = sep.join(map(str, record)).replace("[", "{").replace("]", "}")
                ofile.write(csvline + "\n")

    def delete_IgorIndexedCDR3_Tables(self):
        """
        Method to delete IgorAlignments table.
        """
        tables_list = self.get_list_of_tables_with_name('IgorIndexedCDR3')
        for tablename in tables_list:
            sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
            self.execute_query(sql_cmd_drop)

    # def get_IgorIndexedCDR3_By_seq_index(self, seq_index):
    #     record = self.fetch_IgorIndexedCDR3_By_seq_index(seq_index)
    #     return IgorIndexedSequence.load_FromSQLRecord(record)

    ############################
    def load_IgorModel_Parms(self, mdl_parms): #fln_model_parms):
        ################## TODO: Insert Event_list
        self.execute_query(sqlcmd_ct['MP_Event_list'])
        # self.flnIgorModel_Parms = fln_model_parms
        cur = self.conn.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')

            from .IgorIO import IgorModel_Parms
            # mdl_parms = IgorModel_Parms(model_parms_file=fln_model_parms)
            for event in mdl_parms.Event_list:
                event_dict = event.to_dict()
                # print("-"*50, event.nickname)
                event_realizations_list = event_dict.pop('realizations')
                event_dict['realizations_table'] = "IgorER_" + event_dict['nickname']
                self.insert_IgorRec_Event_FromDict(cur, event_dict)
                # Create database with realizations
                script_table = sqlcmd_ct['ER_event_template'].format(event_dict['realizations_table'])
                # insert realizations from realizations object
                try:
                    # Create table of realizations
                    cur.execute(script_table)
                    # Insert realizations
                    for realization in event_realizations_list:
                        # print(realization.to_dict())
                        self.insert_IgorEvent_realization_FromDict(cur, event_dict['realizations_table'], realization.to_dict())
                except sqlite3.Error as e:
                    print("Problem with realization", realization.to_dict())
                    print(e)

            cur.execute('COMMIT')
            # self.conn.commit()
        except sqlite3.Error as e:
            print("Couldn't create IgorModel_Parms table")
            print(e)

        ################## TODO: Insert Edges XXX
        # Create database if not exist.
        self.execute_query(sqlcmd_ct['MP_Edges'])
        cur = self.conn.cursor()
        sql = 'INSERT INTO IgorMP_Edges(parent_event,child_event) VALUES (?,?)'
        cur.execute('BEGIN TRANSACTION')
        for edge in mdl_parms.Edges:
            # Insert in Edges table
            try:
                # cur = self.conn.cursor()
                data = (mdl_parms.dictNameNickname[ edge[0] ], mdl_parms.dictNameNickname[ edge[1] ])
                print(data)
                cur.execute(sql, data)
                # self.conn.commit()
            except sqlite3.Error as e:
                print(edge, data)
                print(e)
                exit()
                # pass
        cur.execute('COMMIT')

        ################## TODO: Insert ErrorRate
        try:
            self.execute_query(sqlcmd_ct['MP_ErrorRate'])
            cur = self.conn.cursor()
            columns = ', '.join(mdl_parms.ErrorRate_dict.keys())
            placeholders = ':' + ', :'.join(mdl_parms.ErrorRate_dict.keys())
            sql = 'INSERT INTO IgorMP_ErrorRate (%s) VALUES (%s)' % (columns, placeholders)
            cur.execute('BEGIN TRANSACTION')
            # print(mdl_parms.ErrorRate_dict)
            cur.execute(sql, mdl_parms.ErrorRate_dict)
            cur.execute('COMMIT')
        except Exception as e:
            print("Couldn't insert ErrorRate in database")
            print(e)

    def insert_IgorRec_Event_FromDict(self, cur, event_dict: dict):
        # print(event_dict.keys())
        columns = ', '.join(event_dict.keys())
        placeholders = ':' + ', :'.join(event_dict.keys())
        sql = 'INSERT INTO IgorMP_Event_list (%s) VALUES (%s)' % (columns, placeholders)
        # print(sql)
        try:
            # cur = self.conn.cursor()
            cur.execute(sql, event_dict)
            # self.conn.commit()
        except sqlite3.Error as e:
            print(event_dict)
            print(e)
            pass

    def insert_IgorEvent_realization_FromDict(self, cur, event_table, event_realization_dict):
        # TODO: IgorEvent_realization
        columns = ', '.join(event_realization_dict.keys())
        placeholders = ':' + ', :'.join(event_realization_dict.keys())
        tmp_sql = 'INSERT INTO '+event_table+' (%s) VALUES (%s)'
        sql = tmp_sql % (columns, placeholders)
        # print(sql)
        try:
            # cur = self.conn.cursor()
            cur.execute(sql, event_realization_dict)
            # self.conn.commit()
        except sqlite3.Error as e:
            print(event_realization_dict)
            print(e)
            pass

    def load_IgorModel_Marginals(self, mdl_xdata:dict):
        for event_nickname in mdl_xdata.keys():
            da = mdl_xdata[event_nickname]
            lista = list(dict(da.coords).keys())
            # print(lista)
            for dim in da.dims:
                lista.remove(dim)
            da = da.drop(lista)
            # Create table of not real marginals
            if da.attrs['event_type'] == 'DinucMarkov':
                self.execute_query( sqlcmd_ct_Model_Marginals_DinucMarkov(event_nickname, ['x', 'y']) )
                cur = self.conn.cursor()
                try:
                    cur.execute('BEGIN TRANSACTION')
                    df = da.to_dataframe(name='P')
                    df.reset_index(inplace=True)
                    # Rename names of dataframe
                    dicto2rename = dict()
                    for col in df.columns:
                        if not col == 'P':
                            dicto2rename[col] = "id_" + col
                    df.rename(dicto2rename, axis='columns', inplace=True)

                    mm_table = "IgorMM_" + event_nickname
                    columns_list = list(df.columns)
                    str_columns = ",".join(columns_list)
                    placeholders = ':' + ', :'.join(columns_list)
                    sql = 'INSERT INTO ' + mm_table + ' (%s) VALUES (%s);' % (str_columns, placeholders)
                    #print("sql : ", sql)
                    for index, row in df.iterrows():
                        cur.execute(sql, row.to_dict())

                    cur.execute('COMMIT')
                except sqlite3.Error as e:
                    print("-" * 40)
                    print("ERROR: Can't insert in " + event_nickname + " table")
                    print(e)
                    print("-"*40)

            else:
                tmp_list = [event_nickname] + da.attrs['parents']
                da = da.transpose(*tmp_list)
                # Create table for model marginals
                self.execute_query( sqlcmd_ct_Model_Marginals(event_nickname, da.attrs['parents']) )

                cur = self.conn.cursor()
                try:
                    cur.execute('BEGIN TRANSACTION')
                    df = da.to_dataframe(name='P')
                    df.reset_index(inplace=True)

                    # Rename names of dataframe
                    dicto2rename = dict()
                    for col in df.columns:
                        if not col == 'P':
                            dicto2rename[col] = "id_"+col
                    df.rename(dicto2rename, axis='columns', inplace=True)

                    mm_table = "IgorMM_"+event_nickname
                    columns_list = list(df.columns)
                    str_columns = ",".join(columns_list )
                    placeholders = ':' + ', :'.join(columns_list)
                    sql = 'INSERT INTO ' + mm_table + ' (%s) VALUES (%s)' % (str_columns, placeholders)

                    for index, row in df.iterrows():
                        cur.execute(sql, row.to_dict())
                    cur.execute('COMMIT')
                    # self.conn.commit()
                except sqlite3.Error as e:
                    print("-" * 40)
                    print("ERROR: Can't insert in "+event_nickname+" table")
                    print(e)
                    print("-" * 40)

    def insert_IgorModel_Marginals_FromDict(self, cur, event_table, event_realization_dict):
        # TODO: IgorEvent_realization
        columns = ', '.join(event_realization_dict.keys())
        placeholders = ':' + ', :'.join(event_realization_dict.keys())
        tmp_sql = 'INSERT INTO '+event_table+' (%s) VALUES (%s)'
        sql = tmp_sql % (columns, placeholders)
        # print(sql)
        try:
            # cur = self.conn.cursor()
            cur.execute(sql, event_realization_dict)
            # self.conn.commit()
        except sqlite3.Error as e:
            print(event_realization_dict)
            print(e)
            pass

    def load_IgorModel_FromTXT(self, flnIgorModel_Parms, flnIgorModel_Marginals):
        self.flnIgorModel_Parms = flnIgorModel_Parms
        self.flnIgorModel_Marginals = flnIgorModel_Marginals
        from .IgorIO import IgorModel
        mdl = IgorModel(model_parms_file=flnIgorModel_Parms, model_marginals_file=flnIgorModel_Marginals)
        self.load_IgorModel(self, mdl)

    def load_IgorModel(self, mdl):
        print("Loading parms to database")
        self.load_IgorModel_Parms(mdl.parms)
        print("Loading marginals to database")
        self.load_IgorModel_Marginals(mdl.xdata)

    def write_IgorModel_to_TXT(self, flnIgorModel_Parms, flnIgorModel_Marginals):
        print("Saving IGoR's model parms to file: ", flnIgorModel_Parms, " and model marginals to file: ", flnIgorModel_Marginals)
        mdl_parms = self.get_IgorModel_Parms()
        mdl_parms.write_model_parms(flnIgorModel_Parms)
        mdl_marginals = self.get_IgorModel_Marginals()
        mdl_marginals.write_model_marginals(flnIgorModel_Marginals, mdl_parms)

    def write_IgorModel_Parms_to_TXT(self, flnIgorModel_Parms):
        # from .IgorIO import IgorModel_Parms
        print("Saving IGoR's model parms to file: ", flnIgorModel_Parms) #
        mdl_parms = self.get_IgorModel_Parms()
        mdl_parms.write_model_parms(flnIgorModel_Parms)

    def write_IgorModel_Marginals_to_TXT(self, flnIgorModel_Marginals):
        from .IgorIO import IgorModel_Parms, IgorModel_Marginals
        print("Saving IGoR's and model marginals to file: ", flnIgorModel_Marginals)
        mdl_parms = self.get_IgorModel_Parms()
        mdl_marginals = self.get_IgorModel_Marginals()
        mdl_marginals.write_model_marginals(flnIgorModel_Marginals, mdl_parms)


    def delete_IgorModel_Tables(self):
        for tabla in sql_tablename_patterns_dict['model']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            for tablename in tables_list:
                sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
                counts = self.execute_query(sql_cmd_drop)

    def delete_IgorModel_Parms_Tables(self):
        for tabla in ['IgorMP_%', 'IgorER_%']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            for tablename in tables_list:
                sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
                counts = self.execute_query(sql_cmd_drop)

    def delete_IgorModel_Marginals_Tables(self):
        for tabla in ['IgorMM_%']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            for tablename in tables_list:
                sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
                counts = self.execute_query(sql_cmd_drop)


    def load_IgorBestScenarios_FromCSV(self, flnIgorBestScenarios, mdl=None):
        if mdl is None:
            mdl = self.get_IgorModel()

        self.flnIgorBestScenarios = flnIgorBestScenarios
        with open(self.flnIgorBestScenarios, "r") as fp:
            #### CREATE DATABASE
            ## Get best scenarios header file
            igor_bsfileHeader = fp.readline().replace('\n', '').replace('Mismatches', 'mismatches')
            file_header_list = igor_bsfileHeader.split(';')
            bs_events_name_list = file_header_list[3:-1]
            ## Change event_name to event_nickname using dictionary.
            # Get dictionary to change names to nicknames
            events_name_nickname_dict = mdl.parms.get_event_dict('name', 'nickname')
            print("events_name_nickname_dict", events_name_nickname_dict)
            bs_events_nickname_list = [ events_name_nickname_dict[event_name] for event_name in bs_events_name_list]

            # Now generate sql_id, sql_type and foreign_key in the file order
            events_dict_list = list()
            sql_list_ct_table_cols = list()
            sql_list_ct_table_foreign_keys = list()
            sql_list_ins_table = list()
            self.sql_list_of_listtypes = list()
            for ii, event_nickname in enumerate(bs_events_nickname_list) :
                event = mdl.parms.get_Event(event_nickname)
                event_dict = event.to_dict()
                del event_dict['realizations']
                if event.event_type == 'DinucMarkov':
                    event_dict['sql_id'] = event.nickname
                    event_dict['sql_type'] = "text"
                    sql_list_ct_table_cols.append(event_dict['sql_id']+" "+event_dict['sql_type'])
                    sql_list_ins_table.append(event_dict['sql_id'])
                    self.sql_list_of_listtypes.append(3+ii)
                else:
                    event_dict['sql_id'] = 'id_'+event.nickname
                    event_dict['sql_type'] = "integer"
                    event_dict['sql_foreign_key'] = "FOREIGN KEY (id_{}) REFERENCES IgorER_{} (id)".format(event_nickname, event_nickname)
                    sql_list_ct_table_cols.append(event_dict['sql_id'] + " " + event_dict['sql_type'])
                    sql_list_ins_table.append(event_dict['sql_id'])
                    sql_list_ct_table_foreign_keys.append(event_dict['sql_foreign_key'])
                events_dict_list.append(event_dict)

            self.sql_list_of_listtypes.append(3+len(bs_events_nickname_list)) # mismatches position

            # Now make script to create table in database and create table
            self.execute_query(
                sqlcmd_ct['BestScenario_template'].format(','.join(sql_list_ct_table_cols),
                                                          ','.join(sql_list_ct_table_foreign_keys))
            )

            ### INSERT VALUES
            sql_cols_list =['seq_index', 'scenario_rank', 'scenario_proba_cond_seq'] + sql_list_ins_table + ['mismatches', 'mismatcheslen']
            sql_str_cols = ','.join(sql_cols_list)
            sql_str_aux = ','.join(['?' for col in sql_cols_list])

            self.sql_BestScenarios_cols_list = sql_cols_list
            self.sqlcmd_ins_bs = ''' INSERT INTO IgorBestScenarios({}) VALUES({}) '''.format(sql_str_cols, sql_str_aux)

            #print(self.sqlcmd_ins_bs)

            # INSERT VALUES IN DATABASE
            cur = self.conn.cursor()
            try:
                cur.execute('BEGIN TRANSACTION')
                with open(self.flnIgorBestScenarios) as fp:
                    csvline = fp.readline()
                    while csvline:
                        csvline = fp.readline()
                        self.insert_IgorBestScenarios_FromCSVline(cur, csvline)
                cur.execute('COMMIT')
                # self.conn.commit()
            except sqlite3.Error as e:
                print(e)

        self.gen_IgorBestScenarios_cols_list()

    def insert_IgorBestScenarios_FromCSVline(self, cur, csvline):
        # sql = ''' INSERT INTO IgorBestScenarios({}) VALUES({}) '''
        csvline = csvline.replace('\n', '')
        csvlist = csvline.split(";")
        # FIXME: FIND A BETTER WAY TO REPLACE THE TYPE FOR INSERTION
        for ii in range(len(csvlist)):
            if ii in self.sql_list_of_listtypes:
                csvlist[ii] = csvlist[ii].replace("(", "[").replace(")", "]")  # lists vd_dinucl, dj_dinucl, mismatches
            else:
                csvlist[ii] = csvlist[ii].replace("(", "").replace(")", "")

        # print(csvlist)

        # data[-1] convert it to list and
        try:
            data = csvlist + [len(eval(csvlist[-1]))]
            try:
                # cur = self.conn.cursor()
                cur.execute(self.sqlcmd_ins_bs, data)
                # self.conn.commit()
            except sqlite3.Error as e:
                print(data)
                print(e)
                pass
        except Exception as e:
            if _flag_verbose_sqlite:
                print("ERROR: insert_IgorBestScenarios_FromCSVline")
                print(csvlist)
                print(e)
            pass

    def write_IgorBestScenarios_to_CSV(self, flnIgorBestScenarios, mdl=None, sep=';'):
        """
        seq_index;scenario_rank;scenario_proba_cond_seq;GeneChoice_V_gene_Undefined_side_prio7_size168;GeneChoice_J_gene_Undefined_side_prio7_size16;GeneChoice_D_gene_Undefined_side_prio6_size3;Deletion_V_gene_Three_prime_prio5_size21;Deletion_D_gene_Five_prime_prio5_size21;Deletion_D_gene_Three_prime_prio5_size21;Deletion_J_gene_Five_prime_prio5_size21;Insertion_VD_genes_Undefined_side_prio4_size41;DinucMarkov_VD_genes_Undefined_side_prio3_size16;Insertion_DJ_gene_Undefined_side_prio2_size41;DinucMarkov_DJ_gene_Undefined_side_prio1_size16;Mismatches
        99;1;0.0401246;(109);(12);(2);(9);(16);(4);(7);(5);(2,2,0,3,0);(5);(0,2,0,3,3);()
        99;2;0.0401246;(109);(12);(0);(9);(9);(7);(7);(5);(2,2,0,3,0);(5);(0,2,0,3,3);()
        99;3;0.0100351;(109);(12);(2);(10);(13);(7);(7);(0);();(11);(0,2,0,3,3,2,2,2,0,0,3);()
        99;4;0.0100351;(109);(12);(2);(9);(14);(7);(7);(0);();(11);(0,2,0,3,3,2,2,2,0,0,3);()
        """
        print("Saving IGoR's best scenarios to file: ", flnIgorBestScenarios)
        if mdl is None:
            mdl = self.get_IgorModel()

        ### MAKE THE HEADER with name using nicknames
        # Transform the nicknames with id
        db_events_nickname_name_dict = dict()
        for event in mdl.parms.Event_list:
            if event.event_type == 'DinucMarkov':
                db_events_nickname_name_dict[event.nickname] = event.name
            else:
                db_events_nickname_name_dict["id_"+event.nickname] = event.name
        db_events_nickname_name_dict['mismatches'] = 'Mismatches'
        print(db_events_nickname_name_dict)

        list_bs_columns_db = [el[0] for el in self.get_columns_type_of_tables('IgorBestScenarios')[:-1]]
        tmp_list = list()
        for aaa in list_bs_columns_db:
            if aaa in db_events_nickname_name_dict.keys():
                tmp_list.append(db_events_nickname_name_dict[aaa])
            else:
                tmp_list.append(aaa)
        # str_file_header
        # str_file_header = sep.join(list_bs_columns_db)
        str_file_header = sep.join(tmp_list)+"\n"

        str_sql_columns_to_query = ",".join(list_bs_columns_db)
        sqlSelect = "SELECT " + str_sql_columns_to_query + " FROM IgorBestScenarios;"
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        records = cur.fetchall()

        import pathlib
        import os
        directory_name = os.path.dirname(flnIgorBestScenarios)
        pathlib.Path(directory_name).mkdir(parents=True, exist_ok=True)

        # str_file_header = "seq_index;gene_name;score;offset;insertions;deletions;mismatches;length;5_p_align_offset;3_p_align_offset" + "\n"
        with open(flnIgorBestScenarios, "w") as ofile:
            ofile.write(str_file_header)
            for record in records:
                csvline = sep.join(map(str, record)).replace("[", "{").replace("]", "}")
                ofile.write(csvline + "\n")

    def delete_IgorBestScenarios_Tables(self):
        tablenames_list = self.get_list_of_tables_with_name('IgorBestScenarios')
        for tablename in tablenames_list:
            sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
            self.execute_query(sql_cmd_drop)

    def export_IgorBestScenarios_to_AIRR(self, flnAIRR_arrangement, mdl=None, sep='\t'):
        """Export Igor best scenarios to airr rearrangement format"""
        if mdl is None:
            mdl = self.get_IgorModel()

        from .AIRR import AIRR_VDJ_rearrangement
        import airr

        import pathlib
        import os
        directory_name = os.path.dirname(flnAIRR_arrangement)
        pathlib.Path(directory_name).mkdir(parents=True, exist_ok=True)

        # IF VDJ THEN:
        b_D_gene = ( len([ event.event_type  for event in mdl.parms.Event_list if event.nickname == 'd_gene']) > 0 )

        if b_D_gene:
            airr_fields = AIRR_VDJ_rearrangement.list_of_fields()
            airr_rearrangement_writer = airr.create_rearrangement(flnAIRR_arrangement, fields=airr_fields)
            for seq_index, sequence in self.fetch_IgorIndexedSeq_records():
                scenarios_list = self.get_IgorBestScenarios_By_seq_index_IgorModel(seq_index, mdl)
                for scenario in scenarios_list:
                    v_best_aln = self.get_best_IgorAlignment_data_By_seq_index("V", seq_index)

                    airr_rearrangement_dict = mdl.get_AIRR_VDJ_rearragement_dict_from_scenario(scenario, sequence, v_offset=v_best_aln.offset)
                    airr_rearrangement_dict['scenario_rank'] = scenario.scenario_rank
                    airr_rearrangement_dict['scenario_proba_cond_seq'] = scenario.scenario_proba_cond_seq
                    airr_rearrangement_dict['pgen'] = self.fetch_IgorPgen_By_seq_index(seq_index)[1]
                    cdr3_record = self.fetch_IgorIndexedCDR3_By_seq_index(seq_index)
                    airr_rearrangement_dict['junction'] = cdr3_record[3]
                    airr_rearrangement_dict['junction_aa'] = cdr3_record[4]
                    airr_rearrangement_writer.write(airr_rearrangement_dict)

            airr_rearrangement_writer.close()

        else:
            # FIXME: IF VJ THEN:
            airr_fields = AIRR_VDJ_rearrangement.list_of_fields()
            airr_rearrangement_writer = airr.create_rearrangement(flnAIRR_arrangement, fields=airr_fields)
            for seq_index, sequence in self.fetch_IgorIndexedSeq_records():
                scenarios_list = self.get_IgorBestScenarios_By_seq_index_IgorModel(seq_index, mdl)
                for scenario in scenarios_list:
                    v_best_aln = self.get_best_IgorAlignment_data_By_seq_index("V", seq_index)

                    airr_rearrangement_dict = mdl.get_AIRR_VJ_rearragement_dict_from_scenario(scenario, sequence,
                                                                                               v_offset=v_best_aln.offset)
                    airr_rearrangement_dict['scenario_rank'] = scenario.scenario_rank
                    airr_rearrangement_dict['scenario_proba_cond_seq'] = scenario.scenario_proba_cond_seq
                    airr_rearrangement_dict['pgen'] = self.fetch_IgorPgen_By_seq_index(seq_index)[1]
                    cdr3_record = self.fetch_IgorIndexedCDR3_By_seq_index(seq_index)
                    airr_rearrangement_dict['junction'] = cdr3_record[3]
                    airr_rearrangement_dict['junction_aa'] = cdr3_record[4]
                    airr_rearrangement_writer.write(airr_rearrangement_dict)

            airr_rearrangement_writer.close()

        # n_d_5_del = self.mdlParms.Event_dict[strEv].loc[self.id_d_5_del]['value']
        # name_D = self.mdlParms.Event_dict[strEv].loc[self.id_d_gene]['name']
        # v_call d_call c_call v_score d_score j_score v_cigar j_evalue v_identity v_start v_end v_germ_start v_germ_end np1_seq np1_length np2_seq np2_length

        # from pygor3.IgorIO import IgorScenario
        # # OPEN FILE
        # for seq_index in self.fetch_IgorIndexedSeq_indexes():
        #     scenarios = self.get_IgorBestScenarios_By_seq_index(seq_index)
        #     sequence = self.fetch_IgorIndexedSeq_By_seq_index(seq_index)
        #     for scenario in scenarios:
        #         scenario.export_to_AIRR_line()
        #
        #
        # ### MAKE THE HEADER with name using nicknames
        # # Transform the nicknames with id
        # db_events_nickname_name_dict = dict()
        # for event in mdl.parms.Event_list:
        #     if event.event_type == 'DinucMarkov':
        #         db_events_nickname_name_dict[event.nickname] = event.name
        #     else:
        #         db_events_nickname_name_dict["id_" + event.nickname] = event.name
        # db_events_nickname_name_dict['mismatches'] = 'Mismatches'
        # print(db_events_nickname_name_dict)
        #
        # list_bs_columns_db = [el[0] for el in self.get_columns_type_of_tables('IgorBestScenarios')[:-1]]
        # tmp_list = list()
        # for aaa in list_bs_columns_db:
        #     if aaa in db_events_nickname_name_dict.keys():
        #         tmp_list.append(db_events_nickname_name_dict[aaa])
        #     else:
        #         tmp_list.append(aaa)
        # # str_file_header
        # # str_file_header = sep.join(list_bs_columns_db)
        # str_file_header = sep.join(tmp_list) + "\n"
        #
        # str_sql_columns_to_query = ",".join(list_bs_columns_db)
        # sqlSelect = "SELECT " + str_sql_columns_to_query + " FROM IgorBestScenarios;"
        # cur = self.conn.cursor()
        # cur.execute(sqlSelect)
        # records = cur.fetchall()
        #
        # # str_file_header = "seq_index;gene_name;score;offset;insertions;deletions;mismatches;length;5_p_align_offset;3_p_align_offset" + "\n"
        # with open(flnIgorBestScenarios, "w") as ofile:
        #     ofile.write(str_file_header)
        #     for record in records:
        #         csvline = sep.join(map(str, record)).replace("[", "{").replace("]", "}")
        #         ofile.write(csvline + "\n")



    ###### PGEN
    def load_IgorPgen_FromCSV(self, flnIgorPgen):
        print(flnIgorPgen)
        self.execute_query(sqlcmd_ct['Pgen'])
        self.flnIgorPgen = flnIgorPgen
        cur = self.conn.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            with open(self.flnIgorPgen) as fp:
                csvline = fp.readline()
                while csvline:
                    csvline = fp.readline()
                    # print(csvline)
                    self.insert_load_IgorPgen_FromCSVline(cur, csvline)
            cur.execute('COMMIT')
            # self.conn.commit()
        except sqlite3.Error as e:
            print(e)

    def insert_load_IgorPgen_FromCSVline(self, cur, csvline):
        sql = ''' INSERT INTO IgorPgen(seq_index,Pgen_estimate)
                          VALUES(?,?) '''
        csvline = csvline.replace('\n', '')
        data = tuple(csvline.split(";"))
        try:
            cur.execute(sql, data)
        except sqlite3.Error as e:
            print(data)
            print(e)
            pass

    def fetch_IgorPgen(self):
        sqlSelect = "SELECT * FROM IgorPgen;"
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        records = cur.fetchall()
        return records

    def fetch_IgorPgen_By_seq_index(self, seq_index):
        sqlSelect = "SELECT * FROM IgorPgen WHERE seq_index = " + str(seq_index) + ";"
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        record = cur.fetchone()
        return record

    def write_IgorPgen_to_CSV(self, flnIgorPgen, sep=';'):
        print("Writing Pgen to file: ", flnIgorPgen)
        import pathlib
        import os
        directory_name = os.path.dirname(flnIgorPgen)
        pathlib.Path(directory_name).mkdir(parents=True, exist_ok=True)
        str_header = "seq_index"+sep+"Pgen_estimate\n"
        with open(flnIgorPgen, "w") as ofile:
            # seq_index,Pgen_estimate
            ofile.write(str_header)
            for record in self.fetch_IgorPgen():
                try:
                    ofile.write(str(record[0])+sep+str(record[1])+"\n" )
                except Exception as e:
                    print("ERROR: write_IgorPgen_to_CSV ", e)

    def delete_IgorPgen_Tables(self):
        tablenames_list = self.get_list_of_tables_with_name('IgorPgen')
        for tablename in tablenames_list:
            sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
            self.execute_query(sql_cmd_drop)


    ###### return IGoR Model
    def get_IgorModel(self):
        from .IgorIO import IgorModel
        try:
            mdl_parms = self.get_IgorModel_Parms()
            mdl_marginals = self.get_IgorModel_Marginals()
            mdl = IgorModel.load_from_parms_marginals_object( mdl_parms, mdl_marginals )
            return mdl
        except Exception as e:
            # raise e
            e_message = "ERROR: IgorSqliteDB.get_IgorModel " + str(self.igor_fln_db)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def get_IgorModel_Marginals(self):
        print("-"*5, "Marginals", "-"*5)
        from .IgorIO import IgorModel_Marginals
        marginals = IgorModel_Marginals()

        # So there are two different ways to implement this.
        # 1. use a model parms to get all the nicknames and explore the name tables
        # 2. Check query in db the nicknames
        # 3. How to get all the tables with the prefix "IgorMM_"
        # And the best way to do this is to ask the database (2)
        # Get all the list of nicknames:
        nickname_list = self.execute_select_query("SELECT nickname FROM IgorMP_Event_list;")
        nickname_list = [nickname[0] for nickname in nickname_list]

        marginals.network_dict = dict()
        marginals.marginals_dict = dict()

        for event_nickname in nickname_list:
            col_names = self.get_table_colsname_list("IgorMM_"+event_nickname)
            records = self.execute_select_query("SELECT * FROM IgorMM_"+event_nickname+";")

            print(event_nickname) # , records)
            events_tmp = [x[3:] for x in col_names if x != 'P'] # remove the 3 characters 'id_'
            marginals.network_dict[event_nickname] = events_tmp
            df = pd.DataFrame(np.array(records), columns=events_tmp+['P'])
            dimensions = list()
            for event_name in events_tmp:
                dim_ev = len(df[event_name].unique())
                dimensions.append(dim_ev)

            np_marginals = df['P'].to_numpy().reshape(dimensions)

            # FIXME TRANSPOSE WITH DIMENSIONS
            axes_for_transpose = list()
            aux_axes_for_transpose = list()

            new_dict_list = list()
            aux_new_dict_list = list()
            for ii, str_event in enumerate(events_tmp):
                if not str_event == event_nickname:
                    axes_for_transpose.append(ii)
                    new_dict_list.append(str_event)
                else:
                    aux_axes_for_transpose.append(ii)
                    aux_new_dict_list.append(str_event)

            marginals.network_dict[event_nickname] = new_dict_list + aux_new_dict_list
            marginals.marginals_dict[event_nickname] = np.transpose(np_marginals, tuple(axes_for_transpose+aux_axes_for_transpose))

            # FIXME:
            # 1. GEt position of event_nickname in list marginals.network_dict[event_nickname]
            # and put it at the end. marginals.network_dict[event_nickname]

        return marginals

    def get_table_colsname_list(self, tablename):
        sqlcmd_qry = sqlcmd_template_qry_cols.format(tablename=tablename)
        col_names = self.execute_select_query(sqlcmd_qry)
        return [colname[0] for colname in col_names]

    def get_IgorModel_Parms(self):
        try:
            from .IgorIO import IgorModel_Parms
            mdl_parms = IgorModel_Parms()
            # tb_IgorMP_Event_list_columns = self.execute_select_query(
            #     "SELECT name FROM pragma_table_info('IgorMP_Event_list');")
            # tb_IgorMP_Event_list_columns = [aa[0] for aa in tb_IgorMP_Event_list_columns]
            # event_type, seq_type, seq_side, priority, nickname

            mdl_parms.Event_list = self.get_Event_list()
            dict_nickname_name = mdl_parms.get_event_dict('nickname', 'name')
            mdl_parms.Edges = list()
            for edge in self.get_Edges():
                mdl_parms.Edges.append( [dict_nickname_name[edge[0]], dict_nickname_name[edge[1]] ])

            mdl_parms.ErrorRate_dict = self.get_ErrorRate_dict()

            mdl_parms.gen_EventDict_DataFrame()
            return mdl_parms
        except Exception as e:
            e_message = "ERROR: IgorSqliteDB.get_IgorModel_Parms " + str(self.igor_fln_db)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def get_Event_list(self):
        try:
            from .IgorIO import IgorRec_Event
            from .IgorIO import IgorEvent_realization
            sql_cmd = "SELECT event_type, seq_type, seq_side, priority, nickname, realizations_table FROM IgorMP_Event_list;"
            Event_list = list()
            for rec in self.execute_select_query(sql_cmd):
                event = IgorRec_Event(*rec[:-1])
                str_realization_table = rec[-1]
                realization_dbrecords = self.execute_select_query("SELECT * FROM " + str_realization_table + ";")
                for realization_dbrec in realization_dbrecords:
                    realization = IgorEvent_realization()
                    if event.event_type == "GeneChoice":
                        realization.id = int(realization_dbrec[0])
                        realization.value = realization_dbrec[1]
                        realization.name = realization_dbrec[2]
                    elif event.event_type == "DinucMarkov":
                        realization.id = int(realization_dbrec[0])
                        realization.value = realization_dbrec[1]
                    else:
                        realization.id = int(realization_dbrec[0])
                        realization.value = int(realization_dbrec[1])
                    event.add_realization(realization)
                #print(event)
                Event_list.append(event)

            return Event_list
        except Exception as e:
            e_message = "ERROR: IgorSqliteDB.get_Event_list " + str(self.igor_fln_db)
            import sys
            raise type(e)(str(e) + '\n' + e_message).with_traceback(sys.exc_info()[2])

    def get_Edges(self):
        sql_cmd = "SELECT parent_event, child_event FROM IgorMP_Edges;"
        Edges = self.execute_select_query(sql_cmd)
        return Edges

    def get_ErrorRate_dict(self):
        ErrorRate_dict = dict()
        sql_cmd = "SELECT error_type, error_values FROM IgorMP_ErrorRate;"
        ErrorRate_record = self.execute_select_query(sql_cmd)[0]
        if _flag_verbose_sqlite:
            print("ErrorRate_record : ", ErrorRate_record)
        ErrorRate_dict['error_type'] = ErrorRate_record[0]
        ErrorRate_dict['error_values'] = ErrorRate_record[1]
        return ErrorRate_dict

    def get_list_of_tables_with_name(self, table_name_pattern):
        cmd_qry = """
        SELECT name FROM sqlite_master WHERE type ='table' AND name NOT LIKE 'sqlite_%'
	    AND name LIKE '{tablename}';
        """
        sql_cmd_qry = cmd_qry.format(tablename=table_name_pattern)
        records = self.execute_select_query(sql_cmd_qry)
        tables_list = [record[0] for record in records]
        return tables_list

    def get_dict_of_Igortablename_sql(self):
        cmd_qry = """
        SELECT name, sql FROM sqlite_master WHERE type ='table' AND name NOT LIKE 'sqlite_%'
        AND name LIKE 'Igor%';
        """
        records = self.execute_select_query(cmd_qry)
        tablename_sql_dict = dict()
        for record in records:
            tablename_sql_dict[record[0]] = record[1]

        return tablename_sql_dict

    def delete_IgorModel_Tables(self):
        tablepatterns_list = ["IgorER_%", "IgorMM_%", "IgorMP_%"]
        for tablepattern in tablepatterns_list:
            tablenames_list = self.get_list_of_tables_with_name(tablepattern)
            for tablename in tablenames_list:
                sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
                self.execute_query(sql_cmd_drop)

    #########
    def get_columns_type_of_tables(self, tablename):
        select_qry = "SELECT name, type FROM pragma_table_info('{tablename}');".format(tablename=tablename)
        return self.execute_select_query(select_qry)

    def gen_IgorBestScenarios_cols_list(self):
        self.sql_IgorBestScenarios_cols_name_type_list = self.execute_select_query("SELECT name, type FROM pragma_table_info('IgorBestScenarios');")

    def fetch_IgorBestScenarios_By_events_dict(self, events_id_list_dict):
        print(events_id_list_dict)
        sqlselect = "SELECT * FROM IgorBestScenarios"
        sql_tmp_ANDs_list = list()
        print(sqlselect)
        for key, values in events_id_list_dict.items():
            print(key, values)
            sql_tmp_ANDs_list.append( key+ " IN "+ "(" + ",".join(map(str,values))+ ")" )

        if len(sql_tmp_ANDs_list) == 0:
            sql_where = ";"
        else:
            sql_where = " WHERE " +  " AND ".join(sql_tmp_ANDs_list) + ";"

        return self.execute_select_query(sqlselect + sql_where)

    def fetch_IgorBestScenarios_By_seq_index(self, seq_index):
        sqlSelect = "SELECT * FROM IgorBestScenarios WHERE seq_index==" + str(
            seq_index) #+ " ORDER BY score DESC LIMIT 1"
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        records = cur.fetchall()
        return (records)

    def get_IgorBestScenariosDataframe_By_seq_index(self, seq_index):
        bs_records = self.fetch_IgorBestScenarios_By_seq_index(seq_index)
        import pandas as pd
        self.gen_IgorBestScenarios_cols_list()
        columnas = [ item[0] for item in self.sql_IgorBestScenarios_cols_name_type_list]
        return pd.DataFrame.from_records(bs_records, columns=columnas)


    def get_IgorBestScenarios_By_seq_index(self, seq_index):
        from .IgorIO import IgorScenario
        # FIXME: THIS METHOD SHOULD BE EVALUATED ONCE AFTER DATA IS LOADED WITH A TRY IF TABLE EXIST
        #  TO AVOID THE UNNECESSARY RERUNNING.
        self.gen_IgorBestScenarios_cols_list()
        scenarios_list = list()
        for record in self.fetch_IgorBestScenarios_By_seq_index(seq_index):
            scenario = IgorScenario.load_FromSQLRecord(record, self.sql_IgorBestScenarios_cols_name_type_list)
            scenarios_list.append(scenario)
        return scenarios_list

    # FIXME : FIND A BETTER DESING FOR THIS FUNCTION
    def get_IgorBestScenarios_By_seq_index_IgorModel(self, seq_index, mdl):
        scenarios_mdl_list = list()
        scenarios_list = self.get_IgorBestScenarios_By_seq_index(seq_index)
        for scen in scenarios_list:
            scen_dict = scen.realizations_ids_dict
            for event in mdl.parms.Event_list:
                if not (event.event_type == 'DinucMarkov'):
                    scen_dict[event.nickname] = scen_dict.pop('id_' + event.nickname)
            scen.realizations_ids_dict = scen_dict
            scenarios_mdl_list.append(scen)

        return scenarios_mdl_list


    def fetch_AIRR_arrangement_By_seq_index(self, seq_index, mdl=None):
        sqlSelect = "SELECT * FROM IgorBestScenarios WHERE seq_index==" + str(
            seq_index) #+ " ORDER BY score DESC LIMIT 1"

        from .IgorIO import IgorScenario

        # TODO: NOW AIRR SCENARIOS
        # 1. Load model used for scenarios.
        if mdl is None:
            mdl = self.get_IgorModel()

        change_dicto = dict()
        for event in mdl.parms.Event_list:
            if event.event_type == 'DinucMarkov':
                change_dicto[event.nickname] = event.nickname
            else:
                change_dicto['id_'+event.nickname] = event.nickname

        # 2. Get the list of columns
        scenarios_col_list = self.get_columns_type_of_tables('IgorBestScenarios')
        scenarios_db_list = list()
        # So scenarios in database have ids, so for use the model I just want the realizations dont care about the rest
        evento = mdl.parms.get_Event(event.nickname)
        from .IgorIO import IgorRec_Event
        evento = IgorRec_Event()
        evento.realizations[realization_id]
        # for the scenarios I only have
        from collections import namedtuple



        for col, tipo in scenarios_col_list:
            print(col, tipo)
            #check in colname in mdl
            if col in change_dicto:
                change_dicto[col]

        scenarios_list = self.get_IgorBestScenarios_By_seq_index_IgorModel()
        # For each scenario in list pass it to airr
        # genomics_dict = self.get_IgorGenomicDataFrame_dict()
        # # 3. Make a map using the column list with model
        #
        # sqlSelect_VDJ= """
        #     SELECT Tscenarios.seq_index as sequence_id,
        #         Tscenarios.scenario_rank,
        #         Tscenarios.scenario_proba_cond_seq,
        #         Tindexed.sequence,
        #         T_v_choice.name as v_call,
        #         T_d_gene.name as d_call,
        #         T_j_choice.name as j_call
        #     FROM IgorBestScenarios Tscenarios
        #         JOIN IgorIndexedSeq Tindexed,
        #             IgorER_v_choice T_v_choice,
        #             IgorER_d_gene T_d_gene,
        #             IgorER_j_choice T_j_choice
        #     WHERE Tindexed.seq_index==Tscenarios.seq_index
        #         AND Tindexed.seq_index == 0
        #         AND Tscenarios.id_v_choice == T_v_choice.id
        #         AND Tscenarios.id_d_gene == T_d_gene.id
        #         AND Tscenarios.id_j_choice == T_j_choice.id;
        #     """
        #
        #
        # cur = self.conn.cursor()
        # cur.execute(sqlSelect)
        # records = cur.fetchall()
        # return (records)

    def calc_IgorBestScenarios_average_of(self, scenario_function, indices_list=None):
        from .IgorIO import IgorScenario
        self.gen_IgorBestScenarios_cols_list()
        if indices_list is None:
            indices_list = self.fetch_IgorIndexedSeq_indexes()
        print("len: ", len(indices_list))
        function_average = 0
        for indx in indices_list:
            # for indx in indexes_list:
            # aln_data = db.get_IgorAlignment_data_list_By_seq_index('V', indx)
            bs_data_list = self.fetch_IgorBestScenarios_By_seq_index(indx)
            tmp_average = 0
            scenario_norm_factor = 0
            for bs_data in bs_data_list:
                bs = IgorScenario.load_FromSQLRecord(bs_data, self.sql_IgorBestScenarios_cols_name_type_list)
                tmp_average = tmp_average + bs.scenario_proba_cond_seq * scenario_function(bs)
                scenario_norm_factor = scenario_norm_factor + bs.scenario_proba_cond_seq
            # FIXME: PERFORMANCE WILL INCREASE IF WE MAKE THE NORMALIZATION BEFORE?
            tmp_average = tmp_average / scenario_norm_factor
            # print("scenario_norm_factor: ", scenario_norm_factor)
            function_average = function_average + tmp_average

        return (function_average / (len(indices_list)))

    # FIXME: What should be the correct method to calculate the average?
    def FIXME_calc_IgorBestScenarios_average_of(self, scenario_function):
        from .IgorIO import IgorScenario
        self.gen_IgorBestScenarios_cols_list()
        indexes_list = self.fetch_IgorIndexedSeq_indexes()
        print("len: ", len(indexes_list))
        function_average = 0
        seq_pgen_tuple_list = self.fetch_IgorPgen()

        pgen_normalization = 0
        for indx, pgen in seq_pgen_tuple_list:
            # for indx in indexes_list:
            # aln_data = db.get_IgorAlignment_data_list_By_seq_index('V', indx)
            bs_data_list = self.fetch_IgorBestScenarios_By_seq_index(indx)
            tmp_average = 0
            scenario_norm_factor = 0
            pgen_normalization = pgen_normalization + pgen
            for bs_data in bs_data_list:
                bs = IgorScenario.load_FromSQLRecord(bs_data, self.sql_IgorBestScenarios_cols_name_type_list)
                tmp_average = tmp_average + bs.scenario_proba_cond_seq * scenario_function(bs)
                scenario_norm_factor = scenario_norm_factor + bs.scenario_proba_cond_seq
            tmp_average = tmp_average / scenario_norm_factor
            # function_average = function_average + pgen * tmp_average
            function_average = function_average + tmp_average

        return (function_average / (pgen_normalization * len(indexes_list)))

    def calc_IgorBestScenarios_sum_of(self, scenario_function):
        from .IgorIO import IgorScenario
        self.gen_IgorBestScenarios_cols_list()
        indexes_list = self.fetch_IgorIndexedSeq_indexes()
        print("len: ", len(indexes_list))
        function_sum = 0
        for indx in indexes_list:
            bs_data_list = self.fetch_IgorBestScenarios_By_seq_index(indx)
            tmp_average = 0
            for bs_data in bs_data_list:
                bs = IgorScenario.load_FromSQLRecord(bs_data, self.sql_IgorBestScenarios_cols_name_type_list)
                tmp_average = tmp_average + scenario_function(bs)
            function_sum = function_sum + tmp_average
            # function_average = function_average + tmp_average

        return (function_sum )


    ################## DATABASE UTILITIES ##################
    def Q_sequences_in_db(self):
        flag_db = False
        for tabla in sql_tablename_patterns_dict['read_seqs']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            # print(tables_list)
            if len(tables_list) > 0: flag_db = True

        return flag_db

    def Q_ref_genome_in_db(self):
        flag_db = False
        for tabla in sql_tablename_patterns_dict['ref_genome']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            # print(tables_list)
            if len(tables_list) > 0: flag_db = True

        return flag_db

    def Q_CDR3_Anchors_in_db(self, strGene):
        flag_table = False
        tables_list = self.get_list_of_tables_with_name('Igor' + strGene + 'GeneCDR3Anchors')
        if len(tables_list) == 1: flag_table = True
        return flag_table

    def Q_ref_genome_in_db_by_gene(self, strGene):
        flag_table = False
        tables_list = self.get_list_of_tables_with_name('Igor'+strGene+'GeneTemplate')
        if len(tables_list) == 1 : flag_table = True
        return flag_table

    def Q_align_in_db(self):
        flag_db = False
        for tabla in sql_tablename_patterns_dict['align']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            # print(tables_list)
            if len(tables_list) > 0: flag_db = True

        return flag_db

    def Q_align_in_db_by_gene(self, strGene):
        flag_table = False
        tables_list = self.get_list_of_tables_with_name('Igor' + strGene + 'Alignments')
        if len(tables_list) == 1: flag_table = True
        return flag_table


    def Q_model_in_db(self):
        flag_db = False
        for tabla in sql_tablename_patterns_dict['model']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            # print(tables_list)
            if len(tables_list) > 0: flag_db = True

        return flag_db

    def Q_output_in_db(self):
        flag_db = False
        for tabla in sql_tablename_patterns_dict['output']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            # print(tables_list)
            if len(tables_list) > 0: flag_db = True

        return flag_db

    def Q_IgorPgen_in_db(self):
        flag_db = False
        for tabla in ['IgorPgen']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            # print(tables_list)
            if len(tables_list) > 0: flag_db = True
        return flag_db

    def Q_IgorBestScenarios_in_db(self):
        flag_db = False
        for tabla in ['IgorBestScenarios']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            # print(tables_list)
            if len(tables_list) > 0: flag_db = True
        return flag_db

    def list_from_db(self):
        print("=== Sequences tables igor-reads: ")
        for tabla in sql_tablename_patterns_dict['read_seqs']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            for tablename in tables_list:
                counts = self.execute_select_query(
                    "SELECT COUNT(*) FROM {tablename};".format(tablename=tablename))
                print(tablename, " : ", counts[0][0])

        print("=== Genomes References tables igor-genomes: ")
        for tabla in sql_tablename_patterns_dict['ref_genome']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            for tablename in tables_list:
                counts = self.execute_select_query(
                    "SELECT COUNT(*) FROM {tablename};".format(tablename=tablename))
                print(tablename, " : ", counts[0][0])

        print("=== Alignments tables igor-alignments: ")
        for tabla in sql_tablename_patterns_dict['align']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            for tablename in tables_list:
                counts = self.execute_select_query(
                    "SELECT COUNT(*) FROM {tablename};".format(tablename=tablename))
                print(tablename, " : ", counts[0][0])

        print("=== Model tables igor-model: ")
        for tabla in sql_tablename_patterns_dict['model']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            for tablename in tables_list:
                counts = self.execute_select_query(
                    "SELECT COUNT(*) FROM {tablename};".format(tablename=tablename))
                print(tablename, " : ", counts[0][0])

        print("=== Output tables igor-pgen and igor-scenarios: ")
        for tabla in sql_tablename_patterns_dict['output']:
            tables_list = self.get_list_of_tables_with_name(tabla)
            for tablename in tables_list:
                counts = self.execute_select_query(
                    "SELECT COUNT(*) FROM {tablename};".format(tablename=tablename))
                print(tablename, " : ", counts[0][0])


    def copytable_from_source(self, tablename_to_copy, fln_source_db):
        #fln_source_db = igortask.igor_db.flnIgorDB
        sql_ct_from_source_db = """
                        ATTACH '{fln_source_db}' AS sourcedb;
                        INSERT INTO {tablename_to_copy} SELECT * FROM sourcedb.{tablename_to_copy};
                        DETACH sourcedb
                        """
        sql_cmd = sql_ct_from_source_db.format(fln_source_db=fln_source_db, tablename_to_copy=tablename_to_copy)
        self.execute_query(sql_cmd)

    # FIXME: FINISH IT I'M HERE
    def attach_table_from_db(self, flnIgorPgen_to_attach):
        # 1. Attach database flnIgorPgen_to_attach to current connection
        # 2. Create table in flnIgorPgen_to_attach
        # 3. Make a select in base table and insert to new table in flnIgorPgen_to_attach
        # 4. Dettach database flnIgorPgen_to_attach
        """
        ATTACH 'ahorasi.db' as otherdb;
        DETACH otherdb;
        """
        tablename_sql_dict = dict()
        sqlcmd_ct_from_otherdb_template = """
                    SELECT name, sql FROM otherdb.sqlite_master
                    WHERE type ='table' AND name NOT LIKE 'sqlite_%'
                    AND name = '{tablename}';
                    """
        sqlcmd_ct_from_otherdb_template = """
                    SELECT name, sql FROM sqlite_master
                    WHERE type ='table' AND name NOT LIKE 'sqlite_%'
                    AND name = '{tablename}';
                    """


        category = 'sequences'
        if category == 'sequences':
            sql_tablename_patterns_dict['read_seqs']
        if category == 'ref_genome':
            sql_tablename_patterns_dict['ref_genome']
        if category == 'align':
            sql_tablename_patterns_dict['align']
        if category == 'model':
            sql_tablename_patterns_dict['model']

        print( sql_tablename_patterns_dict.keys())
        self.get_dict_of_Igortablename_sql()
        sqlcmd = sqlcmd_ct_from_otherdb_template.format(tablename='IgorPgen')
        record = self.execute_select_query_fetchone(sqlcmd)
        print(record[0])
        print(record[1])


    # def rm_sequences_from_db(self):
    #     for tabla in sql_tablename_patterns_dict['read_seqs']:
    #         tables_list = self.get_list_of_tables_with_name(tabla)
    #         for tablename in tables_list:
    #             sql_cmd_drop = "DROP TABLE IF EXISTS {tablename};".format(tablename=tablename)
    #             self.execute_query(sql_cmd_drop)









