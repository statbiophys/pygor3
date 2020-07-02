#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 16:52:51 2019

@author: alfaceor
"""
#import IgorModel
import pandas as pd
import numpy as np
import sqlite3
from .IgorIO import *

"""
seq_index
scenario_rank
scenario_proba_cond_seq
GeneChoice_V_gene_Undefined_side_prio7_size35
GeneChoice_J_gene_Undefined_side_prio7_size14
GeneChoice_D_gene_Undefined_side_prio6_size2
Deletion_V_gene_Three_prime_prio5_size21
Deletion_D_gene_Five_prime_prio5_size21
Deletion_D_gene_Three_prime_prio5_size21
Deletion_J_gene_Five_prime_prio5_size23
Insertion_VD_genes_Undefined_side_prio4_size31
DinucMarkov_VD_genes_Undefined_side_prio3_size16
Insertion_DJ_gene_Undefined_side_prio2_size31
DinucMarkov_DJ_gene_Undefined_side_prio1_size16
Mismatches
"""




class IgorSqliteDBBestScenariosVDJ:
    """
    Class to create and load table or database with sequences
    """
    def __init__(self, flnIgorSQLBestScenarios):
        #IgorInferenceDB_VDJ.sql
        # FIXME: LOAD DATA FROM PACKAGE ROOT
        self.flnIgorSQLBestScenarios = flnIgorSQLBestScenarios #"IgorDBBestScenariosVDJ.sql"
        self.flnIgorDB  = "" #flnIgorDB
        
        self.conn               = None
    
    def createSqliteDB(self, flnIgorDB):
        """
        Create a SQLite database with the flnIgorDB sql script.
        """
        self.flnIgorDB = flnIgorDB
        self.conn = None

        try:
            self.conn = sqlite3.connect(flnIgorDB)
            qry = open(self.flnIgorSQLBestScenarios, 'r').read()
            cur = self.conn.cursor()
            cur.executescript(qry)
            self.conn.commit()
            cur.close()
            #self.conn.close()
        except sqlite3.Error as e:
            print(e)
    
    
    ###############################################
    ####  IgorBestScenarios Table Methods
    ###############################################
    
    def load_IgorBestScenariosVDJ_FromCSV(self, flnIgorBestScenarios):
        """
        Insert bestScenarios database in database from csv file.
        :param flnIgorBestScenarios:
        :return:
        """
        self.flnIgorBestScenarios = flnIgorBestScenarios        
        
        cur = self.conn.cursor()
        try:
            cur.execute('BEGIN TRANSACTION')
            with open(self.flnIgorBestScenarios) as fp:
                csvline = fp.readline()
                while csvline:
                    csvline = fp.readline()
                    #print(csvline)
                    self.insert_IgorBestScenariosVDJ_FromCSVline(cur, csvline)
            cur.execute('COMMIT')
            #self.conn.commit()
        except sqlite3.Error as e:
            print(e)

    
    def insert_IgorBestScenariosVDJ_FromCSVline(self, cur, csvline):
        """
        Insert IGoR indexed_sequences on Database flnIgorDB
        :param csvline:
        """
        sql = ''' INSERT INTO IgorDBBestScenariosVDJ (
                        seq_index,
                        scenario_rank,
                        scenario_proba_cond_seq,
                        id_v_choice,
                        id_j_choice,
                        id_d_gene,
                        id_v_3_del,
                        id_d_5_del,
                        id_d_3_del,
                        id_j_5_del,
                        id_vd_ins,
                        vd_dinucl,
                        id_dj_ins,
                        dj_dinucl,
                        mismatches,
                        mismatcheslen
                    ) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) '''
        
        csvline = csvline.replace('\n','')
        csvlist = csvline.split(";")
        for ii, csvelem in enumerate(csvlist):
            if ii in [11, 13, 14]:
                csvlist[ii] = csvlist[ii].replace("(", "[").replace(")", "]")  # lists vd_dinucl, dj_dinucl, mismatches
            else:
                csvlist[ii] = csvlist[ii].replace("(", "").replace(")", "")
        
        if len(csvlist) == 15:
            csvlist.append(str(len(eval(csvlist[14])))) # mismatcheslen
            data = tuple(csvlist)
            try:
                #cur = self.conn.cursor()
                cur.execute(sql, data)
                #self.conn.commit()
            except sqlite3.Error as e:
                print(csvlist)
                print(e)
                pass
        

    def fetch_IgorBestScenariosVDJ_By_seq_index(self, seq_index):
        """
        Fetch seq_index and sequence in Igor database.
        :param seq_index: string to specify the type of gene V, D or J
        :return: 
        """

        #print(gene_name)
        sqlSelect = "SELECT * FROM IgorDBBestScenariosVDJ WHERE seq_index = "+str(seq_index)+";"
        #print(sqlSelect)
        cur = self.conn.cursor()
        cur.execute(sqlSelect)
        #if strGene == 'D':
        #    print(sqlSelect)
        record = cur.fetchall()
        return (record)
    
    
    
    
        

#    
#    def __init__(self, flnBestScenarios, model_parms_file=None, model_marginals_file=None):
#        
#        self.seq_index = XXXX
#        self.scenario_rank = XXXX
#        self.scenario_proba_cond_seq = XXXX
#        self.GeneChoice_V_gene_Undefined_side_prio7_size35 = XXXX
#        self.GeneChoice_J_gene_Undefined_side_prio7_size14 = XXXX
#        self.GeneChoice_D_gene_Undefined_side_prio6_size2 = XXXX
#        self.Deletion_V_gene_Three_prime_prio5_size21 = XXXX
#        self.Deletion_D_gene_Five_prime_prio5_size21 = XXXX
#        self.Deletion_D_gene_Three_prime_prio5_size21 = XXXX
#        self.Deletion_J_gene_Five_prime_prio5_size23 = XXXX
#        self.Insertion_VD_genes_Undefined_side_prio4_size31 = XXXX
#        self.DinucMarkov_VD_genes_Undefined_side_prio3_size16 = XXXX
#        self.Insertion_DJ_gene_Undefined_side_prio2_size31 = XXXX
#        self.DinucMarkov_DJ_gene_Undefined_side_prio1_size16 = XXXX
#        self.Mismatches = XXXX
#        
#        
#        self.flnBestScenarios = flnBestScenarios
#        self.mdl = IgorModel.Model(model_parms_file=model_parms_file, model_marginals_file=model_marginals_file)
#        
    def setModelFromFiles(self, model_parms_file=None, model_marginals_file=None):
        self.mdl = IgorModel(model_parms_file=model_parms_file, model_marginals_file=model_marginals_file)
        
    def setModel(self, mdl):
            self.mdl = mdl
    
    def setInputSeqsFile(self, flnInputSeqs):
        self.flnInputSeqs = flnInputSeqs
    
    def setAlignsFile(self, flnAligns):
        self.flnAligns = flnAligns
    
    # TODO: clomplete this if worth it...
    def get_BestScenariosDataFrame(self, flnBestScenarios, flnParms, flnMarginals):
        bestScen = pd.read_csv(flnBestScenarios, sep=';')
        mdl = IgorModel(model_parms_file=flnParms, model_marginals_file=flnMarginals)
        pdBestScen = bestScen.rename(mdl.parms.dictNameNickname, axis=1).set_index('seq_index')
        tmpEvents = ['vd_dinucl', 'dj_dinucl', 'Mismatches']
        tmpEvents02 = ['scenario_rank', 'scenario_proba_cond_seq']
        #strEvent  = 'v_choice'
        for strEvent in pdBestScen.keys():
            if strEvent in tmpEvents:
                pdBestScen[strEvent] = pdBestScen[strEvent].apply(lambda x : x.replace("(", "[").replace(")", "]")).apply(eval)
            elif strEvent in tmpEvents02:
                print(strEvent)
            else:
                #print strEvent
                pdBestScen[strEvent] = pdBestScen[strEvent].apply(eval)

        return pdBestScen

    # Reconstruct the naive sequence from the bestscenarios data.
    def getV_Region(self, seq_id, pdSelected):
        #seq_id=59
        strEv  = 'v_choice'
        id_V   = pdSelected.loc[seq_id][strEv]
        seq_V  = self.mdl.parms.Event_dict[strEv].loc[id_V]['value']#.values    
        name_V = self.mdl.parms.Event_dict[strEv].loc[id_V]['name']
        strEv  = 'v_3_del'
        n_v_3_del = pdSelected.loc[seq_id][strEv]
        #print(name_V, n_v_3_del, seq_V[:-n_v_3_del])
        return (seq_V[:-n_v_3_del])


    def getD_Region(self, seq_id, pdSelected):
        #seq_id=59
        strEv  = 'd_gene'
        id_D   = pdSelected.loc[seq_id][strEv]
        seq_D  = self.mdl.parms.Event_dict[strEv].loc[id_D]['value']#.values    
        name_D = self.mdl.parms.Event_dict[strEv].loc[id_D]['name']
        strEv  = 'd_5_del'
        n_d_5_del = pdSelected.loc[seq_id][strEv]
        strEv  = 'd_3_del'
        n_d_3_del = pdSelected.loc[seq_id][strEv]
        #print(seq_D,len(seq_D))
        #print(seq_id, seq_D,len(seq_D), name_D, n_d_5_del, n_d_3_del, seq_D[n_d_5_del:-n_d_3_del])
        return (seq_D[n_d_5_del:-n_d_3_del])

    def getJ_Region(self, seq_id, pdSelected):
        #seq_id=59
        strEv  = 'j_choice'
        id_J   = pdSelected.loc[seq_id][strEv]
        seq_J  = self.mdl.parms.Event_dict[strEv].loc[id_J]['value']#.values    
        name_J = self.mdl.parms.Event_dict[strEv].loc[id_J]['name']
        strEv  = 'j_5_del'
        n_j_5_del = pdSelected.loc[seq_id][strEv]
        #print(name_J, n_j_5_del, seq_J[n_j_5_del:])
        return (seq_J[n_j_5_del:])

    def getVD_Region(self, seq_id, pdSelected):
        #seq_id=59
        #vd_dinucl
        #vd_ins
        strEv      = 'vd_ins'
        id_VDins   = pdSelected.loc[seq_id][strEv]
        n_VDins  = self.mdl.parms.Event_dict[strEv].loc[id_VDins]['value']#.values    
        #name_VDins = self.mdl.parms.Event_dict[strEv].loc[id_VDins]['name']
        #print(self.mdl.parms.Event_dict[strEv] ) #.loc[id_VDins])
        #print(n_VDins)
        strEv        = 'vd_dinucl'
        id_VD_dinucl = pdSelected.loc[seq_id][strEv]
        seq_VD_dinucl  = self.mdl.parms.Event_dict[strEv].loc[id_VD_dinucl]['value'].values   
        return (''.join(seq_VD_dinucl.tolist()))

    def getDJ_Region(self, seq_id, pdSelected):
        strEv        = 'dj_dinucl'
        id_DJ_dinucl = pdSelected.loc[seq_id][strEv]
        seq_DJ_dinucl  = self.mdl.parms.Event_dict[strEv].loc[id_DJ_dinucl]['value'].values   
        return (''.join(seq_DJ_dinucl.tolist()))
    
    def get_BestScenariosNaiveSeq(self, seq_id, pdBestScen):
        self.getV_Region(seq_id, pdBestScen)
        V_Region  = self.getV_Region(seq_id, pdSelected)
        VD_Region = self.getVD_Region(seq_id, pdSelected)
        D_Region  = self.getD_Region(seq_id, pdSelected)
        DJ_Region = self.getDJ_Region(seq_id, pdSelected)
        J_Region  = self.getJ_Region(seq_id, pdSelected)
        return V_Region+VD_Region+D_Region+DJ_Region+J_Region
        

