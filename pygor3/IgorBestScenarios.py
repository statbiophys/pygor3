#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 10:18:45 2019

@author: alfaceor
"""

from .IgorIO import *

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
        #self.IgorModelParms
        self.strSeq_index = ""
        
#        # To fetch data from datase connection to a database
#        self.IgorDB = ""

    def __str__(self):
        return str(self.to_dict())
    
    def setModel_Parms(self, flnModelParms):
        self.flnModelParms = flnModelParms
        self.mdlParms   = IgorModel_Parms(model_parms_file = self.flnModelParms)

    def to_dict(self):
        dictBestScenario       =  {
            "seq_index"     : self.seq_index     , \
            "scenario_rank" : self.scenario_rank , \
            "scenario_proba_cond_seq" : self.scenario_proba_cond_seq , \
            "v_choice"      : self.id_v_choice   , \
            "j_choice"      : self.id_j_choice   , \
            "d_gene"        : self.id_d_gene     , \
            "v_3_del"       : self.id_v_3_del    , \
            "d_5_del"       : self.id_d_5_del    , \
            "d_3_del"       : self.id_d_3_del    , \
            "j_5_del"       : self.id_j_5_del    , \
            "vd_ins"        : self.id_vd_ins     , \
            "vd_dinucl"     : self.vd_dinucl     , \
            "dj_ins"        : self.id_dj_ins     , \
            "dj_dinucl"     : self.dj_dinucl     , \
            "mismatches"    : self.mismatches    , \
            "mismatcheslen" : self.mismatcheslen
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

    def to_dict_values(self):
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
            "vd_dinucl": self.getVD_Region(), \
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

    # FIXME: CHANGE FOR A BETTER NAME
    def to_dict_export(self):
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
        #seq_index;scenario_rank;scenario_proba_cond_seq;GeneChoice_V_gene_Undefined_side_prio7_size35;GeneChoice_J_gene_Undefined_side_prio7_size14;GeneChoice_D_gene_Undefined_side_prio6_size2;Deletion_V_gene_Three_prime_prio5_size21;Deletion_D_gene_Five_prime_prio5_size21;Deletion_D_gene_Three_prime_prio5_size21;Deletion_J_gene_Five_prime_prio5_size23;Insertion_VD_genes_Undefined_side_prio4_size31;DinucMarkov_VD_genes_Undefined_side_prio3_size16;Insertion_DJ_gene_Undefined_side_prio2_size31;DinucMarkov_DJ_gene_Undefined_side_prio1_size16;Mismatches
        cls = IgorBestScenariosVDJ()
        linesplit = line.split(delimiter)
        linesplit = line.split(";")
        for ii in range(len(linesplit)):
            # TODO: find a better way to do this, if is a list keep it as list
            if (ii in [ 11, 13, 14 ]):
                linesplit[ii] = linesplit[ii]
            else:
                linesplit[ii] = linesplit[ii].replace("(", "").replace(")", "")

        try:
            #1;1;0.123596;(14);(9);(1);(4);(8);(7);(11);(0);();(9);(0,2,0,1,2,3,2,0,0);(122,123,124)
            cls.seq_index      = int(linesplit[ 0])
            cls.scenario_rank  = int(linesplit[ 1])
            cls.scenario_proba_cond_seq = float(linesplit[2])
            #print(linesplit[ 3], type(linesplit[ 3]), len(linesplit[ 3]))
            cls.id_v_choice    = int(linesplit[ 3])
            cls.id_j_choice    = int(linesplit[ 4])
            cls.id_d_gene      = int(linesplit[ 5])
            cls.id_v_3_del     = int(linesplit[ 6])
            cls.id_d_5_del     = int(linesplit[ 7])
            cls.id_d_3_del     = int(linesplit[ 8])
            cls.id_j_5_del     = int(linesplit[ 9])
            cls.id_vd_ins      = int(linesplit[10])
            cls.vd_dinucl      = eval(linesplit[11].replace("(","[").replace(")","]"))
            cls.id_dj_ins      = int(linesplit[12])
            cls.dj_dinucl      = eval(linesplit[13].replace("(","[").replace(")","]"))
            cls.mismatches     = eval(linesplit[14].replace("(","[").replace(")","]"))
            cls.mismatcheslen  = int(len(cls.mismatches))
            
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
            cls.id_v_choice    = int(dictBestScenarios["v_choice"])
            cls.id_j_choice    = int(dictBestScenarios["j_choice"])
            cls.id_d_gene      = int(dictBestScenarios["d_gene"])
            cls.id_v_3_del     = int(dictBestScenarios["v_3_del"])
            cls.id_d_5_del     = int(dictBestScenarios["d_5_del"])
            cls.id_d_3_del     = int(dictBestScenarios["d_3_del"])
            cls.id_j_5_del     = int(dictBestScenarios["j_5_del"])
            cls.id_vd_ins      = int(dictBestScenarios["vd_ins"])
            cls.vd_dinucl      = eval(dictBestScenarios["vd_dinucl"])
            cls.id_dj_ins      = int(dictBestScenarios["dj_ins"])
            cls.dj_dinucl      = eval(dictBestScenarios["dj_dinucl"])
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
            cls.seq_index      = int(sqlRecordBestScenarios[ 0])
            cls.scenario_rank  = int(sqlRecordBestScenarios[ 1])
            cls.scenario_proba_cond_seq = float(sqlRecordBestScenarios[2])
            cls.id_v_choice    = int(sqlRecordBestScenarios[ 3])
            cls.id_j_choice    = int(sqlRecordBestScenarios[ 4])
            cls.id_d_gene      = int(sqlRecordBestScenarios[ 5])
            cls.id_v_3_del     = int(sqlRecordBestScenarios[ 6])
            cls.id_d_5_del     = int(sqlRecordBestScenarios[ 7])
            cls.id_d_3_del     = int(sqlRecordBestScenarios[ 8])
            cls.id_j_5_del     = int(sqlRecordBestScenarios[ 9])
            cls.id_vd_ins      = int(sqlRecordBestScenarios[10])
            cls.vd_dinucl      = eval(sqlRecordBestScenarios[11])
            cls.id_dj_ins      = int(sqlRecordBestScenarios[12])
            cls.dj_dinucl      = eval(sqlRecordBestScenarios[13])
            cls.mismatches     = eval(sqlRecordBestScenarios[14])
            cls.mismatcheslen  = int(sqlRecordBestScenarios[15])

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
        cls = IgorBestScenariosVDJ() #.load_FromSQLRecord(record_bs[ii])
        #cls.setModel_Parms(flnModelParms)
        cls.mdl = mdl
        cls.seq_index = seq_index #59
        cls.strSeq_index  = strSeq_index #db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]
        Event_GeneChoice = ['v_choice', 'j_choice', 'd_gene']
        Event_Deletions  = ['v_3_del', 'd_5_del', 'd_3_del', 'j_5_del']
        Event_Insertions = ['vd_ins', 'dj_ins']
        Event_Dinucl     = ['vd_dinucl', 'dj_dinucl']
        for event_nickname in scenario_dict.keys():
            if event_nickname in Event_GeneChoice:
                pd_event  = cls.mdl.parms.Event_dict[event_nickname]
                gene_name =   scenario_dict[event_nickname] # 'TRBV17*01'
                gene_id   = pd_event.loc[pd_event['name'] == gene_name ].index.values[0]
                if event_nickname == 'v_choice':
                    cls.id_v_choice = gene_id
                elif event_nickname == 'j_choice':
                    cls.id_j_choice = gene_id
                elif event_nickname == 'd_gene':
                    cls.id_d_gene = gene_id
                else:
                    print( "Something vey bad happen with "+str(scenario_dict[event_nickname]) )
                
            elif event_nickname in Event_Deletions:
                pd_event    = cls.mdl.parms.Event_dict[event_nickname]
                realiz_name = scenario_dict[event_nickname]
                realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
                if event_nickname == 'v_3_del':
                    cls.id_v_3_del = realiz_id
                elif event_nickname == 'd_5_del':
                    cls.id_d_5_del = realiz_id
                elif event_nickname == 'd_3_del':
                    cls.id_d_3_del = realiz_id
                elif event_nickname == 'j_5_del':
                    cls.id_j_5_del = realiz_id
                else:
                    print( "Something vey bad happen with "+str(scenario_dict[event_nickname]) )
                    
            elif event_nickname in Event_Insertions:
                pd_event    = cls.mdl.parms.Event_dict[event_nickname]
                realiz_name = scenario_dict[event_nickname]
                realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
                if event_nickname == 'v_3_del':
                    cls.id_vd_ins = realiz_id
                elif event_nickname == 'd_5_del':
                    cls.id_dj_ins = realiz_id
                else:
                    print( "Something vey bad happen with "+str(scenario_dict[event_nickname]) )
                
            elif event_nickname in Event_Dinucl:
                if event_nickname == 'vd_dinucl':
                    pd_event     = cls.mdl.parms.Event_dict[event_nickname]
                    str_sequence = scenario_dict[event_nickname]
                    list_id_seq  = list()
                    for str_nt in str_sequence:
                        realiz_name = str_nt
                        realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
                        list_id_seq.append(realiz_id)
                    cls.vd_dinucl = list_id_seq
                
                elif event_nickname == 'dj_dinucl':
                    pd_event     = cls.mdl.parms.Event_dict[event_nickname]
                    str_sequence = scenario_dict[event_nickname]
                    list_id_seq  = list()
                    for str_nt in str_sequence:
                        realiz_name = str_nt
                        realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
                        list_id_seq.append(realiz_id)
                    cls.dj_dinucl = list_id_seq
                
                else:
                    print("Something wrong with " + str(event_nickname) )

            elif event_nickname in ['mismatches']:
                if isinstance(scenario_dict[event_nickname], list):
                    cls.mismatches = scenario_dict[event_nickname]
                    cls.mismatcheslen = len(cls.mismatches)
                else:
                    print("mismatches aren't list")
            else:
                print( "Something bad happen!" )
        
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
        strScenarioFasta = strScenarioFasta + ">"+str(self.seq_index)+", rank: "+str(self.scenario_rank)+ ", prob: "+str(self.scenario_proba_cond_seq)+"\n"
        strScenarioFasta = strScenarioFasta + self.strSeq_index  + "\n" 
        strScenarioFasta = strScenarioFasta + self.getV_fasta()  + "\n" 
        strScenarioFasta = strScenarioFasta + self.getVD_fasta() + "\n" 
        strScenarioFasta = strScenarioFasta + self.getD_fasta()  + "\n"
        # FIXME: I THINK getDJ_fasta() is wrong (inverse order maybe)
        strScenarioFasta = strScenarioFasta + self.getDJ_fasta() + "\n" 
        strScenarioFasta = strScenarioFasta + self.getJ_fasta()  + "\n" 
        
        # FIXME: TEMPORAL 
        strScenarioFasta = strScenarioFasta + "> v_choice\n"
        strScenarioFasta = strScenarioFasta + self.getV_ntsequence() +"\n"
        
        strScenarioFasta = strScenarioFasta + "> d_gene\n"
        strScenarioFasta = strScenarioFasta + self.getD_ntsequence() +"\n"
        
        strScenarioFasta = strScenarioFasta + "> j_choice\n"
        strScenarioFasta = strScenarioFasta + self.getJ_ntsequence() +"\n"
        
        return strScenarioFasta

    
    #### V region methods
    def getV_fasta(self):
        strV_fasta = ""
        strV_fasta = strV_fasta + ">"+str(self.id_v_choice)+": "+ self.getV_gene_name()+ ", dels 3' = "+str(self.getV_3_dels()) + "\n"
        strV_fasta = strV_fasta + self.getV_Region()  + "\n"
        return strV_fasta
    
    def getV_gene_name(self):
        strEv  = 'v_choice'
        name_V = self.mdlParms.Event_dict[strEv].loc[self.id_v_choice]['name']
        return name_V
    
    def getV_ntsequence(self):
        strEv  = 'v_choice'
        seq_V  = self.mdlParms.Event_dict[strEv].loc[self.id_v_choice]['value']
        return seq_V
    
    def getV_3_dels(self):
        strEv  = 'v_3_del'
        n_v_3_del = self.mdlParms.Event_dict[strEv].loc[self.id_v_3_del]['value']
        return n_v_3_del
    
    def getV_Region(self):
        #seq_id=59
        strEv  = 'v_choice'
        seq_V  = self.mdlParms.Event_dict[strEv].loc[self.id_v_choice]['value']
        n_v_3_del = self.getV_3_dels()
        if n_v_3_del == 0 :
            return (seq_V)
        # FIXME: ADD palindromic insertions
        elif n_v_3_del < 0 :
            return (seq_V+'X'*n_v_3_del)
        else:
            return (seq_V[:-n_v_3_del])
   
    #### J region methods
    def getJ_fasta(self):
        strJ_fasta = ""
        strJ_fasta = strJ_fasta + ">"+str(self.id_j_choice)+": "+ self.getJ_gene_name() + ", dels 5' = "+str(self.getJ_5_dels()) + "\n"
        strJ_fasta = strJ_fasta + self.getJ_Region()  + "\n"
        return strJ_fasta
    
    def getJ_gene_name(self):
        strEv  = 'j_choice'
        name_J = self.mdlParms.Event_dict[strEv].loc[self.id_j_choice]['name']
        return name_J
    
    def getJ_ntsequence(self):
        strEv  = 'j_choice'
        seq_J  = self.mdlParms.Event_dict[strEv].loc[self.id_j_choice]['value']
        return seq_J
    
    def getJ_5_dels(self):
        strEv  = 'j_5_del'
        n_j_5_del = self.mdlParms.Event_dict[strEv].loc[self.id_j_5_del]['value']
        return n_j_5_del
    
    def getJ_Region(self):
        #seq_id=59
        strEv  = 'j_choice'
        seq_J  = self.mdlParms.Event_dict[strEv].loc[self.id_j_choice]['value']#.values            
        n_j_5_del = self.getJ_5_dels()
        if n_j_5_del == 0 :
            return (seq_J)
        # FIXME: ADD palindromic insertions
        elif n_j_5_del < 0 :
            return (seq_J+'X'*n_j_5_del)
        else:
            return (seq_J[n_j_5_del:])
    
    
    #### D region methods
    def getD_fasta(self):
        strD_fasta = ""
        strD_fasta = strD_fasta + ">"+str(self.id_d_gene)+": "+ self.getD_gene_name() \
                + ", "+str(self.id_d_5_del)+" dels 5' = "+str(self.getD_5_dels()) \
                + ", "+str(self.id_d_3_del)+" dels 3' = "+str(self.getD_3_dels()) + "\n"
        strD_fasta = strD_fasta + self.getD_Region()  + "\n"
        return strD_fasta
    
    def getD_gene_name(self):
        strEv  = 'd_gene'
        name_D = self.mdlParms.Event_dict[strEv].loc[self.id_d_gene]['name']
        return name_D
    
    def getD_ntsequence(self):
        strEv  = 'd_gene'
        seq_D  = self.mdlParms.Event_dict[strEv].loc[self.id_d_gene]['value']
        return seq_D
    
    def getD_5_dels(self):
        strEv  = 'd_5_del'
        n_d_5_del = self.mdlParms.Event_dict[strEv].loc[self.id_d_5_del]['value']
        return n_d_5_del
    
    def getD_3_dels(self):
        strEv  = 'd_3_del'
        n_d_3_del = self.mdlParms.Event_dict[strEv].loc[self.id_d_3_del]['value']
        return n_d_3_del
    
    def getD_Region(self):
        strEv  = 'd_gene'
        seq_D  = self.mdlParms.Event_dict[strEv].loc[self.id_d_gene]['value']#.values            
        n_d_5_del = self.getD_5_dels()
        n_d_3_del = self.getD_3_dels()
        if n_d_3_del == 0 :
            return (seq_D[n_d_5_del:])
        # FIXME: ADD palindromic insertions
        elif n_d_3_del < 0 :
            return (seq_D[n_d_5_del:]+'X'*n_d_3_del)
        else:
            return (seq_D[n_d_5_del:-n_d_3_del])
    
    
    #### VD region methods
    def getVD_fasta(self):
        strVD_fasta = ""
        strVD_fasta = strVD_fasta + ">"+str(self.id_vd_ins)+", VD insertions = "+str(self.getVD_ins()) + "\n"
        strVD_fasta = strVD_fasta + self.getVD_Region()  + "\n"
        return strVD_fasta
    
    def getVD_ins(self):
        strEv  = 'vd_ins'
        n_vd_ins = self.mdlParms.Event_dict[strEv].loc[self.id_vd_ins]['value']
        return n_vd_ins
    
    def getVD_Region(self):
        strEv        = 'vd_dinucl'
        seq_VD_dinucl  = self.mdlParms.Event_dict[strEv].loc[self.vd_dinucl]['value'].values   
        return (''.join(seq_VD_dinucl.tolist())) 
    
    #### DJ region methods
    def getDJ_fasta(self):
        strDJ_fasta = ""
        strDJ_fasta = strDJ_fasta + ">"+str(self.id_dj_ins)+", DJ insertions = "+str(self.getDJ_ins()) + "\n"
        strDJ_fasta = strDJ_fasta + self.getDJ_Region()  + "\n"
        return strDJ_fasta
    
    def getDJ_ins(self):
        strEv  = 'dj_ins'
        n_dj_ins = self.mdlParms.Event_dict[strEv].loc[self.id_dj_ins]['value']
        return n_dj_ins
    
    def getDJ_Region(self):
        strEv        = 'dj_dinucl'
        seq_DJ_dinucl  = self.mdlParms.Event_dict[strEv].loc[self.dj_dinucl]['value'].values   
        return (''.join(seq_DJ_dinucl.tolist())) 
    
    # FIXME: CHANGE NAME TO get_ScenarioProb.
    def get_EventProb(self):
        ###### v_choice
        strEvent = 'v_choice'
        da_event = self.mdl.xdata[strEvent]
        p_V = da_event[{strEvent: self.id_v_choice}]
        #print("p_V = ", p_V)
        #p_V = da_event.where(da_event['lbl__'+strEvent] == bs.getV_gene_name() ,drop=True)
        
        ###### j_choice
        strEvent = 'j_choice'
        da_event = self.mdl.xdata[strEvent]
        p_J = da_event[{strEvent: self.id_j_choice}]
        #print("p_J = ", p_J)
        #p_J = da_event.where(da_event['lbl__'+strEvent] == bs.getJ_gene_name() ,drop=True)
        
        ###### d_gene
        strEvent = 'd_gene'
        da_event = self.mdl.xdata[strEvent]
        p_DgJ = da_event[{'d_gene':self.id_d_gene, 'j_choice':self.id_j_choice}]
        #print("p_DgJ = ", p_DgJ)
        #p_DgJ
        
        
        ###### v_3_del
        strEvent = 'v_3_del'
        da_event = self.mdl.xdata[strEvent]
        p_V_3_del = da_event[{'v_3_del':self.id_v_3_del, 'v_choice':self.id_v_choice}]
        #print("p_V_3_del = ", p_V_3_del)
        #p_V_3_del
        
        
        ###### j_5_del
        strEvent = 'j_5_del'
        da_event = self.mdl.xdata[strEvent]
        p_J_5_del = da_event[{'j_5_del':self.id_j_5_del, 'j_choice':self.id_j_choice}]
        #print("p_J_5_del = ", p_J_5_del)
        #p_J_5_del
        
        ###### d_5_del
        strEvent = 'd_5_del'
        da_event = self.mdl.xdata[strEvent]
        p_D_5_del = da_event[{'d_5_del':self.id_d_5_del, 'd_gene':self.id_d_gene}]
        #print("p_D_5_del = ", p_D_5_del)
        #p_D_5_del
        
        ###### d_3_del
        strEvent = 'd_3_del'
        da_event = self.mdl.xdata[strEvent]
        p_D_3_del = da_event[{'d_3_del':self.id_d_3_del, 'd_5_del':self.id_d_5_del, 'd_gene':self.id_d_gene}]
        #print("p_D_3_del = ", p_D_3_del)
        #p_D_3_del
        
        ###### vd_ins
        strEvent = 'vd_ins'
        da_event = self.mdl.xdata[strEvent]
        p_VD_ins = da_event[{'vd_ins':self.id_vd_ins}]
        #print("p_VD_ins = ", p_VD_ins)
        
        
        ###### vd_dinucl
        strEvent = 'vd_dinucl'
        da_event = self.mdl.xdata[strEvent]
        # Get the last nucleotide of V region (after deletions)
        
        str_prev_nt = self.getV_Region()[-1]
        pd_tmp = self.mdl.parms.Event_dict[strEvent]
        prev_nt = pd_tmp.loc[pd_tmp['value'] == str_prev_nt].index.values[0]
        
        # for each nucleotide on inserted list
        Num_nt = 4 # 4 nucleotides A, C, G, T 
        p_VD_dinucl = 1
        for curr_nt in self.vd_dinucl:
            id_dinucl   = prev_nt*Num_nt + curr_nt
            prob_tmp    = da_event[{'vd_dinucl':id_dinucl}]
            p_VD_dinucl = p_VD_dinucl * prob_tmp
            #print(prev_nt, curr_nt, id_dinucl, prob_tmp, p_VD_dinucl)
            prev_nt = curr_nt
        
        #
        #print("p_VD_dinucl = ", p_VD_dinucl)
        
        ###### dj_ins
        strEvent = 'dj_ins'
        da_event = self.mdl.xdata[strEvent]
        p_DJ_ins = da_event[{'dj_ins':self.id_dj_ins}]
        #print("p_DJ_ins = ", p_DJ_ins)
        
        ###### dj_dinucl
        strEvent = 'dj_dinucl'
        da_event = self.mdl.xdata[strEvent]
        # Get the last nucleotide of V region (after deletions)
        
#        str_prev_nt = (self.getV_Region() + self.getVD_Region() + self.getD_Region() )[-1]
        str_prev_nt = (self.getJ_Region() )[0]
        pd_tmp  = self.mdl.parms.Event_dict[strEvent]
        prev_nt = pd_tmp.loc[pd_tmp['value'] == str_prev_nt].index.values[0]
        #print("prev_nt : ", prev_nt)
        
        # for each nucleotide on inserted list
        Num_nt = 4 # 4 nucleotides A, C, G, T 
        p_DJ_dinucl = 1
        #self.dj_dinucl = self.dj_dinucl[::-1]
        for curr_nt in self.dj_dinucl:
            id_dinucl   = prev_nt*Num_nt + curr_nt
            prob_tmp    = da_event[{'dj_dinucl':id_dinucl}]
            p_DJ_dinucl = p_DJ_dinucl * prob_tmp
            #print(prev_nt, curr_nt, id_dinucl, prob_tmp, p_DJ_dinucl)
            prev_nt = curr_nt
        
        #print("p_DJ_dinucl = ", p_DJ_dinucl)
        
        
        p_vecE = p_V*p_J*p_DgJ*p_V_3_del*p_J_5_del*p_D_5_del*p_D_3_del*\
                p_VD_ins*p_VD_dinucl*p_DJ_ins*p_DJ_dinucl
        
        return p_vecE.values
    
    def get_DictNicknameProbs(self):
        dictNicknameProbs       =  dict()
        {
            "v_choice"      : self.id_v_choice   , \
            "j_choice"      : self.id_j_choice   , \
            "d_gene"        : self.id_d_gene     , \
            "v_3_del"       : self.id_v_3_del    , \
            "d_5_del"       : self.id_d_5_del    , \
            "d_3_del"       : self.id_d_3_del    , \
            "j_5_del"       : self.id_j_5_del    , \
            "vd_ins"        : self.id_vd_ins     , \
            "vd_dinucl"     : self.vd_dinucl     , \
            "dj_ins"        : self.id_dj_ins     , \
            "dj_dinucl"     : self.dj_dinucl     , \
            "mismatches"    : self.mismatches    , \
            "mismatcheslen" : self.mismatcheslen
        }
        ###### v_choice
        strEvent = 'v_choice'
        da_event = self.mdl.xdata[strEvent]
        p_V = da_event[{strEvent: self.id_v_choice}]
        dictNicknameProbs[strEvent] = p_V 
        #print("p_V = ", p_V)
        #p_V = da_event.where(da_event['lbl__'+strEvent] == bs.getV_gene_name() ,drop=True)
        
        ###### j_choice
        strEvent = 'j_choice'
        da_event = self.mdl.xdata[strEvent]
        p_J = da_event[{strEvent: self.id_j_choice}]
        dictNicknameProbs[strEvent] = p_J 
        #print("p_J = ", p_J)
        #p_J = da_event.where(da_event['lbl__'+strEvent] == bs.getJ_gene_name() ,drop=True)
        
        ###### d_gene
        strEvent = 'd_gene'
        da_event = self.mdl.xdata[strEvent]
        p_DgJ = da_event[{'d_gene':self.id_d_gene, 'j_choice':self.id_j_choice}]
        dictNicknameProbs[strEvent] = p_DgJ
        #print("p_DgJ = ", p_DgJ)
        #p_DgJ
        
        
        ###### v_3_del
        strEvent = 'v_3_del'
        da_event = self.mdl.xdata[strEvent]
        p_V_3_del = da_event[{'v_3_del':self.id_v_3_del, 'v_choice':self.id_v_choice}]
        dictNicknameProbs[strEvent] = p_V_3_del
        #print("p_V_3_del = ", p_V_3_del)
        #p_V_3_del
        
        
        ###### j_5_del
        strEvent = 'j_5_del'
        da_event = self.mdl.xdata[strEvent]
        p_J_5_del = da_event[{'j_5_del':self.id_j_5_del, 'j_choice':self.id_j_choice}]
        dictNicknameProbs[strEvent] = p_J_5_del
        
        ###### d_5_del
        strEvent = 'd_5_del'
        da_event = self.mdl.xdata[strEvent]
        p_D_5_del = da_event[{'d_5_del':self.id_d_5_del, 'd_gene':self.id_d_gene}]
        dictNicknameProbs[strEvent] = p_D_5_del
        
        ###### d_3_del
        strEvent = 'd_3_del'
        da_event = self.mdl.xdata[strEvent]
        p_D_3_del = da_event[{'d_3_del':self.id_d_3_del, 'd_5_del':self.id_d_5_del, 'd_gene':self.id_d_gene}]
        dictNicknameProbs[strEvent] = p_D_3_del
        
        ###### vd_ins
        strEvent = 'vd_ins'
        da_event = self.mdl.xdata[strEvent]
        p_VD_ins = da_event[{'vd_ins':self.id_vd_ins}]
        dictNicknameProbs[strEvent] = p_VD_ins
        
        
        ###### vd_dinucl
        strEvent = 'vd_dinucl'
        da_event = self.mdl.xdata[strEvent]
        # Get the last nucleotide of V region (after deletions)
        
        str_prev_nt = self.getV_Region()[-1]
        pd_tmp = self.mdl.parms.Event_dict[strEvent]
        prev_nt = pd_tmp.loc[pd_tmp['value'] == str_prev_nt].index.values[0]
        
        # for each nucleotide on inserted list
        Num_nt = 4 # 4 nucleotides A, C, G, T 
        p_VD_dinucl = 1
        for curr_nt in self.vd_dinucl:
            id_dinucl   = prev_nt*Num_nt + curr_nt
            prob_tmp    = da_event[{'vd_dinucl':id_dinucl}]
            p_VD_dinucl = p_VD_dinucl * prob_tmp
            #print(prev_nt, curr_nt, id_dinucl, prob_tmp, p_VD_dinucl)
            prev_nt = curr_nt
        
        dictNicknameProbs[strEvent] = p_VD_dinucl
        
        ###### dj_ins
        strEvent = 'dj_ins'
        da_event = self.mdl.xdata[strEvent]
        p_DJ_ins = da_event[{'dj_ins':self.id_dj_ins}]
        dictNicknameProbs[strEvent] = p_DJ_ins
        
        ###### dj_dinucl
        strEvent = 'dj_dinucl'
        da_event = self.mdl.xdata[strEvent]
        # Get the last nucleotide of V region (after deletions)
        
        str_prev_nt = (self.getV_Region() + self.getVD_Region() + self.getD_Region() )[-1]
        pd_tmp  = self.mdl.parms.Event_dict[strEvent]
        prev_nt = pd_tmp.loc[pd_tmp['value'] == str_prev_nt].index.values[0]
        #print("prev_nt : ", prev_nt)
        
        # for each nucleotide on inserted list
        Num_nt = 4 # 4 nucleotides A, C, G, T 
        p_DJ_dinucl = 1
        for curr_nt in self.dj_dinucl:
            id_dinucl   = prev_nt*Num_nt + curr_nt
            prob_tmp    = da_event[{'dj_dinucl':id_dinucl}]
            p_DJ_dinucl = p_DJ_dinucl * prob_tmp
            #print(prev_nt, curr_nt, id_dinucl, prob_tmp, p_DJ_dinucl)
            prev_nt = curr_nt
        
        dictNicknameProbs[strEvent] = p_DJ_dinucl
        
        return dictNicknameProbs

    def get_ErrorProb(self):
        r = float(self.mdl.parms.ErrorRate['SingleErrorRate'])
        L = len(self.strSeq_index)
        print("error rate: ", r, "n mismatches", self.mismatcheslen)
        #return r**(self.mismatcheslen)
        return (r/3)**(self.mismatcheslen) # * (1-r)**( L - self.mismatcheslen)
    
