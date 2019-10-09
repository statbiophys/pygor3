#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 16:52:51 2019

@author: alfaceor
"""
import IgorModel
import pandas as pd
import numpy as np


class IgorBestScenarios:
    """
    Class to load Best scenarios file and attach the model parms and marginals
    """
    def __init__(self, flnBestScenarios, model_parms_file=None, model_marginals_file=None):
        self.flnBestScenarios = flnBestScenarios
        self.mdl = IgorModel.Model(model_parms_file=model_parms_file, model_marginals_file=model_marginals_file)
        
    def setModelFromFiles(self, model_parms_file=None, model_marginals_file=None):
        self.mdl = IgorModel.Model(model_parms_file=model_parms_file, model_marginals_file=model_marginals_file)
        
    def setModel(self, mdl):
            self.mdl = mdl
    
    def setInputSeqsFile(self, flnInputSeqs):
        self.flnInputSeqs = flnInputSeqs
    
    def setAlignsFile(self, flnAligns):
        self.flnAligns = flnAligns
    
    def get_BestScenariosDataFrame(self, flnBestScenarios, flnParms, flnMarginals):
        bestScen = pd.read_csv(flnBestScenarios, sep=';')
        mdl = IgorModel.Model(model_parms_file=flnParms, model_marginals_file=flnMarginals)
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
        

