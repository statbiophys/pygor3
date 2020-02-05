#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 16:03:42 2019

@author: alfaceor
"""

import IgorModel
import IgorSqliteDB
import IgorSqliteDBBestScenarios
import IgorAlignment_data
import numpy as np
from TemporalUtils import *
#import subpro

# IGoR run parameters
batchname = "DebugGood9"
flnIgorDB = "DebugGood9.db"
strWD="uchicago/selected_sequences/"
flnIgorIndexedSeq = strWD+"aligns/"+batchname+"_indexed_sequences.csv"

IgorSpecie    = "mouse"
IgorChain     = "tcr_beta"
IgorModelPath = "../../IGoR/models/"+IgorSpecie+"/"+IgorChain+"/"
IgorRefGenomePath = IgorModelPath+"ref_genome/"

flnVGeneTemplate = IgorRefGenomePath+"genomicVs.fasta"
flnDGeneTemplate = IgorRefGenomePath+"genomicDs.fasta"
flnJGeneTemplate = IgorRefGenomePath+"genomicJs.fasta"

flnVGeneCDR3Anchors = IgorRefGenomePath+"V_gene_CDR3_anchors.csv"
flnJGeneCDR3Anchors = IgorRefGenomePath+"J_gene_CDR3_anchors.csv"

### IGoR Alignments files
flnVAlignments = strWD+"aligns/"+batchname+"_V_alignments.csv"
flnDAlignments = strWD+"aligns/"+batchname+"_D_alignments.csv"
flnJAlignments = strWD+"aligns/"+batchname+"_J_alignments.csv"

### IGoR ouptut files
flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
flnIgorBestScenarios = strWD+batchname+"_output/best_scenarios_counts.csv"

### load IGoR sequences database
db = IgorSqliteDB.IgorSqliteDB()
db.createSqliteDB(flnIgorDB)
db.load_VDJ_Database(flnIgorIndexedSeq, \
                     flnVGeneTemplate, flnDGeneTemplate, flnJGeneTemplate, \
                     flnVAlignments, flnDAlignments, flnJAlignments)

### load IGoR model parms and marginals.
mdl = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
mdlParms = IgorModel.Model_Parms(flnModelParms) # mdl.parms
mdlMargs = IgorModel.Model_Marginals(flnModelMargs) # mdl.marginals

# load IGoR best scenarios file.
db_bs = IgorSqliteDBBestScenarios.IgorSqliteDBBestScenariosVDJ()
db_bs.createSqliteDB("newbatcho_bs.db")
db_bs.load_IgorBestScenariosVDJ_FromCSV(flnIgorBestScenarios)

#***************** BEGIN Get best scenarios with insertions 11 *****************#
sqlSelect = "SELECT * FROM IgorDBBestScenariosVDJ WHERE id_dj_ins = "+str(11)+";"
cur = db_bs.conn.cursor()
cur.execute(sqlSelect)
record_ins_11 = cur.fetchall()

#### FROM ALL THE RECORDS WITH 11 insertions get one.
seq_index = record_ins_11[0][0]
alnDataListV = db.appendList_IgorAlignments_data_By_seq_index("V", seq_index)
alnDataListD = db.appendList_IgorAlignments_data_By_seq_index("D", seq_index)
alnDataListJ = db.appendList_IgorAlignments_data_By_seq_index("J", seq_index)

print( alnDataListV[0].to_dict() ) #.deletions
print( alnDataListD[0].to_dict() )
print( alnDataListJ[0].to_dict() )
seq_index = record_ins_11[1][0]
strSeq   = db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]
print("SELECTED SEQUENCE TO ANALYZE")
print(" seq_index   : ", seq_index )
print(" sequence    : ", strSeq )
print(" seq. lenght : ", len(strSeq) )

################ NOW I WANT THE BEST CASE SCENARIOS FOR THE SEQ_INDEX

alnDataListV = db.appendList_IgorAlignments_data_By_seq_index("V", seq_index)
alnDataListD = db.appendList_IgorAlignments_data_By_seq_index("D", seq_index)
alnDataListJ = db.appendList_IgorAlignments_data_By_seq_index("J", seq_index)

print(alnDataListV[0].to_dict())
print(alnDataListD[1].to_dict())
print(alnDataListJ[0].to_dict())


import IgorBestScenarios

record_bs = db_bs.fetch_IgorBestScenariosVDJ_By_seq_index(seq_index)
for record in record_bs:
    print(record[2])
mdl = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

probs_Event = list()
probs_Scena = list()
num_DJ_Ins  = list()
num_VD_Ins  = list()
ranks_Scena = range(len(record_bs))
print(len(record_bs))
for ii in ranks_Scena:
    bs = IgorBestScenarios.IgorBestScenariosVDJ.load_FromSQLRecord(record_bs[ii])
    bs.setModel_Parms(flnModelParms)
    bs.strSeq_index  = db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]
    bs.save_scenario_fasta("scen__seqIndex_"+str(seq_index)+"_rank_"+str(ii)+".fasta")
    bs.mdl = mdl
    print(bs.scenario_rank)
    probs_Event.append( bs.get_EventProb() )
    #bs.get_ErrorProb()
    probs_Scena.append( bs.scenario_proba_cond_seq )
    num_VD_Ins.append( bs.getVD_ins())
    num_DJ_Ins.append( bs.getDJ_ins())
    
print(bs)


##### The event I think is the best one
my_bs = IgorBestScenarios.IgorBestScenariosVDJ() #.load_FromSQLRecord(record_bs[ii])
my_bs.setModel_Parms(flnModelParms)
my_bs.mdl = mdl
my_bs.seq_index = seq_index
my_bs.strSeq_index  = db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]

# V gene
pd_event  = mdl.parms.Event_dict['v_choice']
gene_name = 'TRBV17*01'
gene_id   = pd_event.loc[pd_event['name'] == gene_name ].index.values[0]
my_bs.id_v_choice = gene_id

# J gene
pd_event  = mdl.parms.Event_dict['j_choice']
gene_name = 'TRBJ2-3*01'
gene_id   = pd_event.loc[pd_event['name'] == gene_name ].index.values[0]
my_bs.id_j_choice = gene_id

# D gene
pd_event  = mdl.parms.Event_dict['d_gene']
gene_name = 'TRBD1*01'
gene_id   = pd_event.loc[pd_event['name'] == gene_name ].index.values[0]
my_bs.id_d_gene = gene_id

# v_3_del 
pd_event    = mdl.parms.Event_dict['v_3_del']
pd_event['value']
realiz_name = 2
realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
my_bs.id_v_3_del = realiz_id


# d_5_del 
pd_event    = mdl.parms.Event_dict['d_5_del']
pd_event['value']
realiz_name = 6
realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
my_bs.id_d_5_del = realiz_id


# d_3_del 
pd_event    = mdl.parms.Event_dict['d_3_del']
pd_event['value']
realiz_name = 1
realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
my_bs.id_d_3_del = realiz_id

# j_5_del 
pd_event    = mdl.parms.Event_dict['j_5_del']
pd_event['value']
realiz_name = 0
realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
my_bs.id_j_5_del = realiz_id


# vd_ins 
pd_event    = mdl.parms.Event_dict['vd_ins']
pd_event['value']
realiz_name = 1
realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
my_bs.id_vd_ins = realiz_id



# dj_ins 
pd_event    = mdl.parms.Event_dict['dj_ins']
pd_event['value']
realiz_name = 3
realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
my_bs.id_dj_ins = realiz_id

my_bs.vd_dinucl = [1]
my_bs.dj_dinucl = [3, 1, 3]

        # v_3_del = 2 
        # d_5_del = 6
        # d_3_del = 1
        # vd_ins  = 1 and should be a "C"
        # dj_ins  = 3 and should be a "TCT"

# FIXME: TMP set to zero ins and dels
my_bs.mismatcheslen = 0

my_bs.getV_Region()

print(my_bs)

my_bs.save_scenario_fasta("My_scen__seqIndex_"+str(seq_index)+"_rank_"+str(ii)+".fasta")
# so the correction to this scenario must be
# v_3_del = 2 
# d_5_del = 6
# d_3_del = 1
# vd_ins  = 1 and should be a "C"
# dj_ins  = 3 and should be a "TCT"


my_bs.get_EventProb() #1.16859699e-10

##### Plotting probabilities
import matplotlib.pyplot as plt
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18} 
plt.rc('font', **font)
plt.rc('text', usetex=True)
fig, ax = plt.subplots(2,1, figsize=(10,8))
fig.suptitle("$a = "+str(seq_index)+"$")

ax2 = ax[0].twinx()  # instantiate a second axes that shares the same x-axis


ax[0].set_ylabel("$P(E_i | R^{a}, \\theta_c)$", size=20)
ax[1].set_ylabel("$DJ_{ins}$", size=20)


ax[0].plot(ranks_Scena, probs_Scena, 'o-')
ax2.set_ylabel("$P(E_i| , \\theta_c)$", size=20)
ax2.yaxis.label.set_color('red')
ax2.spines['right'].set_color('red')
ax2.tick_params(axis='y', colors='red')
ax2.plot(ranks_Scena, probs_Event, marker='o', color='r')

#ax[1].plot(ranks_Scena, num_VD_Ins, 'o-')
ax[1].plot(ranks_Scena, num_DJ_Ins, 'o-')
ax[1].set_xlabel("$E_i$", size=18)

fig.tight_layout()
fig.savefig("ProbsEvents02.pdf")

###########  CHECK POSSIBLE SCENARIOS GIVEN THE ALIGNMENTS   #############

alnDataListV = db.appendList_IgorAlignments_data_By_seq_index("V", seq_index)
alnDataListD = db.appendList_IgorAlignments_data_By_seq_index("D", seq_index)
alnDataListJ = db.appendList_IgorAlignments_data_By_seq_index("J", seq_index)




######### CHECKING PROBABILITIES


import matplotlib.pyplot as plt
mdl.xdata.keys()
mdl.xdata['dj_ins'].plot(marker='o')
mdl.xdata['v_choice']


My_bs = IgorBestScenarios.IgorBestScenariosVDJ()

#['v_choice', 'd_gene', 'j_choice', 'v_3_del', 'd_3_del', 'd_5_del', 'j_5_del', 'vd_ins', 'dj_ins', 'vd_dinucl', 'dj_dinucl'])
