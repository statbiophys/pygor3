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
batchname = "TRbeta"
flnIgorDB = "chicagoMouse.db"
strWD="uchicago/"
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
db_bs.createSqliteDB("chicagoMouse_bs.db")
db_bs.load_IgorBestScenariosVDJ_FromCSV(flnIgorBestScenarios)

#***************** Get best scenarios with insertions 11 *****************#
sqlSelect = "SELECT * FROM IgorDBBestScenariosVDJ WHERE id_dj_ins = "+str(11)+";"
cur = db_bs.conn.cursor()
cur.execute(sqlSelect)
record_ins_11 = cur.fetchall()

#### FROM ALL THE RECORDS WITH 11 insertions get one.
seq_index = record_ins_11[0][0]
strSeq   = db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]
print("SELECTED SEQUENCE TO ANALYZE")
print(" seq_index   : ", seq_index)
print(" sequence    : ", strSeq)
print(" seq. lenght : ", len(strSeq) )




### begin VDJ ###
alnDataListVDJ = db.appendList_IgorAlignments_data_By_seq_index("V", seq_index)
db.appendList_IgorAlignments_data_By_seq_index("D", seq_index, alnDataList=alnDataListVDJ)
db.appendList_IgorAlignments_data_By_seq_index("J", seq_index, alnDataList=alnDataListVDJ)
alnDataListVDJ_pd = create_alnDataList_pandas(seq_index, strSeq, alnDataListVDJ)  
alnDataListVDJ_pd = addInsertionGaps2alnDataFrame(alnDataListVDJ_pd)
alnDataListVDJ_pd = create_alnDataList_pandas(seq_index, strSeq, alnDataListVDJ)  
alnDataListVDJ_pd = addInsertionGaps2alnDataFrame(alnDataListVDJ_pd)
writeAlignmentsFastaOnlyInsertions(alnDataListVDJ_pd, "VDJonlyIns.fasta")
### end VDJ ###

#seq_index│scenario_rank│scenario_proba_cond_seq│GeneChoice_V_gene_Undefined_side_prio7_size35│GeneChoice_J_gene_Undefined_side_prio7_size14│GeneC
#51│1│0.687667│(3)│(5)│(1)│(9)│(4)│(8)│(15)│(1)│(1)│(4)│(0,0,0,0)│(125,126,128,131,132,133,134,135,136,140,142,150,157,158,159,163,164)

### specific genes
alnDataList = db.appendList_IgorAlignments_data_By_seq_index("J", seq_index)
alnDataList_pd = create_alnDataList_pandas(seq_index, strSeq, alnDataList)  
alnDataList_pd = addInsertionGaps2alnDataFrame(alnDataList_pd)

ajam = alnDataList_pd['score'] > 20
ajam.loc[0] = True
select_alnDataList_pd = alnDataList_pd.loc[ajam ]
#select_alnDataList_pd.append(alnDataList_pd.loc[0])
#writeAlignmentsFasta(select_alnDataList_pd, "alignNoDels.fasta")
writeAlignmentsFastaOnlyInsertions(select_alnDataList_pd, "alignNoDels.fasta")

# one complete alignment
print(alnDataList_pd )
id_pd = 2
strSeqWithDels = addTemporalDelsInSeq(alnDataList_pd, id_pd)
alnDataList_pd['seq_with_ins_offset'].loc[0] = strSeqWithDels
writeAlignmentsFastaOnlyInsertions(alnDataList_pd.loc[[0,id_pd]], "example.fasta")



################# CHECKING ALIGNMENTS SCORE.
alnDataListVDJ_pd.keys()
alnDataListVDJ_pd[['description','seq_no_align']]

alnDataListVDJ_pd.loc[0]['seq_no_align']
alnDataListVDJ_pd.loc[2]['seq_no_align']

################ NOW I WANT THE BEST CASE SCENARIOS FOR THE SEQ_INDEX



#mdl.xdata['v_choice']['lbl__v_choice'=='TRBV2*01']
#
#mdl.xdata['v_choice'][{'lbl__v_choice': 'TRBV17*01'}]

alnDataListV = db.appendList_IgorAlignments_data_By_seq_index("V", seq_index)
alnDataListD = db.appendList_IgorAlignments_data_By_seq_index("D", seq_index)
alnDataListJ = db.appendList_IgorAlignments_data_By_seq_index("J", seq_index)

from TemporalUtils import *

def writeFastaPlease(seq_index, strSeq, alnDataList, flnAlignsFasta):
    alnDataList_pd = create_alnDataList_pandas(seq_index, strSeq, alnDataList)  
    alnDataList_pd = addInsertionGaps2alnDataFrame(alnDataList_pd)
    writeAlignmentsFastaOnlyInsertions(alnDataList_pd, flnAlignsFasta)


alnDataListV[0].strGene_name
writeFastaPlease(seq_index, strSeq, alnDataListJ, "id_"+str(seq_index)+"__J_onlyIns.fasta")

for alnData in alnDataListV:
    print(type(alnData))

print(len(alnDataListV))
print(len(alnDataListD))
print(len(alnDataListJ))


alnDataListVDJ[3].to_dict()

print(alnDataListVDJ[0])




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

for alnV in alnDataListV:
    print (alnV.to_dict())
    print([alnV.offset_5_p
len(alnV.strGene_seq[0:102])

alnV = alnDataListV[0]
alnD = alnDataListD[1]
alnJ = alnDataListJ[0]
alnV.to_dict()
alnD.to_dict()
alnJ.to_dict()
print("V :", alnV.offset_5_p, alnV.offset_3_p)
print("D :", alnD.offset_5_p, alnD.offset_3_p)
print("J :", alnJ.offset_5_p, alnJ.offset_3_p)
######### CHECKING PROBABILITIES


import matplotlib.pyplot as plt
mdl.xdata.keys()
mdl.xdata['dj_ins'].plot(marker='o')
mdl.xdata['v_choice']


My_bs = IgorBestScenarios.IgorBestScenariosVDJ()

#['v_choice', 'd_gene', 'j_choice', 'v_3_del', 'd_3_del', 'd_5_del', 'j_5_del', 'vd_ins', 'dj_ins', 'vd_dinucl', 'dj_dinucl'])


nuc44 = [5,-14,-14,-14,-14,2,-14,2,2,-14,-14,1,1,1,0,\
	        -14,5,-14,-14,-14,2,2,-14,-14,2,1,-14,1,1,0,\
	        -14,-14,5,-14,2,-14,2,-14,2,-14,1,1,-14,1,0,\
	        -14,-14,-14,5,2,-14,-14,2,-14,2,1,1,1,-14,0,\
	        -14,-14,2,2,1.5,-14,-12,-12,-12,-12,1,1,-13,-13,0,\
	        2,2,-14,-14,-14,1.5,-12,-12,-12,-12,-13,-13,1,1,0,\
	        -14,2,2,-14,-12,-12,1.5,-14,-12,-12,1,-13,-13,1,0,\
	        2,-14,-14,2,-12,-12,-14,1.5,-12,-12,-13,1,1,-13,0,\
	        2,-14,2,-14,-12,-12,-12,-12,1.5,-14,-13,1,-13,1,0,\
	        -14,2,-14,2,-12,-12,-12,-12,-14,1.5,1,-13,1,-13,0,\
	        -14,1,1,1,1,-13,1,-13,-13,1,0.5,-12,-12,-12,0,\
	        1,-14,1,1,1,-13,-13,1,1,-13,-12,0.5,-12,-12,0,\
	        1,1,-14,1,-13,1,-13,1,-13,1,-12,-12,0.5,-12,0,\
	        1,1,1,-14,-13,1,1,-13,1,-13,-12,-12,-12,0.5,0,\
	        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ]
nuc44arr = np.array(nuc44).reshape(15,15)

nuc44pd = pd.DataFrame(nuc44arr, columns=["A","C","G","T","R","Y","K","M","S","W","B","D","H","V","N"])
nuc44pd["nt"]= ["A","C","G","T","R","Y","K","M","S","W","B","D","H","V","N"]
#nuc44pd.set_index("nt")

nuc44pd.loc[1]

# So the things I want to solve is how to show the alignments.
# Because the simple answer will be how to check the alignments pair to pair first.
#from Bio import SeqIO
#
#Bio.SeqIO


#csvline = "2197;1;0.357694;(1);(9);(0);(7);(10);(7);(11);(0);();(9);(0,2,0,1,2,3,2,0,3);()\n"
