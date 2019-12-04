#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 05:13:57 2019

@author: alfaceor
"""


batchname = "TRbeta"
strWD="uchicago/"
IgorSpecie    = "mouse"
IgorChain     = "tcr_beta"

IgorModelPath = "../../IGoR/models/"+IgorSpecie+"/"+IgorChain+"/"
IgorRefGenomePath = IgorModelPath+"ref_genome/"


flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"

flnVGeneTemplate = IgorRefGenomePath+"genomicVs.fasta"
flnDGeneTemplate = IgorRefGenomePath+"genomicDs.fasta"
flnJGeneTemplate = IgorRefGenomePath+"genomicJs.fasta"


### load IGoR model parms and marginals.
import IgorModel
mdl = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
mdlParms = IgorModel.Model_Parms(flnModelParms) # mdl.parms
mdlMargs = IgorModel.Model_Marginals(flnModelMargs) # mdl.marginals


# 1. Plot model graph

# 2. Plot marginals
#mdl.xdata['v_choice'].plot()
import matplotlib.pyplot as plt
mdl.xdata['vd_ins'].plot()
mdl.xdata['dj_ins'].plot()
mdl.xdata['j_5_del'].plot()
mdl.xdata['j_5_del'][{'j_choice' : 6}].plot()
mdl.xdata['j_choice'][{'j_choice' : 12}]



### IGoR ouptut files

### load IGoR sequences database
import IgorSqliteDB
flnIgorDB = batchname+".db"

### IGoR Indexed_sequences
flnIgorIndexedSeq = strWD+"aligns/"+batchname+"_indexed_sequences.csv"

### IGoR Alignments files
flnVAlignments = strWD+"aligns/"+batchname+"_V_alignments.csv"
flnDAlignments = strWD+"aligns/"+batchname+"_D_alignments.csv"
flnJAlignments = strWD+"aligns/"+batchname+"_J_alignments.csv"


db = IgorSqliteDB.IgorSqliteDB()
db.createSqliteDB(flnIgorDB)
db.load_VDJ_Database(flnIgorIndexedSeq, \
                     flnVGeneTemplate, flnDGeneTemplate, flnJGeneTemplate, \
                     flnVAlignments, flnDAlignments, flnJAlignments)



import IgorSqliteDBBestScenarios

flnIgorDB_bs = batchname+"_bs.db"
!rm TRbeta_bs.db
flnIgorBestScenarios = strWD+batchname+"_output/best_scenarios_counts.csv"

db_bs = IgorSqliteDBBestScenarios.IgorSqliteDBBestScenariosVDJ()
db_bs.createSqliteDB(flnIgorDB_bs)
db_bs.load_IgorBestScenariosVDJ_FromCSV(flnIgorBestScenarios)

#***************** BEGIN Get best scenarios with insertions 11 *****************#
import numpy as np
sqlSelect = "SELECT seq_index FROM IgorDBBestScenariosVDJ WHERE id_dj_ins = "+str(11)+";"
sqlSelect = "SELECT seq_index FROM IgorDBBestScenariosVDJ WHERE id_dj_ins = "+str(11)+" AND scenario_rank = "+str(1)+";"


sqlSelect = "SELECT * FROM IgorDBBestScenariosVDJ WHERE seq_index = "+str(59)+";"


cur = db_bs.conn.cursor()
cur.execute(sqlSelect)
seq_indexes = cur.fetchall()
seq_indexes = np.unique(np.array(seq_indexes))
print(len(seq_indexes))

#sqlSelect = "SELECT * FROM IgorDBBestScenariosVDJ WHERE id_dj_ins = "+str(11)+";"
sqlSelect = "SELECT * FROM IgorDBBestScenariosVDJ WHERE id_dj_ins = "+str(11)+" AND scenario_rank = "+str(1)+";"

cur = db_bs.conn.cursor()
cur.execute(sqlSelect)
record_bs_ins_11 = cur.fetchall()

import IgorBestScenarios
import IgorAlignment_data
len(record_bs_ins_11)
set_indexes_ins_11 = set()
#indexes = list()
for record_bs in record_bs_ins_11:
    bs_djins11 = IgorBestScenarios.IgorBestScenariosVDJ.load_FromSQLRecord(record_bs)
    set_indexes_ins_11.add(bs_djins11.seq_index)
    print(bs_djins11.seq_index, bs_djins11.id_dj_ins, bs_djins11.scenario_rank)
    for strGene in ["V", "D", "J"]:
        records_aln = db.fetch_IgorAlignments_By_seq_index(strGene, bs_djins11.id_dj_ins)
        #print(strGene, len(records_aln))
        for record_aln in records_aln:
            aln = IgorAlignment_data.IgorAlignment_data.load_FromSQLRecord(record_aln)
            #print(strGene, aln.gene_id)
            #aln.gene_id
            if not len(aln.insertions) == 0:
                print(strGene+" alns ins -> ", aln.insertions)
            if not len(aln.deletions) == 0:
                print(strGene+" alns dels -> ", aln.deletions)


bs_djins11.setModel_Parms(flnModelParms)
bs_djins11.mdl = mdl
bs_djins11.strSeq_index = db.fetch_IgorIndexedSeq_By_seq_index(bs_djins11.seq_index)[1]
bs_djins11.str_scenario_fasta()

bs_djins11.getJ_ntsequence()
bs_djins11.save_scenario_fasta("aln_bs__1665__r_1.fasta")
print(bbb)

'CTCCTATGAACAGTACTTCGGTCCCGGCACCAGGCTCACGGTTTTAG'
       'GAACAGTACTTCGGTCCCGGCACCAGGCTCACGGTTTT'

'GGGACAGGGGGC'
    'GGGGGGG'
    
'GGAGCACTCGTCTATCAATATCCCAGAAGAACCATCTGTAAGAGTGGAACTTCCATGAGGATGGAGTGTCAAGCTGTGGGTTTTCAGGCAACTTCTGTAGCTTGGTATCGTCAATCGCCTCAAAAGACATTTGAACTGATAGCACTTTCTACTGTGAACTCAGCAATCAAATATGAACAAAATTTTACCCAGGAAAAATTTCCCATCAGTCATCCCAACTTATCCTTTTCATCTATGACAGTTTTAAATGCATATCTTGAAGACAGAGGCTTATATCTCTGTGGTGCTAGGGA'    
'GGAGCACTCGTCTATCAATATCCCAGAAGAACCATCTGTAAGAGTGGAACTTCCATGAGGATGGAGTGTCAAGCTGTGGGTTTTCAGGCAACTTCTGTAGCTTGGTATCGTCAATCGCCTCAAAAGACATTTGAACTGATAGCACTTTCTACTGTGAACTCAGCAATCAAATATGAACAAAATTTTACCCAGGAAAAATTTCCCATCAGTCATCCCAACTTATCCTTTTCATCTATGACAGTTTTAAATGCATATCTTGAAGACAGAGGCTTATATCTCTGTGGTGC'

# How many sequence I have with indels.
sqlSelect = "SELECT * FROM IgorVAlignments;" # WHERE NOT insertions = [];"
cur = db.conn.cursor()
cur.execute(sqlSelect)
records = cur.fetchall()
set_indexes = set()
for record in records:
    aln = IgorAlignment_data.IgorAlignment_data.load_FromSQLRecord(record)
    if not len(aln.insertions) == 0:
        print(aln.seq_index, aln.insertions, aln.deletions)
        set_indexes.add(aln.seq_index)

set_indexes.intersection(set_indexes_ins_11)




aln.to_dict()

db.fetch_IgorIndexedSeq_By_seq_index(59)


bs_djins11.get_EventProb()

aaa = bs_djins11.to_dict()
aaa = bs_djins11.to_dict_names()
aaa = bs_djins11.to_dict_ntsequences()

aaa['v_choice']

aaa['v_3_del']
aaa['d_gene']
aaa['d_5_del']
aaa['d_3_del']
aaa['j_choice']
aaa['j_5_del']
aaa['dj_ins']
aaa['dj_dinucl']
aaa['mismatches']

mdl.parms.Event_dict['d_gene']


################### Mixcr sequence.

dict_scenario_mixcr = {
"v_choice" : 23,
"j_choice" : 13,
"d_gene" : 1,
"v_3_del" : 4,
"j_5_del" : 11,
"d_5_del" : 10,
"d_3_del" : 5,
"vd_ins" : 5,
"vd_dinucl" : "[1,2,2,2,2]",
"dj_ins" : 0,
"dj_dinucl" : "[]"
}

bs_mixcr = IgorBestScenarios.IgorBestScenariosVDJ.load_FromDict(dict_scenario_mixcr)
bs_mixcr.setModel_Parms(flnModelParms)
bs_mixcr.mdl = mdl

bs_mixcr.seq_index = bs_djins11.seq_index
bs_mixcr.strSeq_index = db.fetch_IgorIndexedSeq_By_seq_index(bs_mixcr.seq_index)[1]
bs_mixcr.save_scenario_fasta("aln_bs__1665__mixcr.fasta")


bs_mixcr.get_EventProb()

zzz = bs_mixcr.get_DictNicknameProbs()
zzz['v_choice']
zzz['d_gene']
zzz['j_choice']

zzz['v_3_del']
zzz['j_5_del']
zzz['d_5_del']
zzz['d_3_del']
zzz['vd_ins']

zzz['dj_ins']

mdl.parms.Event_dict['d_3_del']
mdl.xdata['d_3_del'][{'d_gene' : 1, 'd_5_del' : 11}].plot(marker='o')

mdl.xdata['d_3_del'][{'d_gene' : 1, 'd_3_del' : 5}].plot(marker='o')


bs_mixcr.getJ_Region()
bs_mixcr.getJ_5_dels()





