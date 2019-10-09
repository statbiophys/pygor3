#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 16:03:42 2019

@author: alfaceor
"""

import IgorSqliteDB
import IgorAlignment_data
import numpy as np
from TemporalUtils import *
#import subpro

batchname = "TRbeta"
flnIgorDB = "chicagoMouse.db"
strWD="uchicago/"
flnIgorIndexedSeq = strWD+"aligns/"+batchname+"_indexed_sequences.csv"

IgorSpecie    = "mouse"
IgorChain     = "tcr_beta"
IgorModelPath = "../../IGoR/models/"+IgorSpecie+"/"+IgorChain+"/ref_genome/"

flnVGeneTemplate = IgorModelPath+"genomicVs.fasta"
flnDGeneTemplate = IgorModelPath+"genomicDs.fasta"
flnJGeneTemplate = IgorModelPath+"genomicJs.fasta"

flnVGeneCDR3Anchors = IgorModelPath+"V_gene_CDR3_anchors.csv"
flnJGeneCDR3Anchors = IgorModelPath+"J_gene_CDR3_anchors.csv"

flnVAlignments = strWD+"aligns/"+batchname+"_V_alignments.csv"
flnDAlignments = strWD+"aligns/"+batchname+"_D_alignments.csv"
flnJAlignments = strWD+"aligns/"+batchname+"_J_alignments.csv"

### load IGoR database
db = IgorSqliteDB.IgorSqliteDB()
db.createSqliteDB(flnIgorDB)
db.load_VDJ_Database(flnIgorIndexedSeq, flnVGeneTemplate, flnDGeneTemplate, flnJGeneTemplate, flnVAlignments, flnDAlignments, flnJAlignments)

#sqlSelect = "SELECT * FROM Igor"+strGene.upper()+"Alignments WHERE seq_index=="+str(seq_index)+" ORDER BY score DESC"
sqlSelect = "SELECT * FROM IgorVAlignments WHERE (NOT insertions='[]') AND (NOT deletions='[]') " #ORDER BY score DESC"
cur = db.conn.cursor()
cur.execute(sqlSelect)
record = cur.fetchall()
#2
print(record[10])
#db.conn.commit()

# sequences with insertions and deletions 160, 474
seq_index=51
strSeq   = db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]

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
alnDataList = db.appendList_IgorAlignments_data_By_seq_index("D", seq_index)
#db.appendList_IgorAlignments_data_By_seq_index("D", seq_index, alnDataList=alnDataList)
#db.appendList_IgorAlignments_data_By_seq_index("J", seq_index, alnDataList=alnDataList)
alnDataList_pd = create_alnDataList_pandas(seq_index, strSeq, alnDataList)  
alnDataList_pd = addInsertionGaps2alnDataFrame(alnDataList_pd)

#select_alnDataList_pd = pd.DataFrame(alnDataList_pd.loc[0])
#select_alnDataList_pd = pd.concat(select_alnDataList_pd, )
#select_alnDataList_pd.append( alnDataList_pd.loc[ alnDataList_pd['score'] > 16 ] )
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


(alnDataList_pd.loc[1].score)
(alnDataList_pd.loc[2].score)


#ofile = open("example.fasta","w")
#
#ofile.write(">"+str(seq_index)+"\n")
#ofile.write(strAlnSeq+"\n")
#ofile.write(">"+strName_V+"\n")
#ofile.write(strAlnSeq_V+"\n")
#ofile.write(">"+strName_D+"\n")
#ofile.write(strAlnSeq_D+"\n")
#ofile.write(">"+strName_J+"\n")
#ofile.write(strAlnSeq_J+"\n")
#ofile.write(">"+strName_J2+"\n")
#ofile.write(strAlnSeq_J2+"\n")
#
#ofile.close()


# get all alnData.offset in a numpy 
#aln_offsets = list(map(lambda x: x.offset, alnDataList) )
#aln_offsets.append(seq_offset)
#min(aln_offsets)
#aln_offsets
#aaa= np.array(list(map(lambda x: x.offset, alnDataList)).append(seq_offset)).min()
#aaa
#print(len(alnDataList))
#alnDataList[0].strGene_name
#print(v_aligns)
#print(d_aligns)
#print(j_aligns)
#
#print(j_aligns[1])
#
#



#v_alignsSQLrecords = db.fetch_IgorAlignments_By_seq_index("V", seq_index)
#d_alignsSQLrecords = db.fetch_IgorAlignments_By_seq_index("D", seq_index)
#j_alignsSQLrecords = db.fetch_IgorAlignments_By_seq_index("J", seq_index)

#v_offset   = v_alignsSQLrecords[0][3]
#d_offset   = d_alignsSQLrecords[0][3]
#j_offset   = j_alignsSQLrecords[0][3]
#j2_offset   = j_alignsSQLrecords[1][3]
#
#seq_offset = 0
#min_offset = np.array([seq_offset, v_offset, d_offset, j_offset]).min()
#
#strSeq   = db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]
#
#strSeq_V = db.fetch_IgorGeneTemplate_By_gene_id("V", v_alignsSQLrecords[0][1])[2]
#strSeq_D = db.fetch_IgorGeneTemplate_By_gene_id("D", d_alignsSQLrecords[0][1])[2]
#strSeq_J = db.fetch_IgorGeneTemplate_By_gene_id("J", d_alignsSQLrecords[0][1])[2]
#strSeq_J2 = db.fetch_IgorGeneTemplate_By_gene_id("J", d_alignsSQLrecords[1][1])[2]
#
#
#strName_V = db.fetch_IgorGeneTemplate_By_gene_id("V", v_alignsSQLrecords[0][1])[1]
#strName_D = db.fetch_IgorGeneTemplate_By_gene_id("D", d_alignsSQLrecords[0][1])[1]
#strName_J = db.fetch_IgorGeneTemplate_By_gene_id("J", d_alignsSQLrecords[0][1])[1]
#strName_J2 = db.fetch_IgorGeneTemplate_By_gene_id("J", d_alignsSQLrecords[1][1])[1]
#
#strAlnSeq   = (seq_offset-min_offset)*"-"+ strSeq
#strAlnSeq_V = (v_offset-min_offset)*"-"  + strSeq_V
#strAlnSeq_D = (d_offset-min_offset)*"-"  + strSeq_D
#strAlnSeq_J = (j_offset-min_offset)*"-"  + strSeq_J
#strAlnSeq_J2 = (j2_offset-min_offset)*"-"  + strSeq_J2


# So I need a function that recieve the objects to 
# First deal with deletions and insertions
# so the mainn reference is the sequence.
# Given that the main reference is the sequence, insertions are the last part to modify
# on all sequences.
# There must be a somekind of function where  I pass 
# the alignments like



