#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  3 16:03:42 2019

@author: alfaceor
"""

import IgorSqliteDB

flnIgorDB = "joder.db"
flnIgorIndexedSeq = "pytest/aligns/Barb02_indexed_sequences.csv"

flnVGeneTemplate = "../../IGoR/models/human/tcr_beta/ref_genome/genomicVs.fasta"
flnDGeneTemplate = "../../IGoR/models/human/tcr_beta/ref_genome/genomicDs.fasta"
flnJGeneTemplate = "../../IGoR/models/human/tcr_beta/ref_genome/genomicJs.fasta"

flnVAlignments = "pytest/aligns/Barb02_V_alignments.csv"
flnDAlignments = "pytest/aligns/Barb02_D_alignments.csv"
flnJAlignments = "pytest/aligns/Barb02_J_alignments.csv"

db = IgorSqliteDB.IgorSqliteDB()
db.createSqliteDB(flnIgorDB)
db.load_VDJ_Database(flnIgorIndexedSeq, flnVGeneTemplate, flnDGeneTemplate, flnJGeneTemplate, flnVAlignments, flnDAlignments, flnJAlignments)

seq_index=929
aaa = db.fetch_IgorIndexedSeq_By_seq_index(seq_index)

v_align = db.fetch_IgorAlignments_By_seq_index("V", seq_index)
d_align = db.fetch_IgorAlignments_By_seq_index("D", seq_index)
j_align = db.fetch_IgorAlignments_By_seq_index("J", seq_index)

print(v_align[1])
print(d_align[1])
print(j_align[1])

print(db.fetch_IgorGeneTemplate_By_gene_id("V", v_align[1]) )

#query = "SELECT vgene_id, sequence FROM IgorVGeneTemplate WHERE seq_index==929 ORDER BY score DESC"

#query = '''SELECT seq_index, IgorVGeneTemplate.sequence 
#        FROM IgorVAlignments 
#        INNER JOIN IgorVGeneTemplate ON 
#        IgorVAlignments.vgene_id = IgorVGeneTemplate.vgene_id
#        '''


#query = '''SELECT seq_index, IgorVGeneTemplate.sequence 
#        FROM IgorVAlignments 
#        INNER JOIN IgorVGeneTemplate ON 
#        IgorVAlignments.vgene_id = IgorVGeneTemplate.vgene_id
#        '''

#cur.execute(query)
#results = cur.fetchall()
#print(results)
#import pandas as pd
#colnames=["seq_index", "gene_id", "score", "offset", "insertions", "deletions", "mismatches", "length", "offset_5_p_align", "offset_3_p_align" ] 
#print(len(colnames))
#df = pd.DataFrame(results, columns=colnames)
#df['offset'].loc[0]
#df.loc[0]

