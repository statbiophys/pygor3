
##### The event I think is the best one
my_bs = IgorBestScenarios.IgorBestScenariosVDJ() #.load_FromSQLRecord(record_bs[ii])
my_bs.setModel_Parms(flnModelParms)
my_bs.mdl = mdl
my_bs.seq_index = 59
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
realiz_name = 0
realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
my_bs.id_v_3_del = realiz_id


# d_5_del 
pd_event    = mdl.parms.Event_dict['d_5_del']
pd_event['value']
realiz_name = 0
realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
my_bs.id_d_5_del = realiz_id

# d_3_del 
pd_event    = mdl.parms.Event_dict['d_3_del']
pd_event['value']
realiz_name = 0
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
realiz_name = 0
realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
my_bs.id_vd_ins = realiz_id

# dj_ins 
pd_event    = mdl.parms.Event_dict['dj_ins']
pd_event['value']
realiz_name = 0
realiz_id   = pd_event.loc[pd_event['value'] == realiz_name ].index.values[0]
my_bs.id_dj_ins = realiz_id


# FIXME: TMP set to zero ins and dels
my_bs.mismatcheslen = 0

my_bs.mdlParms.Event_dict['v_choice'].loc[my_bs.id_v_choice]['value']
my_bs.getV_Region()





























###### v_choice
strEvent = 'v_choice'
da_event = mdl.xdata[strEvent]
p_V = da_event[{strEvent: bs.id_v_choice}]
p_V
#p_V = da_event.where(da_event['lbl__'+strEvent] == bs.getV_gene_name() ,drop=True)

###### j_choice
strEvent = 'j_choice'
da_event = mdl.xdata[strEvent]
p_J = da_event[{strEvent: bs.id_j_choice}]
p_J
#p_J = da_event.where(da_event['lbl__'+strEvent] == bs.getJ_gene_name() ,drop=True)

###### d_gene
strEvent = 'd_gene'
da_event = mdl.xdata[strEvent]
p_DgJ = da_event[{'d_gene':bs.id_d_gene, 'j_choice':bs.id_j_choice}]
p_DgJ


###### v_3_del
strEvent = 'v_3_del'
da_event = mdl.xdata[strEvent]
p_V_3_del = da_event[{'v_3_del':bs.id_v_3_del, 'v_choice':bs.id_v_choice}]
p_V_3_del


###### j_5_del
strEvent = 'j_5_del'
da_event = mdl.xdata[strEvent]
p_J_5_del = da_event[{'j_5_del':bs.id_j_5_del, 'j_choice':bs.id_j_choice}]
p_J_5_del

###### d_5_del
strEvent = 'd_5_del'
da_event = mdl.xdata[strEvent]
p_D_5_del = da_event[{'d_5_del':bs.id_d_5_del, 'd_gene':bs.id_d_gene}]
p_D_5_del

###### d_3_del
strEvent = 'd_3_del'
da_event = mdl.xdata[strEvent]
p_D_3_del = da_event[{'d_3_del':bs.id_d_3_del, 'd_5_del':bs.id_d_5_del, 'd_gene':bs.id_d_gene}]
p_D_3_del

###### vd_ins
strEvent = 'vd_ins'
da_event = mdl.xdata[strEvent]
p_VD_ins = da_event[{'vd_ins':bs.id_vd_ins}]
p_VD_ins


###### vd_dinucl
strEvent = 'vd_dinucl'
da_event = mdl.xdata[strEvent]
# Get the last nucleotide of V region (after deletions)

str_prev_nt = bs.getV_Region()[-1]
pd_tmp = mdl.parms.Event_dict[strEvent]
prev_nt = pd_tmp.loc[pd_tmp['value'] == str_prev_nt].index.values[0]

# for each nucleotide on inserted list
Num_nt = 4 # 4 nucleotides A, C, G, T 
p_VD_dinucl = 1
for curr_nt in bs.vd_dinucl:
    id_dinucl   = prev_nt*Num_nt + curr_nt
    prob_tmp    = da_event[{'vd_dinucl':id_dinucl}]
    p_VD_dinucl = p_VD_dinucl * prob_tmp
    print(prev_nt, curr_nt, id_dinucl, prob_tmp, p_VD_dinucl)
    prev_nt = curr_nt

p_VD_dinucl

###### dj_ins
strEvent = 'dj_ins'
da_event = mdl.xdata[strEvent]
p_DJ_ins = da_event[{'dj_ins':bs.id_dj_ins}]
p_DJ_ins

###### dj_dinucl
strEvent = 'dj_dinucl'
da_event = mdl.xdata[strEvent]
# Get the last nucleotide of V region (after deletions)

str_prev_nt = (bs.getV_Region() + bs.getVD_Region() + bs.getD_Region() )[-1]
pd_tmp = mdl.parms.Event_dict[strEvent]
prev_nt = pd_tmp.loc[pd_tmp['value'] == str_prev_nt].index.values[0]

# for each nucleotide on inserted list
Num_nt = 4 # 4 nucleotides A, C, G, T 
p_DJ_dinucl = 1
for curr_nt in bs.vd_dinucl:
    id_dinucl   = prev_nt*Num_nt + curr_nt
    prob_tmp    = da_event[{'dj_dinucl':id_dinucl}]
    p_DJ_dinucl = p_DJ_dinucl * prob_tmp
    print(prev_nt, curr_nt, id_dinucl, prob_tmp, p_DJ_dinucl)
    prev_nt = curr_nt

p_DJ_dinucl



mdl.parms.G.edges()
mdl.parms.G.nodes()

perror = float(mdl.parms.ErrorRate['SingleErrorRate'])

bs.scenario_proba_cond_seq
bs.mismatcheslen

p_vecE = p_V*p_J*p_DgJ*p_V_3_del*p_J_5_del*p_D_5_del*p_D_3_del*p_VD_ins*p_VD_dinucl*p_DJ_ins*p_DJ_dinucl*(perror*bs.mismatcheslen)
p_vecE.values

p_cond =   p_V*(1-p_V)*\
            p_J*(1-p_J)*\
            p_DgJ*(1-p_DgJ)*\
            p_V_3_del*(1-p_V_3_del)*\
            p_J_5_del*(1-p_J_5_del)*\
            p_D_5_del*(1-p_D_5_del)*\
            p_D_3_del*(1-p_D_3_del)*\
            p_VD_ins*(1-p_VD_ins)*\
            p_VD_dinucl*(1-p_VD_dinucl)*\
            p_DJ_ins*(1-p_DJ_ins)*\
            p_DJ_dinucl*(1-p_DJ_dinucl)
            #*(perror*bs.mismatcheslen)

p_cond.values



p_cond02 =   (1-p_V)*\
            (1-p_J)*\
            (1-p_DgJ)*\
            (1-p_V_3_del)*\
            (1-p_J_5_del)*\
            (1-p_D_5_del)*\
            (1-p_D_3_del)*\
            (1-p_VD_ins)*\
            (1-p_VD_dinucl)*\
            (1-p_DJ_ins)*\
            (1-p_DJ_dinucl)

p_cond02.values

pgen = 5.07785e-16

pgen/p_vecE.values
bs.scenario_proba_cond_seq








bs0 = IgorBestScenarios.IgorBestScenariosVDJ.load_FromSQLRecord(record_bs[0])
bs1 = IgorBestScenarios.IgorBestScenariosVDJ.load_FromSQLRecord(record_bs[1])

bs0.strSeq   = db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]
bs1.strSeq   = db.fetch_IgorIndexedSeq_By_seq_index(seq_index)[1]

    

























aaa = IgorAlignment_data.IgorAlignment_data.load_FromSQLRecord(j_aligns[1])
aaa.seq_index
aaa.strGene_name


print(aaa.deletions)
print(aaa.mismatches)

#def writeIgorAlignments2Fasta(strSeq, listOfAlignmentsObjects ):


# TODO: make a print of the sequence like an aln file or fasta
ofile = open("example.fasta","w")

ofile.write(">"+str(seq_index)+"\n")
ofile.write(strAlnSeq+"\n")
ofile.write(">"+strName_V+"\n")
ofile.write(strAlnSeq_V+"\n")
ofile.write(">"+strName_D+"\n")
ofile.write(strAlnSeq_D+"\n")
ofile.write(">"+strName_J+"\n")
ofile.write(strAlnSeq_J+"\n")
ofile.write(">"+strName_J2+"\n")
ofile.write(strAlnSeq_J2+"\n")

ofile.close()





#
#alnDataList_pd['deletionsAUX'] = alnDataList_pd['deletions'].apply(np.array)
#alnDataList_pd['deletionsAUX']
#alnDataList_pd['seq_with_ins'].loc[0]


#addIns2Seq(alnDataList_pd['seq_no_align'], alnDataList_pd['insertions'])

alnDataList_pd.keys()
alnDataList_pd['deletions']
alnDataList_pd['offset']
alnDataList_pd['relative_offset']

# TODO: 
# 1. For seq_no_align not align I have the deletions
# 2. but I need to positions of the deletions in seq_with_ins_offset
# 3. Make a new list of deletions using seq_with_ins.
#  3.1. How? by checking how many gaps have been inserted at the position
#  3.2. Example

# original_deletions [1, 4, 8, 10]
# tgatatccctgat
# tga-ta-tccc-tgat
deletionsOrig = [9, 12, 16, 18]
#strSeq     =        "tgaAtaGtcccGtgat"
strSeq          =        "taAtGtccGtataggg"
strSeqSemiFinal =        "t-aAt-Gtcc-Gt-ataggg"
offset          =  -8
min_offset      = -8 #-12
strOriginal     ="ggatggcctgatatccctgat"
strWithGaps     ="ggatggcctga-ta-tccc-tgat"

org_delPos = deletionsOrig[3]
print((0-min_offset)*"-"+strSeq)
print((0-min_offset)*"-"+strSeqSemiFinal)

print(strOriginal[org_delPos])
print(strOriginal[:org_delPos])
new_delPos = org_delPos + strWithGaps[:org_delPos].count('-')
print(strWithGaps[:new_delPos])
print(strWithGaps[new_delPos])

import re
strWithGaps
strGap='-'
gapsPos = np.array([m.start() for m in re.finditer(strGap, strWithGaps)])
gapsPos
deletionsOrigArr = np.array(deletionsOrig)
gapsPos
deletionsAUX = deletionsOrigArr
print("original: ", deletionsAUX)
for gap in gapsPos:
    ppp = deletionsAUX[deletionsAUX <= gap ]
    qqq = deletionsAUX[deletionsAUX >  gap ]+1
    deletionsAUX = np.concatenate((ppp, qqq))
    print(deletionsAUX)

deletionsOrig[0] in gapsPos


x = np.array([0.2, 6.4, 3.0, 1.6])
bins = np.array([0.0, 1.0, 2.5, 4.0, 10.0])
inds = np.digitize(deletionsOrig[0], gapsPos)
inds


strWithGaps[19]

print(len(strOriginal[:delPos]))
deletionsWithIns = list()
for dels in deletionsOrig:
    deletionsWithIns.append(dels)
    
    print(strWithGaps[dels])
    print(strOriginal[dels])
    

alnDataList_pd[['seq_with_ins_offset', 'offset']]

# Write a Fasta file for alignments



alnDataList_pd['relative_offset']
alnDataList_pd['deletions'].map(len)

alnDataList_pd['seq_no_align']
alnDataList_pd.keys()
print(strSeq_Orig)
joder = addIns2Seq(strSeq_Orig, insList)
print(joder)
if iii == sorted(iii):
    for i in iii:
        print(i)






#
#import IgorAlignment_data
#seq_index = 51
#dataList = list()
#strGene = "V"
#alignsSQLrecords = db.fetch_IgorAlignments_By_seq_index(strGene, seq_index)
#for alignSQLrecord in alignsSQLrecords:
#    #print(align)
#    gene_id = alignSQLrecord[1]
#    strGene_name = db.fetch_IgorGeneTemplate_By_gene_id(strGene, gene_id)[1]
#    aln_data = IgorAlignment_data.IgorAlignment_data.load_FromSQLRecord(alignSQLrecord, strGene_name=strGene_name)
#    aln_data.strGene_class = strGene
#    dataList.append(aln_data)
#    #if not len(aln_data.insertions) == 0:
#    print(aln_data.strGene_name, aln_data.score, aln_data.offset, aln_data.insertions)
#
#dataList[0].strGene_class
#joder = list(map(lambda x: x.strGene_name, dataList))
#joder
#joder[0]
#
#print(dataList[0].strGene_name)


#
#
#print(db.fetch_IgorGeneTemplate_By_gene_id("V", v_aligns[0][1])[2] )
#print((d_offset-min_offset)*" "+db.fetch_IgorGeneTemplate_By_gene_id("D", d_aligns[0][1])[2] )
#print((j_offset-min_offset)*" "+db.fetch_IgorGeneTemplate_By_gene_id("J", j_aligns[0][1])[2] )
#
#
#
#print("V", db.fetch_IgorGeneTemplate_By_gene_id("V", v_aligns[0][1])[2] )
#print("D", db.fetch_IgorGeneTemplate_By_gene_id("D", d_aligns[0][1])[2] )
#print("J", db.fetch_IgorGeneTemplate_By_gene_id("J", j_aligns[0][1])[2] )

# Now we need a pandas to manipulate the alignment data.


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
