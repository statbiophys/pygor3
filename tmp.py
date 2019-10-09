

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
