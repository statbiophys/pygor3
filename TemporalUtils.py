# TODO: Create a pandas Dataframe with the following columns: 
# description, no_align_seq, offset, align_seq, insertions, deletions
# TODO: FIRST put the mismatches with CAPITAL LETTERS
import pandas as pd
import numpy as np
import re





def create_alnDataList_pandas(seq_index, strSeq, alnDataList):
    seq_offset = 0
    alnDataList_pd = pd.DataFrame(columns=['description', 'seq_no_align', 'offset'])
    alnDataList_pd = alnDataList_pd.append(\
                   {'description' : str(seq_index), \
                    'seq_no_align': strSeq.lower(), \
                    'offset'      : seq_offset ,\
                    'insertions'  : list(), \
                    'deletions'   : list(), \
                    'mismatches'  : list() ,\
                    'score'       : -1
                    }, ignore_index=True)

    # insert alignments in a pandas dataframe.
    for alnData in alnDataList:
        alnDataList_pd = alnDataList_pd.append(\
                       {'description'  : alnData.strGene_name, \
                        'seq_no_align' : alnData.strGene_seq.lower(), \
                        'offset'       : alnData.offset, \
                        'insertions'   : alnData.insertions, \
                        'deletions'    : alnData.deletions, \
                        'mismatches'   : alnData.mismatches, \
                        'score'   : alnData.score \
                        }, ignore_index=True)
    
    return alnDataList_pd


def addGap2Seq(strSeq_Orig, insPosList, offset):
    if len(insPosList) == 0:
        return strSeq_Orig
    else:
        strGap = "-"
        lastInsPos = 0
        strSeq_withIns = ""
        #print("antes", insPosList)
        insPosList = [i - offset for i in insPosList]
        #print("despues",insPosList)
        for insPos in insPosList:
            strSeq_withIns = strSeq_withIns+strSeq_Orig[lastInsPos:insPos]+strGap
            lastInsPos = insPos
        strSeq_withIns = strSeq_withIns+strSeq_Orig[lastInsPos:]
        #print(strSeq_withIns)
        return strSeq_withIns

def addInsertionGaps2alnDataFrame(alnDataList_pd):
    min_offset = alnDataList_pd['offset'].min()
    alnDataList_pd['relative_offset']  = alnDataList_pd['offset'] - min_offset
    alnDataList_pd['seq_align_offset'] = alnDataList_pd['relative_offset']*"-"+ alnDataList_pd['seq_no_align']

    alnDataList_pd['seq_with_ins'] = alnDataList_pd.apply(lambda x: addGap2Seq(x['seq_no_align'], x['insertions'], x['offset']) , axis=1)
    alnDataList_pd['seq_with_ins_offset'] = alnDataList_pd['relative_offset']*"-"+ alnDataList_pd['seq_with_ins']
    
    alnDataList_pd['seq_with_ins'] = alnDataList_pd.apply(lambda x: addGap2Seq(x['seq_no_align'], x['insertions'], x['offset']) , axis=1)

    return alnDataList_pd


#def addDeletionsGaps2alnDataFrame(alnDataList_pd):
#    min_offset = alnDataList_pd['offset'].min()
#    alnDataList_pd['deletions']

def modifyDeletions(strWithGaps, deletionsOrig, relative_offset):
    strGap='-'
    gapsPos = np.array([m.start() for m in re.finditer(strGap, strWithGaps)])
    deletionsOrigArr = np.array(deletionsOrig)
    deletionsAUX = deletionsOrigArr
    print("original: ", deletionsAUX)
    for gap in gapsPos:
        ppp = deletionsAUX[deletionsAUX <= gap ]
        qqq = deletionsAUX[deletionsAUX >  gap ]+1
        deletionsAUX = np.concatenate((ppp, qqq))
        #print(deletionsAUX)
    deletionsAUX = deletionsAUX + relative_offset
    return list(deletionsAUX)


# FIXME: FIND A BETTER SOLUTIONS TO ADD THE DELETIONS ON THE INPUT SEQUENCE
# A lot of work to do here.
def addTemporalDelsInSeq(alnDataList_pd, id_pd):
    #alnDataList_pd['seq_with_ins'].loc[0]
    newDelPosList = modifyDeletions( alnDataList_pd['seq_with_ins'].loc[id_pd], alnDataList_pd['deletions'].loc[id_pd], alnDataList_pd['relative_offset'].loc[id_pd])
    #modifyDeletions( alnDataList_pd['seq_with_ins'].loc[0], alnDataList_pd['deletions'].loc[0], alnDataList_pd['relative_offset'].loc[0])
    # FIXME: this is wrong! correct the offset
    strSeqWithDels =  addGap2Seq(alnDataList_pd['seq_with_ins_offset'].loc[0], newDelPosList, 0)
    #alnDataList_pd['seq_final'].loc[id_pd] =  addGap2Seq(alnDataList_pd['seq_with_ins_offset'].loc[id_pd], newDelPosList, 0)
    return strSeqWithDels


def writeAlignmentsFastaOnlyInsertions(alnDataList_pd, flnFastaOutput):
    ofile = open(flnFastaOutput,"w")
    for index, row in alnDataList_pd.iterrows():
        #print(row['description'], row['seq_align_offset'])
        ofile.write(">"+row['description']+", score = "+str(row['score'])+"\n")
        ofile.write(row['seq_with_ins_offset']+"\n")
    ofile.close()



