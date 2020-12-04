#!/usr/bin/python

import pandas as pd
import numpy as np
#import subprocess
import os

#from  optparse import OptionParser
#parser = OptionParser()
#(options, args) = parser.parse_args()

strSpecie = "musmusculus"
strChain  ="TRA"

fln_InputMixcrR1 = "ILN_L_M11_MVA_D6_TCR_R1.t5.fastq.gz" # args[0]
fln_InputMixcrR2 = "ILN_L_M11_MVA_D6_TCR_R2.t5.fastq.gz" # args[1]
basename = fln_InputMixcrR1.split(".")[0][:-3]
print basename
fln_MixcrAlignments = basename+".vdjca" #"ILN_L_M11_MVA_D6_TCR.t5.vdjca"
fln_MixcrAlignmentsTXT = basename+".txt" #"ILN_L_M11_MVA_D6_TCR.t5.txt"



cmd_mixcrAlign = "mixcr align -f -s "+strSpecie+" "+fln_InputMixcrR1+" " +fln_InputMixcrR2+" " +fln_MixcrAlignments 
#-r aligmentReport.txt -OreadsLayout=Collinear \
#-OvParameters.geneFeatureToAling=VTranscript "
cmd_mixcrExport = "mixcr exportAlignments -f -chains --preset full " +fln_MixcrAlignments + " "+ fln_MixcrAlignmentsTXT
# mixcr exportAlignments -f -chains --preset full alignments.vdjca alignments.txt
cmd = cmd_mixcrAlign+" ; " + cmd_mixcrExport
print cmd
os.system(cmd)
##p = subprocess.call(cmd, shell=True)
##p.wait()

df = pd.read_csv(fln_MixcrAlignmentsTXT, delimiter='\t')

dfChain = df.loc[ df['chains'] == strChain, ['targetSequences', 'nSeqCDR3', 'aaSeqCDR3'] ].dropna()
dfChain
booliansChain = []
for cdr3 in dfChain['aaSeqCDR3']:
    if ("*" in cdr3) or ("_" in cdr3) :
        booliansChain.append(True)
    else:
        booliansChain.append(False)

dfChain[booliansChain]
fln_InputSequenceChain  = basename+ "_"+strChain+".csv"
#dfTRA[boolians].loc[:, ['targetSequences', 'nSeqCDR3'] ].to_csv("AA.txt", index=False, sep=';')
#dfTRA[booliansTRA].loc[:, ['targetSequences'] ].to_csv(fln_InputSequenceTRA, index=False, sep=';')
dfChain[booliansChain].loc[:, ['targetSequences'] ].applymap(lambda x : x.split(",")[0]).to_csv(fln_InputSequenceChain, index=False, sep=';')

print strChain

dfChain[booliansChain].loc[:, 'aaSeqCDR3'].nunique()

