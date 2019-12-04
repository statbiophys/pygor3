#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 05:13:57 2019

@author: alfaceor
"""

import IgorBestScenarios
import IgorModel

strWD="uchicago/selected_sequences/"
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
mdl = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
mdlParms = IgorModel.Model_Parms(flnModelParms) # mdl.parms
mdlMargs = IgorModel.Model_Marginals(flnModelMargs) # mdl.marginals

# 1. Plot model graph

# 2. Plot marginals
#mdl.xdata['v_choice'].plot()
mdl.xdata['dj_ins'].plot()

mdl.xdata['dj_dinucl'].plot()

### IGoR ouptut files
flnIgorBestScenarios = strWD+batchname+"_output/best_scenarios_counts.csv"


# Good9
scenario_good9 ={
"vd_ins" : 12,
"dj_ins" : 8,
"j_5_del" : 11,
"vd_dinucl" : "[2,1,0,2,1,3,1,1,3,3,2,0]",
"d_gene" : 0,
"dj_dinucl" : "[1,0,1,0,0,0,1,2]",
"v_choice" : 3,
"j_choice" : 0,
"v_3_del" : 14,
"d_5_del" : 1,
"d_3_del" : 19
}


# Bad
scenario = {
"vd_ins" : 12,
"dj_ins" : 11,
"j_5_del" : 8,
"vd_dinucl" : "[1,1,0,2,1,2,2,3,2,0,3,2]",
"d_gene" : 1,
"dj_dinucl" : "[1,2,3,2,0,0,2,2,2,2,2]",
"v_choice" : 14,
"j_choice" : 9,
"v_3_del" : 16,
"d_5_del" : 8,
"d_3_del" : 12
}

######
scenario_igor_bad = {
"vd_ins" : 12,
"dj_ins" : 11,
"j_5_del" : 8,
"vd_dinucl" : "[1,1,0,2,1,2,2,3,2,0,3,2]",
"d_gene" : 1,
"dj_dinucl" : "[1,2,3,2,0,0,2,2,2,2,2]",
"v_choice" : 14,
"j_choice" : 9,
"v_3_del" : 16,
"d_5_del" : 8,
"d_3_del" : 12
}

line_scenario_igor_bad = "1;1;0.123596;(14);(9);(1);(4);(8);(7);(11);(0);();(9);(0,2,0,1,2,3,2,0,0);(122,123,124)"


bs_igor_good9 = IgorBestScenarios.IgorBestScenariosVDJ.load_FromDict(scenario_good9)
bs_igor_good9.setModel_Parms(flnModelParms)
bs_igor_good9.mdl = mdl
bs_igor_good9.get_EventProb()
print(bs_igor_good9)


bs_igor_bad = IgorBestScenarios.IgorBestScenariosVDJ.load_FromLineBestScenario(line_scenario_igor_bad) #.load_FromSQLRecord(record_bs[ii])
bs_igor_bad.setModel_Parms(flnModelParms)
bs_igor_bad.mdl = mdl
bs_igor_bad.get_EventProb()
print(bs_igor_bad)


bs_igor_bad2 = IgorBestScenarios.IgorBestScenariosVDJ.load_FromDict(scenario_igor_bad)#.load_FromSQLRecord(record_bs[ii])
bs_igor_bad2.setModel_Parms(flnModelParms)
bs_igor_bad2.mdl = mdl
bs_igor_bad2.get_EventProb()
print(bs_igor_bad)


dict_scenario_mixcr_bad= {
"v_choice" : 14,
"j_choice" : 9,
"d_gene" : 1,
"v_3_del" : 4,
"j_5_del" : 4,
"d_5_del" : 9,
"d_3_del" : 3,
"vd_ins" : 1,
"vd_dinucl" : "[]",
"dj_ins" : 1,
"dj_dinucl" : "[0]"
}


bs_mixcr_bad = IgorBestScenarios.IgorBestScenariosVDJ.load_FromDict(dict_scenario_mixcr_bad) #.load_FromSQLRecord(record_bs[ii])
bs_mixcr_bad.setModel_Parms(flnModelParms)
bs_mixcr_bad.mdl = mdl
bs_mixcr_bad.get_EventProb()

bs_mixcr_bad.






# calculate prob of Seq 9

#nickname: dj_dinucl;proba_contribution : 7.2938e-06;new_scenario_proba : 1.21359e-26
#==== Rec_Event::iterate_wrap_up
#scenario_proba : 1.21359e-26


print(scenario)

my_bs = IgorBestScenarios.IgorBestScenariosVDJ.load_FromDict(scenario) #.load_FromSQLRecord(record_bs[ii])

my_bs.setModel_Parms(flnModelParms)
my_bs.mdl = mdl

my_bs.getD_Region()
my_bs.getVD_Region()
P_scenario = my_bs.get_EventProb()

(my_bs.getV_Region() + my_bs.getVD_Region() + my_bs.getD_Region() )[-1]
print(P_scenario)

my_bs.get_ErrorProb()
Probas = my_bs.get_DictNicknameProbs()
print(Probas  )
Probas['v_choice'] * \
Probas['j_choice'] * \
Probas['d_gene'] * \
Probas['v_3_del'] * \
Probas['d_5_del'] * \
Probas['d_3_del'] * \
Probas['j_5_del'] * \
Probas['vd_ins'] * \
Probas['vd_dinucl'] * \
Probas['dj_ins'] * \
Probas['dj_dinucl']


my_bs.getV_Region() + my_bs.getVD_Region() + my_bs.getD_Region()


#print(my_bs)
