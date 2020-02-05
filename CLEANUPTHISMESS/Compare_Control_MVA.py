#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  3 09:30:58 2019

@author: alfaceor
"""
import IgorModel

IgorSpecie    = "mouse"
IgorChain     = "tcr_beta"

IgorModelPath = "../../IGoR/models/"+IgorSpecie+"/"+IgorChain+"/"
IgorRefGenomePath = IgorModelPath+"ref_genome/"

flnVGeneTemplate = IgorRefGenomePath+"genomicVs.fasta"
flnDGeneTemplate = IgorRefGenomePath+"genomicDs.fasta"
flnJGeneTemplate = IgorRefGenomePath+"genomicJs.fasta"


#### Control data
#strWD="ReDo_Control/"
strWD="ReDo_Ctrl/"
batchname = "TRBCtrl"
lblOld = "Ctrl"
IgorModelInferencePath = strWD + batchname + "_inference/"
flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"

mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#### MVA_D6 data
strWD="ReDo_MVA_D6/"
batchname = "TRBNew"
lblNew = "MVA D6"
IgorModelInferencePath = strWD + batchname + "_inference/"
flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"

mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

flnFigPrefix="Ctrl__MVA_D6__"


### load IGoR model parms and marginals.
import matplotlib.pyplot as plt
###### BEGIN PLOT DECORATION VARIABLES
#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 20} 
#plt.rc('font', **font)
#plt.rc('text', usetex=True)
###### END PLOT DECORATION VARIABLES


########### GENE CHOICE ###########
##### v_choice #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'v_choice'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"


daOld = mdlOld.xdata[strEvent]
daNew = mdlNew.xdata[strEvent]

daOld.plot(ax=ax, marker = 'o', label=lblOld)
daNew.plot(ax=ax, marker = 'o', label=lblNew)

ax.legend(bbox_to_anchor=(0.85, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values, rotation=90)
fig.tight_layout()

fig.savefig(flnFigPrefix+strEvent+".pdf")


##### j_choice #####\
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'j_choice'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"


daOld = mdlOld.xdata[strEvent]
daNew = mdlNew.xdata[strEvent]

daOld.plot(ax=ax, marker = 'o', label=lblOld)
daNew.plot(ax=ax, marker = 'o', label=lblNew)

ax.legend(bbox_to_anchor=(0.85, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values, rotation=90)
fig.tight_layout()

fig.savefig(flnFigPrefix+strEvent+".pdf")



##### d_gene #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'd_gene'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"

daNew = mdlNew.xdata[strEvent]*mdlNew.xdata['j_choice']
daNew = daNew.sum(dim='j_choice')

daOld = mdlOld.xdata[strEvent]*mdlOld.xdata['j_choice']
daOld = daOld.sum(dim='j_choice')

daOld.plot(ax=ax, marker = 'o', label=lblOld)
daNew.plot(ax=ax, marker = 'o', label=lblNew)

ax.legend(bbox_to_anchor=(0.85, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values, rotation=90)
fig.tight_layout()

fig.savefig(flnFigPrefix+strEvent+".pdf")


########### DELETIONS ###########
##### v_3_del #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'v_3_del'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"


daNew = mdlNew.xdata[strEvent]*mdlNew.xdata['v_choice']
daNew = daNew.sum(dim='v_choice')

daOld = mdlOld.xdata[strEvent]*mdlOld.xdata['v_choice']
daOld = daOld.sum(dim='v_choice')

daOld.plot(ax=ax, marker = 'o', label=lblOld)
daNew.plot(ax=ax, marker = 'o', label=lblNew)

ax.legend(bbox_to_anchor=(0.85, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values)
fig.tight_layout()

fig.savefig(flnFigPrefix+strEvent+".pdf")


##### j_5_del #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'j_5_del'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"

daNew = mdlNew.xdata[strEvent]*mdlNew.xdata['j_choice']
daNew = daNew.sum(dim='j_choice')

daOld = mdlOld.xdata[strEvent]*mdlOld.xdata['j_choice']
daOld = daOld.sum(dim='j_choice')

daOld.plot(ax=ax, marker = 'o', label=lblOld)
daNew.plot(ax=ax, marker = 'o', label=lblNew)

ax.legend(bbox_to_anchor=(0.85, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values)
fig.tight_layout()

fig.savefig(flnFigPrefix+strEvent+".pdf")



##### d_5_del #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'd_5_del'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"

daNew = mdlNew.xdata['d_gene']*mdlNew.xdata['j_choice']
daNew = daNew.sum(dim='j_choice')
daNew = mdlNew.xdata[strEvent]*daNew
daNew = daNew.sum(dim='d_gene')

daOld = mdlOld.xdata['d_gene']*mdlOld.xdata['j_choice']
daOld = daOld.sum(dim='j_choice')
daOld = mdlOld.xdata[strEvent]*daOld
daOld = daOld.sum(dim='d_gene')

daOld.plot(ax=ax, marker = 'o', label=lblOld)
daNew.plot(ax=ax, marker = 'o', label=lblNew)

ax.legend(bbox_to_anchor=(0.85, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values)
fig.tight_layout()

fig.savefig(flnFigPrefix+strEvent+".pdf")


##### d_3_del #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'd_3_del'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"

daNew = mdlNew.xdata['d_gene']*mdlNew.xdata['j_choice']
daNew = daNew.sum(dim='j_choice') # P(D)
daNew = mdlNew.xdata['d_5_del']*daNew # P(delD5|D)*P(D)
daNew = mdlNew.xdata[strEvent]*daNew # P(delD3|delD5,D)*P(delD5|D)*P(D)
daNew = daNew.sum(dim=('d_5_del', 'd_gene'))


daOld = mdlOld.xdata['d_gene']*mdlOld.xdata['j_choice']
daOld = daOld.sum(dim='j_choice') # P(D)
daOld = mdlOld.xdata['d_5_del']*daOld # P(delD5|D)*P(D)
daOld = mdlOld.xdata[strEvent]*daOld # P(delD3|delD5,D)*P(delD5|D)*P(D)
daOld = daOld.sum(dim=('d_5_del', 'd_gene'))

daOld.plot(ax=ax, marker = 'o', label=lblOld)
daNew.plot(ax=ax, marker = 'o', label=lblNew)

ax.legend(bbox_to_anchor=(0.85, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values)
fig.tight_layout()

fig.savefig(flnFigPrefix+strEvent+".pdf")



########### INSERTIONS ###########
##### vd_ins #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'vd_ins'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"


daNew = mdlNew.xdata[strEvent]
daOld = mdlOld.xdata[strEvent]

daOld.plot(ax=ax, marker = 'o', label=lblOld)
daNew.plot(ax=ax, marker = 'o', label=lblNew)

ax.legend(bbox_to_anchor=(0.85, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values)
fig.tight_layout()

fig.savefig(flnFigPrefix+strEvent+".pdf")



##### dj_ins #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'dj_ins'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"


daNew = mdlNew.xdata[strEvent]
daOld = mdlOld.xdata[strEvent]

daOld.plot(ax=ax, marker = 'o', label=lblOld)
daNew.plot(ax=ax, marker = 'o', label=lblNew)

ax.legend(bbox_to_anchor=(0.85, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values)
fig.tight_layout()

fig.savefig(flnFigPrefix+strEvent+".pdf")



"""
fig, ax = plt.subplots(figsize=(12,6))

strEvent  = 'd_gene'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "initial_model.txt"
flnModelMargs = IgorModelInferencePath + "initial_marginals.txt"
mdl = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

da = mdl.xdata[strEvent]*mdl.xdata['j_choice']
da = da.sum(dim='j_choice')
da.plot(ax=ax, marker = 'o', label="Initial")


for i in [1, 2, 3, 4, 5]:
    flnModelParms = IgorModelInferencePath + "iteration_"+str(i)+"_parms.txt"
    flnModelMargs = IgorModelInferencePath + "iteration_"+str(i)+".txt"
    
    mdl = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
    #da = mdl.xdata[strEvent]
    da = mdl.xdata[strEvent]*mdl.xdata['j_choice']
    da = da.sum(dim='j_choice')
    da.plot(ax=ax, marker = 'o', label=str(i))

ax.legend(bbox_to_anchor=(0.80, 0.52), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdl.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdl.xdata[strEvent]['lbl__'+strEvent].values, rotation=90)
fig.tight_layout()

#fig.savefig(strEvent+".pdf")
"""






