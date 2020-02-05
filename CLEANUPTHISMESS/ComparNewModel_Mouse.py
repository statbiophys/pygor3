#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 16:22:13 2019

@author: alfaceor
"""

batchname = "TRbeta"
strWD="ReDoTRbeta/"
IgorSpecie    = "mouse"
IgorChain     = "tcr_beta"

IgorModelPath = "../../IGoR/models/"+IgorSpecie+"/"+IgorChain+"/"
IgorRefGenomePath = IgorModelPath+"ref_genome/"

flnVGeneTemplate = IgorRefGenomePath+"genomicVs.fasta"
flnDGeneTemplate = IgorRefGenomePath+"genomicDs.fasta"
flnJGeneTemplate = IgorRefGenomePath+"genomicJs.fasta"

#flnModelParms = IgorModelPath + "models/model_parms.txt"
#flnModelMargs = IgorModelPath + "models/model_marginals.txt"


### load IGoR model parms and marginals.
import IgorModel
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
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"
mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#da = mdl.xdata[strEvent]
daNew = mdlNew.xdata[strEvent] #*mdlNew.xdata['j_choice']
daOld = mdlOld.xdata[strEvent] #*mdlOld.xdata['j_choice']
#daNew = daNew.sum(dim='j_choice')
daOld.plot(ax=ax, marker = 'o', label="Current")
daNew.plot(ax=ax, marker = 'o', label="Inferred")

ax.legend(bbox_to_anchor=(0.85, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values, rotation=90)
fig.tight_layout()

fig.savefig(strEvent+"_compar.pdf")

##### j_choice #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'j_choice'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"
mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#da = mdl.xdata[strEvent]
daNew = mdlNew.xdata[strEvent] #*mdlNew.xdata['j_choice']
daOld = mdlOld.xdata[strEvent] #*mdlOld.xdata['j_choice']
daOld.plot(ax=ax, marker = 'o', label="Current")
daNew.plot(ax=ax, marker = 'o', label="Inferred")

ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values, rotation=90)
fig.tight_layout()

fig.savefig(strEvent+"_compar.pdf")


##### d_gene #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'd_gene'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"
mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#da = mdl.xdata[strEvent]
daNew = mdlNew.xdata[strEvent]*mdlNew.xdata['j_choice']
daNew = daNew.sum(dim='j_choice')

daOld = mdlOld.xdata[strEvent]*mdlOld.xdata['j_choice']
daOld = daOld.sum(dim='j_choice')

daOld.plot(ax=ax, marker = 'o', label="Current")
daNew.plot(ax=ax, marker = 'o', label="Inferred")

ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values, rotation=90)
fig.tight_layout()

fig.savefig(strEvent+"_compar.pdf")

########### DELETIONS ###########
##### v_3_del #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'v_3_del'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"
mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#da = mdl.xdata[strEvent]
daNew = mdlNew.xdata[strEvent]*mdlNew.xdata['v_choice']
daNew = daNew.sum(dim='v_choice')

daOld = mdlOld.xdata[strEvent]*mdlOld.xdata['v_choice']
daOld = daOld.sum(dim='v_choice')

daOld.plot(ax=ax, marker = 'o', label="Current")
daNew.plot(ax=ax, marker = 'o', label="Inferred")

ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values) #, rotation=90)
fig.tight_layout()

fig.savefig(strEvent+"_compar.pdf")

##### j_5_del #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'j_5_del'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"
mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#da = mdl.xdata[strEvent]
daNew = mdlNew.xdata[strEvent]*mdlNew.xdata['j_choice']
daNew = daNew.sum(dim='j_choice')

daOld = mdlOld.xdata[strEvent]*mdlOld.xdata['j_choice']
daOld = daOld.sum(dim='j_choice')

daOld.plot(ax=ax, marker = 'o', label="Current")
daNew.plot(ax=ax, marker = 'o', label="Inferred")

ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values) #, rotation=90)
fig.tight_layout()

fig.savefig(strEvent+"_compar.pdf")

##### d_5_del #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'd_5_del'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"
mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#da = mdl.xdata[strEvent]
#daNew = mdlNew.xdata[strEvent]*mdlNew.xdata['d_gene']*mdlNew.xdata['j_choice']
#daNew = daNew.sum(dim='d_gene')
daNew = mdlNew.xdata['d_gene']*mdlNew.xdata['j_choice']
daNew = daNew.sum(dim='j_choice')
daNew = mdlNew.xdata[strEvent]*daNew
daNew = daNew.sum(dim='d_gene')

#daOld = mdlOld.xdata[strEvent]*mdlOld.xdata['d_gene']
#daOld = daOld.sum(dim='d_gene')
daOld = mdlOld.xdata['d_gene']*mdlOld.xdata['j_choice']
daOld = daOld.sum(dim='j_choice')
daOld = mdlOld.xdata[strEvent]*daOld
daOld = daOld.sum(dim='d_gene')


daOld.plot(ax=ax, marker = 'o', label="Current")
daNew.plot(ax=ax, marker = 'o', label="Inferred")

ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values) #, rotation=90)
fig.tight_layout()

fig.savefig(strEvent+"_compar.pdf")

##### d_3_del #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'd_3_del'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"
mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#da = mdl.xdata[strEvent]
#daNew = mdlNew.xdata[strEvent]*mdlNew.xdata['d_gene']*mdlNew.xdata['j_choice']
#daNew = daNew.sum(dim='d_gene')
daNew = mdlNew.xdata['d_gene']*mdlNew.xdata['j_choice'] # P(D|J) x P(J)
daNew = daNew.sum(dim='j_choice') # P(D) = \sum_J P(D|J) x P(J)
daNew = mdlNew.xdata['d_5_del']*daNew
daNew = daNew.sum(dim='d_gene')
daNew = mdlNew.xdata[strEvent]

mdlNew.xdata['d_5_del']

#daOld = mdlOld.xdata[strEvent]*mdlOld.xdata['d_gene']
#daOld = daOld.sum(dim='d_gene')
daOld = mdlOld.xdata['d_gene']*mdlOld.xdata['j_choice']
daOld = daOld.sum(dim='j_choice')
daOld = mdlOld.xdata[strEvent]*daOld
daOld = daOld.sum(dim='d_gene')


daOld.plot(ax=ax, marker = 'o', label="Current")
daNew.plot(ax=ax, marker = 'o', label="Inferred")

ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values) #, rotation=90)
fig.tight_layout()




########### INSERTIONS ###########
##### vd_ins #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'vd_ins'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"
mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#da = mdl.xdata[strEvent]
daNew = mdlNew.xdata[strEvent]

daOld = mdlOld.xdata[strEvent]

daOld.plot(ax=ax, marker = 'o', label="Current")
daNew.plot(ax=ax, marker = 'o', label="Inferred")

ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values) #, rotation=90)
fig.tight_layout()

fig.savefig(strEvent+"_compar.pdf")


##### dj_ins #####
fig, ax = plt.subplots(figsize=(12,8))

strEvent  = 'dj_ins'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"
mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#da = mdl.xdata[strEvent]
daNew = mdlNew.xdata[strEvent]

daOld = mdlOld.xdata[strEvent]

daOld.plot(ax=ax, marker = 'o', label="Current")
daNew.plot(ax=ax, marker = 'o', label="Inferred")

ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values) #, rotation=90)
fig.tight_layout()

fig.savefig(strEvent+"_compar.pdf")




mdlOld.xdata[strEvent]
mdlOld.xdata.keys()


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






