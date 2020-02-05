#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 11:01:59 2019

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

fig.savefig(strEvent+".pdf")
mdl.xdata[strEvent] #.plot(marker = 'o')



########### DELETIONS ###########
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

fig.savefig(strEvent+".pdf")
mdl.xdata[strEvent] #.plot(marker = 'o')









"""

fig, ax = plt.subplots(figsize=(8,6))

strEvent  = 'v_choice'
lblEvent  = strEvent.replace("_", " ")
xEtiqueta = lblEvent
yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

flnModelParms = IgorModelInferencePath + "initial_model.txt"
flnModelMargs = IgorModelInferencePath + "initial_marginals.txt"
mdl = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

mdl.xdata[strEvent].plot(ax=ax, marker = 'o', label="Initial")
mdl.xdata[strEvent].sum()
for i in [1, 2, 3, 4, 5]:
    flnModelParms = IgorModelInferencePath + "iteration_"+str(i)+"_parms.txt"
    flnModelMargs = IgorModelInferencePath + "iteration_"+str(i)+".txt"
    
    mdl = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
    mdl.xdata[strEvent].plot(ax=ax, marker = 'o', label=str(i))

ax.legend()
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdl.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdl.xdata[strEvent]['lbl__'+strEvent].values, rotation=90)

fig.savefig(strEvent+".pdf")
mdl.xdata[strEvent] #.plot(marker = 'o')

"""




# Genechoice
