#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 08:56:57 2019

@author: alfaceor
Compare insertions
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


########### INSERTIONS ###########
##### vd_ins #####
fig, ax = plt.subplots(figsize=(12,8))

yEtiqueta = "$ P $"
IgorModelInferencePath = strWD + batchname + "_inference/"

#################################
flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
mdlOld = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#--- 1. Load Old model events
strEvent  = 'vd_ins'
lblEvent  = strEvent.replace("_", " ")
daOld = mdlOld.xdata[strEvent]
daOld.plot(ax=ax, marker = 'o', label="Current, "+lblEvent)

strEvent  = 'dj_ins'
lblEvent  = strEvent.replace("_", " ")
daOld = mdlOld.xdata[strEvent]
daOld.plot(ax=ax, marker = 'o', label="Current, "+lblEvent)

#################################
flnModelParms = IgorModelInferencePath + "final_parms.txt"
flnModelMargs = IgorModelInferencePath + "final_marginals.txt"
mdlNew = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)

#--- 2. Load New model events.
strEvent  = 'vd_ins'
lblEvent  = strEvent.replace("_", " ")
daNew = mdlNew.xdata[strEvent]
daNew.plot(ax=ax, marker = 'x', label="Inferred, "+lblEvent)

strEvent  = 'dj_ins'
lblEvent  = strEvent.replace("_", " ")
daNew = mdlNew.xdata[strEvent]
daNew.plot(ax=ax, marker = 'x', label="Inferred, "+lblEvent)


#--- 3. Figure details
ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xticks(mdlOld.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlOld.xdata[strEvent]['lbl__'+strEvent].values) #, rotation=90)
xEtiqueta = "Insertions"
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

fig.tight_layout()



##### dj_ins #####
#fig, ax = plt.subplots(figsize=(12,8))

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

daOld.plot(ax=ax, marker = 'o', label="Current, "+lblEvent)
daNew.plot(ax=ax, marker = 'o', label="Inferred, "+lblEvent)

ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)
#mdl.xdata[strEvent]['lbl__'+strEvent]
ax.set_xticks(mdlNew.xdata[strEvent][strEvent].values)
ax.set_xticklabels(mdlNew.xdata[strEvent]['lbl__'+strEvent].values) #, rotation=90)
fig.tight_layout()

fig.savefig("Insertions_compar.pdf")















ax.legend(bbox_to_anchor=(0.8, 0.72), loc="lower left", prop={'size':15})
ax.set_xlabel(xEtiqueta)
ax.set_ylabel(yEtiqueta)

