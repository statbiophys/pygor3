#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:16:35 2019

@author: alfaceor
"""

from Bio import SeqIO

IgorSpecie    = "mouse"
IgorChain     = "tcr_beta"
IgorModelPath = "../../IGoR/models/"+IgorSpecie+"/"+IgorChain+"/"
IgorRefGenomePath = IgorModelPath+"ref_genome/"

#strWD="uchicago/selected_sequences/"
#

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"

import IgorModel
mdl = IgorModel.Model(model_parms_file=flnModelParms, model_marginals_file=flnModelMargs)
mdl.parms.write_model_parms()

evento = mdl.parms.Event_list[0]

evento.seq_type

evento.event_typeIgorModel
evento.seq_type
evento.seq_side
evento.priority
evento.nickname

flnVGeneTemplate = IgorRefGenomePath+"genomicVs.fasta"
flnDGeneTemplate = IgorRefGenomePath+"genomicDs.fasta"
flnJGeneTemplate = IgorRefGenomePath+"genomicJs.fasta"




ev = IgorModel.Rec_Event('Genechoice', 'V_gene', 'Undefined_side', 7, 'v_choice')

ev.add_realization()
