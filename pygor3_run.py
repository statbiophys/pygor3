#!/home/alfaceor/anaconda3/envs/pygor3/bin/python
import numpy as np
import subprocess

import pygor3
import panel as pn
pn.extension()

Titulo = "<h1> Run IGoR<\h1>"
#igor_specie="human"disabled
#igor_chain="tcr_beta"
igor_species_list = ["human", "mouse"]

igor_option_path_dict={
    "alpha": "tcr_alpha", 
    "beta" : "tcr_beta", 
    "light": "light", 
    "heavy_naive" : "bcr_heavy", 
    "heavy_memory": "bcr_heavy"
}

import subprocess
defaults=dict()
defaults['igor_exec_path'] = subprocess.check_output(["which", "igor"]).decode("utf-8").replace('\n', '')
#defaults['igor_data_dir'] = subprocess.check_output(["which", "igor"]).decode("utf-8").replace('\n', '')

wd_igor_exec = pn.widgets.TextInput(placeholder=defaults['igor_exec_path'], name="Igor exec path", disabled=True)

wd_igor_wd = pn.widgets.TextInput(placeholder=defaults['igor_exec_path'], name="Working directory")
wd_igor_batchname = pn.widgets.TextInput(name="Batchname")

wd_igor_input_seqs = pn.widgets.TextInput(name="Input sequences path")

wd_igor_specie = pn.widgets.Select(options=igor_species_list, name="Specie")
wd_igor_chain = pn.widgets.Select(options=igor_option_path_dict, name="Chain")
wd_igor_Run = pn.widgets.Button(name="Infer", button_type='success')

pn_header = Titulo
pn_first = pn.Row(wd_igor_exec, wd_igor_batchname, wd_igor_wd)
pn_input_sequences = pn.Row(wd_igor_input_seqs)
pn_model_selection = pn.Column(wd_igor_specie, wd_igor_chain)
pn_last = pn.Row(wd_igor_Run)
layout = pn.Column(pn_header, pn_first, pn_input_sequences, pn_model_selection, pn_last)


"""
#subprocess.call(command, shell=True)
def InferRun():
    igor_command="echo "


    
def showOutput():
    #update layout
    
"""
layout.show()



