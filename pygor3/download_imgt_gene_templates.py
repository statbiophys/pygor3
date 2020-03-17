import importlib
import pygor3 as p3
importlib.reload(p3)

# get data from IMGT and generate a model

# print( p3.imgt.imgt_params )
species = p3.imgt.get_species_list()
print(species)
flnGenome = p3.imgt.download_gene_template('Homo+sapiens', 'TRBV')

igor_specie="human"
igor_chain="tcr_beta"
mdl = p3.IgorModel.load_default(igor_specie, igor_chain)

mdl.parms.plot_Graph()
mdl.get_events_nicknames_list()
mdl.plot_Event_Marginal('v_3_del')
mdl.plot_Event_Marginal('v_choice')

mdl.xdata['v_choice']


#mdl.plot_Event_Marginal('v_choice')
# From the loaded model I want a list of event_types


"""
mdl.xdata['d_gene']
mdl.xdata['j_choice']

mdl.xdata['d_gene'].dot( mdl.xdata['j_choice'] )

mdl.xdata['d_3_del'].shape
mdl.xdata['d_3_del'].coords
mdl.xdata['d_3_del'].dims
mdl.xdata['d_3_del']*mdl.xdata['d_5_del']
(mdl.xdata['d_3_del'].dot(mdl.xdata['d_5_del'])).dims
mdl.parms.G['d_3_del']

elist = mdl.parms.Event_list
print(str(elist[0]))
nelist = sorted(elist, key=lambda x:x.priority, reverse=False)

mdl.parms.dictNicknameName.keys()
#0 Define the templates and define model.
# p3.imgt.download_gene(specie, gene)


#0.1 Get templates from IMGT from a list of available templates.
#1 Choose an specie and chain, this means define a model.

iseq = p3.IgorIndexedSequence(10,'ATCTGAGCT')
print(iseq)


mdl = p3.IgorModel.load_default("mouse", "tcr_beta")


#2 Input a sequence
#3 Align the sequence given the model, show the alignment
#3 Give me the pgen of the sequence probability of this sequence.


bs = p3.IgorBestScenariosVDJ()
print(bs)

bs.
"""