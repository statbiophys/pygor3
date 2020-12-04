#Example: load default model
import pygor3 as p3

igor_specie="human"
igor_chain="tcr_beta"
mdl0 = p3.IgorModel.load_default(igor_specie, igor_chain)
mdl0.parms.plot_Graph()
mdl0.marginals

mdl0.xdata

mdl0.plot_Event_Marginal('v_choice')
mdl0.plot_Event_Marginal('j_choice')
mdl0.plot_Event_Marginal('d_gene')


mdl0.xdata['v_choice']

mdl0.xdata['v_choice'].plot()
mdl0.plot_Event_Marginal('v_choice')
mdl0.xdata['j_choice'].plot()
mdl0.plot_Event_Marginal('j_choice')
mdl0.xdata['d_gene'][{'j_choice' : 6}].plot()
mdl0.xdata['d_gene'][{'j_choice' : 6, 'v_choice' : 0}].plot()


mdl = p3.IgorModel()
mdl.parms.G = mdl0.parms.G
mdl.parms.plot_Graph()
print( list( mdl.parms.G ) )

vj_dict_default = {
    'v_choice': {},
    'j_choice': {'v_choice'},
    'd_gene': {'v_choice', 'j_choice'},
    'v_3_del': {'v_choice'},
    'd_5_del': {'d_gene'},
    'd_3_del': {'d_gene', 'd_5_del'},
    'j_5_del': {'j_choice'},
    'vd_ins': {},
    'vd_dinucl': {},
    'dj_ins': {},
    'dj_dinucl': {}
}

vj_dict_parms = {
    'v_choice': {},
    'j_choice': {'v_choice'},
    'd_gene': {'v_choice', 'j_choice'},
    'v_3_del': {'v_choice'},
    'd_5_del': {'d_gene'},
    'd_3_del': {'d_gene', 'd_5_del'},
    'j_5_del': {'j_choice'},
    'vd_ins': {},
    'vd_dinucl': {},
    'dj_ins': {},
    'dj_dinucl': {}
}

# Given a dictionary create Bayes network
mdl1_parms = p3.IgorModel_Parms.from_network_dict(vj_dict_default)
# this should generate a minimal default model.
# And what I mean with that? Event_list without realizations
print(mdl1_parms.Event_list[0])
# Now add the realizations to the in Parms class
mdl1_parms.Event_dict['v_choice']

# So the idea is that
pathGenomic = p3.rcParams['paths.igor_models']+ "/human/tcr_beta/"+"ref_genome/"

vj_dict_parms = {
    'v_choice': {"filename" : pathGenomic+"genomicVs.fasta"},
    'j_choice': {"filename" : pathGenomic+"genomicJs.fasta"},
    'd_gene': {"filename" : pathGenomic+"genomicDs.fasta"},
    'v_3_del': {"limits" :(-4,20)},
    'd_5_del': {"limits" : (-4,20)},
    'd_3_del': {"limits" : (-4,20)},
    'j_5_del': {"limits" : (-4,20)},
    'vd_ins': {"limits" : (0,24)},
    'vd_dinucl': {},
    'dj_ins': {"limits" : (0,24)},
    'dj_dinucl': {}
}



mdl1_parms.load_from_dict()
mdl1_parms.load_GeneChoice_realizations_by_nickname('v_choice', pathGenomic+"genomicVs.fasta")
mdl1_parms.load_GeneChoice_realizations_by_nickname('j_choice', pathGenomic+"genomicJs.fasta")
mdl1_parms.load_GeneChoice_realizations_by_nickname('d_gene', pathGenomic+"genomicDs.fasta")
mdl1_parms.load_Deletion_realizations_by_nickname('v_3_del')
mdl1_parms.load_Deletion_realizations_by_nickname('d_5_del')
mdl1_parms.load_Deletion_realizations_by_nickname('d_3_del')
mdl1_parms.load_Deletion_realizations_by_nickname('j_5_del')
mdl1_parms.load_Insertion_realizations_by_nickname('vd_ins')
mdl1_parms.load_Insertion_realizations_by_nickname('dj_ins')
mdl1_parms.load_DinucMarkov_realizations_by_nickname('vd_dinucl')
mdl1_parms.load_DinucMarkov_realizations_by_nickname('dj_dinucl')

ev = mdl1_parms.get_Event('v_3_del')
ev.realizations
list( map(lambda x: x.nickname, mdl1_parms.Event_list) )
print( mdl1_parms.get_Event('v_choice') )

ev = mdl0.parms.get_Event('vd_dinucl')
print(ev)
ev.realizations
ev.get_realization_vector()

print(mdl0.parms.get_Event('v_choice') )


mdl0.parms.get_Event('vd_dinucl').realizations