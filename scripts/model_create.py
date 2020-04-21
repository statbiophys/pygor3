#!/usr/bin/env python3
import pygor3 as p3

from  optparse import OptionParser
def main():
    parser = OptionParser()
    parser.add_option("-s", "--species", dest="species", help='Igor species')
    parser.add_option("-c", "--chain", dest="chain", help='Igor chain')
    parser.add_option("-b", "--batch", dest="batch", help='Batchname to identify run. If not set random name is generated')
    parser.add_option("-o", "--output", dest="output", help='filename of csv file to export data')

    (options, args) = parser.parse_args()

    igor_specie="human"
    igor_chain="tcr_beta"
    mdl0 = p3.IgorModel.load_default(igor_specie, igor_chain)
    #mdl0.parms.plot_Graph()
    #print( p3.Igor_vj_dict_default)
    #print( p3.IgorRec_Event_default_dict.keys() )
    # load Model_Parms

    mdl_parms = p3.IgorModel_Parms()
    aaa = p3.IgorRec_Event_default_dict['v_choice']
    print(aaa)
    event = p3.IgorRec_Event.from_dict(aaa)
    mdl_parms.Event_list.append(event)

    mdl_parms.Event_list.append(p3.IgorRec_Event.from_dict(p3.IgorRec_Event_default_dict['v_choice']))
    print(mdl_parms.Event_list) #load_events_from_dict({})
    for ev in mdl_parms.Event_list:
        print(ev)
    # evento = mdl_parms.Event_list[0]
    # print(evento)

    """
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
    """


if __name__ == "__main__":
    main()