#!/usr/bin/python3
"""
Created on Tue Jan 14 16:33:37 2020

@author: alfaceor
"""
import importlib
import pygor3 as p3
importlib.reload(p3)

igor_specie="human"
igor_chain="tcr_beta"
mdl = p3.IgorModel.load_default(igor_specie, igor_chain)
#print(mdl)
#mdl.parms.plot_Graph()
print( mdl.get_events_nicknames_list() )

strEvent = 'd_3_del'
Grev = (mdl.parms.G.reverse())
GG = mdl.parms.G
GGG = GG.subgraph(GG.predecessors(strEvent))

GGG.nodes
GG
print(GG)
import networkx as nx
nx.draw(GGG, with_labels=True)
list(nx.isolates(GGG))

for node in GGG.nodes:
    print(node)
    pred = list( GGG.predecessors(node) )

    if len(pred) == 0:
        print("Ho!")
    else:
        print(pred)

print(mdl.parms.G.predecessors('d_3_del'))

# 1. Find the dependencies of the event
#print(Grev[strEvent])
def p_nota(str_ev:str, lista):
    return "P("+str_ev+"|"+str(lista)+")"

# FIXME: get_marginal should be a function of a list of strEvents.
def get_marginal(mdl, Grev, strEvent):
    lista =list( Grev.neighbors(strEvent) )
    len_lista = len(lista)
    str_marg = "P("+strEvent+") = "
    if len_lista == 0 :
        #print(mdl.xdata[strEvent])
        return str_marg+strEvent
    else:
        return str_marg+p_nota(strEvent, lista)

        # tmp_marg = p_nota(strEvent)
        # print(strEvent, 'lista : ', lista)
        # for str_ev in lista:
        #     tmp_marg = tmp_marg+ " * "+ p_nota(str_ev) +"*"+get_marginal(mdl, Grev, str_ev)
        # return tmp_marg


str_marginal = get_marginal(mdl, Grev, strEvent)
final_str = "P("+strEvent+")"
print(final_str+" = "+ str_marginal)
print(get_marginal(mdl, Grev, 'd_gene'))

"""
def get_marginal(mdl, Grev, strEvent):
    lista =list( Grev.neighbors(strEvent) ) 
    len_lista = len(lista)
    if len_lista == 0 :
        print(mdl.xdata[strEvent])
    else:
        print( len_lista, lista )

"""




#mdl.plot_Event_Marginal('v_choice')
# From the loaded model I want a list of event_types



