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

strEvent = 'v_choice'
strEvent = 'd_3_del'
Grev = (mdl.parms.G.reverse())

# 1. Find the dependencies of the event
#print(Grev[strEvent])
def p_nota(str_ev:str):
    return "P("+str_ev+"|)"

# FIXME: get_marginal should be a function of a list of strEvents.
def get_marginal(mdl, Grev, strEvent):
    lista =list( Grev.neighbors(strEvent) )
    len_lista = len(lista)
    if len_lista == 0 :
        #print(mdl.xdata[strEvent])
        return strEvent
    else:
        tmp_marg = ''
        print(strEvent, 'lista : ', lista)
        for str_ev in lista:
            tmp_marg = tmp_marg+ "\n,* "+p_nota(str_ev) +"*"+get_marginal(mdl, Grev, str_ev)
        return tmp_marg


str_marginal = get_marginal(mdl, Grev, strEvent)
final_str = "P("+strEvent+")"
print(final_str+" = "+ str_marginal)

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



