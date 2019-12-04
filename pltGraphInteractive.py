#!/usr/bin/python
import networkx as nx

# 1. Load IGoR data
import matplotlib.pyplot as plt
#fig, ax = plt.subplots()

import numpy as np
import pygor

batchname = "TRbeta"
flnIgorDB = "chicagoMouse.db"
strWD="uchicago/"
flnIgorIndexedSeq = strWD+"aligns/"+batchname+"_indexed_sequences.csv"

IgorSpecie    = "mouse"
IgorChain     = "tcr_beta"
IgorModelPath = "../../IGoR/models/"+IgorSpecie+"/"+IgorChain+"/"
IgorRefGenomePath = IgorModelPath+"ref_genome/"

flnModelParms = IgorModelPath + "models/model_parms.txt"
flnModelMargs = IgorModelPath + "models/model_marginals.txt"
flnIgorBestScenarios = strWD+batchname+"_output/best_scenarios_counts.csv"




model = pygor.models.genmodel.GenModel(model_parms_file=flnModelParms, marginals_file=flnModelMargs)
#model = pygor.models.genmodel.GenModel(model_parms_file="model_parms.txt", marginals_file= "model_marginals.txt")
#model.marginals

model.edges['Deletion_D_gene_Five_prime_prio5_size21'].parents

#model.get_event('v_choice', by_nickname=True)
for recEv in model.events:
  print recEv.nickname, recEv.priority
#recomEvent.priority

probData = model.marginals[0]
probData['v_choice']
# 2. Create a Graph
# Get the graph properties.
graphData = model.marginals[1]
# 2.1. Get the nodes
nodes = graphData.keys() # Give me how is the netwok.
# 2.2. Get the edges
graphData
G = nx.DiGraph()


G.nodes
G.nodes.data()
for node in nodes:
  G.add_node(node)
#print G.number_of_nodes()

for node in nodes:
  for node2link in graphData[node]:
    G.add_edge(node2link, node)
    #print node, node2link


# 3. Plot Graph
pos = nx.spring_layout(G)


# lo que quiero es algo asi:
# prio_dict = {prio1: [event1, event2], prio2: [event3, event4, event5], prio3:[evento6]}
prio_dict = dict()
for event in model.events:
  if not (event.priority in prio_dict):
    prio_dict[event.priority] = list()
  prio_dict[event.priority].append(event)
  
print prio_dict
# ahora definir la posicion del coso este
xwidth = 240
yfactor= 40
for key in prio_dict:
  lenKey = len(prio_dict[key])
  if lenKey == 1:
    pos[prio_dict[key][0].nickname] = np.array([float(xwidth)/2.0, float(key)*yfactor])
  else:
    xx = np.linspace(0,xwidth,lenKey)
    for ii, ev in enumerate(prio_dict[key]):
      xpos = xx[ii] #float(xwidth)*float(ii)/float(lenKey)
      pos[ev.nickname] = np.array([xpos, float(key)*yfactor])


print "*"*20
print pos



"""

#def get_event_same_prio(events):
#  #model.events
prioNicknameList = []
for event in model.events:
  prioNicknameList.append([float(event.priority), event])

# ordered list by priority
sortedList = sorted(prioNicknameList, key=lambda x : x[0])

sortedList
model.events

lastPrio = 0 #sortedList[0][0]
xscaleFactor=1
yscaleFactor=10
xpos0 = 0
xpos  = xpos0
for a in sortedList:
  currPrio = a[0]
  if currPrio == lastPrio :
    xpos = xpos + xscaleFactor
  else:
    xpos = xpos0
    lastPrio = currPrio

  pos[a[1].nickname] = np.array([xpos, a[0]*yscaleFactor])

#print pos.keys()
# center the graph
xvalues = []
yset = set()
for key in pos.keys():
  xvalues.append(pos[key][0])
  yset.add(pos[key][1])
xvalues = np.array(xvalues)
print xvalues, xvalues.max()
# 2. know look for the guys with the same priority or y position
y_dict=dict()
for yy in yset:
  y_dict[yy] = list() #yset.add(pos[key][1])
print y_dict

for aa in y_dict:
  for key in pos:
    if aa == pos[key][1]:
      y_dict[aa].append(key)

print y_dict
for aa in y_dict:
  nn = len(y_dict[aa])
  if nn ==1 :
    pos[y_dict[aa][0]][0] = float(xvalues.max())/2.
  else:
    for eventlbl in y_dict[aa]:
      pos[y_dict[aa]][0] = pos[y_dict[aa]][0]*float(xvalues.max())
"""

#y_dict= dict()
#for key in pos.keys():
#    if pos[key][1] = 

"""
{'d_5_del': array([ 2., 50.]), 'vd_dinucl': array([ 0., 30.]), 'v_choice': array([ 0., 70.]), 'dj_ins': array([ 0., 20.]), 'j_choice': array([ 1., 70.]), 'vd_ins': array([ 0., 40.]), 'd_gene': array([ 0., 60.]), 'd_3_del': array([ 1., 50.]), 'v_3_del': array([ 0., 50.]), 'dj_dinucl': array([ 0., 10.]), 'j_5_del': array([ 3., 50.])}
"""



fig, ax = plt.subplots()
ax.set_aspect('equal')
#nx.draw(G, pos, with_labels=True, font_weight='bold', nodesize=2000) # FIXME: make a better plot: cutting edges.
nx.draw(G, pos, ax=ax, with_labels=True, arrowsize=20,
        node_size=800, font_size=10, font_weight='bold') # FIXME: make a better plot: cutting edges.

plt.show()
#nx.draw_circular(G, with_labels=True, arrowsize=20, arrowstyle='fancy', node_size=4000, font_weight='bold') # FIXME: make a better plot: cutting edges.


# So now I want to know what happend when I choose a recombination event
# For instance for a given priority If I choose v_choice

#pickedEvent='v_choice'
#probData[pickedEvent]
#for i in probData[pickedEvent]:
#  print i
#  #dim = len(i)
#  #print dim
#  if True:
#    fig2, ax2 = plt.subplots()
#
#    ax2.set_title(pickedEvent)
#    ax2.set_ylabel("P_event")
#    #ax2.set_ylim(1e-4,1)
#    ax2.set_yscale("log")
#    ax2.plot(probData[pickedEvent], '-o')
#    ordIndx = np.argsort(probData[pickedEvent])[::-1]
#    ax2.plot(probData[pickedEvent][ordIndx], '-o')
#
#
#
#plt.show()


