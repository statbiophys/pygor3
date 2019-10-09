#      Author: Quentin Marcou & Carlos Olivares
#
#  This source code is distributed as part of the IGoR software.
#  IGoR (Inference and Generation of Repertoires) is a versatile software to
#  analyze and model immune receptors generation, selection, mutation and all
#  other processes.
#   Copyright (C) 2019  Quentin Marcou & Carlos Olivares
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

import numpy as np
import pandas as pd
import xarray as xr
import networkx as nx

#import networkx as nx

# Generation of label for a simple identification of genomic template sequence.
def genLabel(strName):
    aaa = strName.split("|")
    if len(aaa) > 1 :
        return aaa[1]
    else:
        return strName


class Model:
    def __init__(self, model_parms_file=None, model_marginals_file=None):
        self.parms     = Model_Parms()
        self.marginals = Model_Marginals()
        self.xdata      = dict()

        # check input files
        flag_parms     = (model_parms_file is not None)
        flag_marginals = (model_marginals_file is not None)
        flag_xdata     = (flag_parms and flag_marginals)

        if flag_parms:
            self.parms.read_model_parms(model_parms_file)
        if flag_marginals:
            self.marginals.read_model_marginals(model_marginals_file)
        if flag_xdata:
            self.generate_xdata()

    def generate_xdata(self):
        for key in self.marginals.marginals_dict:
            self.xdata[key] = xr.DataArray(self.marginals.marginals_dict[key], dims=tuple(self.marginals.network_dict[key]))
            #print "key: ", key, self.xdata[key].dims
            strCoord = "priority"
            self.xdata[key][strCoord] = self.parms.get_Event(key).priority
            for strDim in self.xdata[key].dims:
                self.xdata[key][strDim] = range(len(self.xdata[key][strDim]))
                if strDim in ['v_choice', 'j_choice', 'd_gene']:
                    #print strDim
                    labels = self.parms.Event_dict[strDim]['name'].map(genLabel).values
                    #print type(labels)
                    strCoord = 'lbl__'+strDim
                    self.xdata[key][strCoord] = (strDim, labels) # range(len(self.xdata[key][coord]))
                elif not (strDim in ['vd_dinucl', 'dj_dinucl', 'vj_dinucl'] ):
                    labels = self.parms.Event_dict[strDim]['value'].values
                    #print strDim
                    #print labels
                    strCoord = 'lbl__'+strDim
                    self.xdata[key][strCoord] = (strDim, labels) # range(len(self.xdata[key][coord]))
                    

class Model_Parms:
    """
    Class to get a list of Events directly from the *_parms.txt
    """
    def __init__(self, model_parms_file=None):
        self.Event_list = list() # list of Rec_event
        self.Edges      = list()
        self.ErrorRate  = list()
        self.Event_dict = dict()
        self.Edges_dict = dict()
        self.dictNameNickname=dict()
        self.dictNicknameName=dict()
        self.G          = nx.DiGraph()
        self.preMarginalDF = pd.DataFrame()

        if model_parms_file is not None:
            print(model_parms_file)
            self.read_model_parms(model_parms_file)
            #self.get_EventDict_DataFrame()

    def __str__(self):
        tmpstr = "{ 'len Event_list': "+str(len(self.Event_list)) \
                +", 'len Egdes': "+str(len(self.Edges)) \
                +", 'len ErrorRate': "+str(len(self.ErrorRate))+" }"
        return tmpstr
        #return "{ Event_list, Egdes, ErrorRate}"
        
    
#    def __eq__(self, other):
#        if isinstance(self, other.__class__):
#            return (self.Event_list == other.Event_list) and (self.Edges == other.Edges) and (self.ErrorRate == other.ErrorRate)
#        else:
#            return NotImplemented
#
#    def __hash__(self):
#        return 3;
    def read_model_parms(self, filename):
        """Reads a model graph structure from a model params file.
        Note that for now this method does not read the error rate information.
        """
        with open(filename, "r") as ofile:
            # dictionary containing recarrays?
            line = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)

            if strip_line == "@Event_list":
                self.read_Event_list(ofile)

            line = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
            if strip_line == "@Edges":
                self.read_Edges(ofile)

            # FIXME: ErrorRate added
            line = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
            if strip_line == "@ErrorRate" :
                self.read_ErrorRate(ofile)
        self.get_EventDict_DataFrame()

    # save in Event_list
    def read_Event_list(self, ofile):
        lastPos    = ofile.tell()
        line = ofile.readline()
        strip_line = line.rstrip('\n')  # Remove end of line character
        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        #event = Rec_Event()

        while strip_line[0] == '#':
            # get the metadata of the event list
            event_metadata = strip_line[1:].split(";") #GeneChoice;V_gene;Undefined_side;7;v_choice
            event_metadata[3] = int(event_metadata[3]) # change priority to integer
            event = Rec_Event( *event_metadata )
            #self.G.add_node(event.nickname)
            # Now read the realizations (or possibilities)
            lastPos    = ofile.tell()
            line       = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)


            while strip_line[0] == '%':
                realization = Event_realization()
                realizData  = strip_line[1:].split(";")
                if event.event_type == "GeneChoice":
                    realization.name  = realizData[0]
                    realization.value = realizData[1]
                    realization.index = int(realizData[2])
                elif event.event_type == "DinucMarkov":
                    realization.value = realizData[0]
                    realization.index = int(realizData[1])
                else:
                    realization.value = int(realizData[0])
                    realization.index = int(realizData[1])

                event.add_realization(realization)
                # next line
                lastPos  = ofile.tell()
                line     = ofile.readline()
                strip_line = line.rstrip('\n')  # Remove end of line character
                strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)

            self.Event_list.append(event)
        ofile.seek(lastPos)


    def read_Edges(self, ofile):
        #print "read_Edges"
        lastPos  = ofile.tell()
        line = ofile.readline()
        strip_line = line.rstrip('\n')  # Remove end of line character
        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        while strip_line[0] == '%':
            edge = strip_line[1:].split(';')
            self.Edges.append(edge)
            #read nextline
            lastPos  = ofile.tell()
            line = ofile.readline()
            strip_line = line.rstrip('\n')  # Remove end of line character
            strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        ofile.seek(lastPos)


    def read_ErrorRate(self, ofile):
        lastPos  = ofile.tell()
        line = ofile.readline()
        strip_line = line.rstrip('\n')  # Remove end of line character
        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
        while strip_line[0] == '#':
            if 'SingleErrorRate' == strip_line[1:] :
                lastPos  = ofile.tell()
                line = ofile.readline()
                strip_line = line.rstrip('\n').rstrip()
                error = strip_line
                self.ErrorRate = {"SingleErrorRate" : error }
        ofile.seek(lastPos)

    def get_Event(self, event_nickname, by_nickname=True):
        """Returns the RecEvent with corresponding name or nickname."""
        if by_nickname:
            for ev in self.Event_list:
                if ev.nickname == event_nickname:
                    return ev
            raise Exception(
                'RecEvent with nickname \"' + event_nickname + "\" not found.")
        else:
            for ev in self.Event_list:
                if ev.name == event_nickname:
                    return ev
            raise Exception(
                'RecEvent with name \"' + event_nickname + "\" not found.")

    def get_EventDict_DataFrame(self):
        self.Event_dict = dict()
        self.dictNameNickname = dict()
        #dictio = dict()
        for event in self.Event_list:
            #dictio[event.nickname] = event.get_realization_DataFrame()
            self.Event_dict[event.nickname] = event.get_realization_DataFrame()
            self.dictNameNickname[event.name] = event.nickname
        #return dictio
        self.dictNicknameName = {v: k for k, v in self.dictNameNickname.items()}
        self.getBayesGraph()
    
    def getBayesGraph(self):
        self.G = nx.DiGraph()
        for rec_event in self.Event_list:
            self.G.add_node(rec_event.nickname)
            self.Edges_dict[rec_event.nickname] = list()
            self.dictNameNickname[rec_event.name] = rec_event.nickname
        
        for edge in self.Edges:
            # Graph to get the dependecies
            self.G.add_edge(self.dictNameNickname[edge[0]], self.dictNameNickname[edge[1]])
            self.Edges_dict[self.dictNameNickname[edge[1]]].append(self.dictNameNickname[edge[0]])
        self.G = self.G.reverse()
        
    def genPreMarginalDF(self):
        data=[]
        for event in self.Event_list:    
            #print (parms.dictNameNickname[event.name])
            #parms.Edges
            #tmpDict = dict()
            lista = []
            for edge in self.Edges:
                if edge[1] == event.name:
                    #print(parms.dictNameNickname[edge[0]])
                    #print(edge[0])
                    lista.append(self.dictNameNickname[edge[0]])
            tmpDict = {'event': event.nickname, 'priority': event.priority, 'Edges': lista }
            data.append(tmpDict)
        
        self.preMarginalDF = pd.DataFrame(data) #.set_index('event')
        self.preMarginalDF['nEdges'] = self.preMarginalDF['Edges'].map(len)
        self.preMarginalDF.sort_values(['priority', 'nEdges'], ascending=[False,True])


    def genMarginalFile(self, model_marginals_file=None):
        self.genPreMarginalDF()
        #self.preMarginalDF
        if model_marginals_file == None:
            model_marginals_file = "model_marginals.txt"
        ofile = open(model_marginals_file, "w")
        for index, row in self.preMarginalDF.iterrows():
            nickname = row['event']
            ofile.write("@"+nickname+"\n")
            #DimEvent = len(parms.Event_dict[event.nickname])
            #DimEdges = len(parms.Edges_dict[event.nickname])
            
            DimEvent = len(self.Event_dict[nickname])
            strDimLine = "$Dim["
            DimList = []
            if row['nEdges'] == 0:
                strDimLine = strDimLine +str(DimEvent)
                strDimLine = strDimLine +"]"
            else:
                for evNick in row['Edges']: #parms.Edges_dict[event.nickname]:
                    Dim = len(self.Event_dict[evNick])
                    strDimLine = strDimLine +str(Dim)+","
                    DimList.append(Dim)
                strDimLine = strDimLine + str(DimEvent)
                strDimLine = strDimLine +"]"
            ofile.write(strDimLine+"\n")
            
            
            lista = row['Edges'] # self.Event_dict[nickname]
            for indices in np.ndindex(tuple(DimList)):
                #print indices
                strTmp = "#"
                for ii in range(len(lista)):
                    strTmp = strTmp+"["+lista[ii]+","+str(indices[ii])+"]"
                    if not (ii == len(lista)-1):
                        strTmp = strTmp + ","
                ofile.write(strTmp+"\n")
                ofile.write("%")
                unifProb = (1./DimEvent)
                for jj in range(DimEvent):
                    ofile.write(str(unifProb))
                    if not (jj == DimEvent-1):
                        ofile.write(",")
                ofile.write("\n")
                    
            #print (parms.G[event.nickname])
            #print (len(parms.G[event.nickname]))
            
            # TODO: 
#            if len(DimList) == 0:
#                ofile.write("#\n")
#                # TODO: Escribe todas las prob.
#            else:
#                for dd in range(len(DimList)):
#                    for 'v_choice' == 0 range(DimList[dd]):
#                        for 'j_choice' in 1 to len('j_choice'):
#                    print(dd)
#                for evNick in row['Edges']: #parms.Edges_dict[event.nickname]:
#                    # TODO: Now make the list for each subevent
#                    for ii in range(len(self.Event_dict[evNick]))
#                    ofile.write("#["+evNick+","+str()+"]\n")
#                    Dim = len(self.Event_dict[evNick])
#                    strDimLine = strDimLine +","+str(Dim)
#                strDimLine = strDimLine +"]"
#            ofile.write(strDimLine+"\n")
#            
                               
        ofile.close()


class Rec_Event:
    """Recombination event class containing event's name, type, realizations,
    etc... Similar to IGoR's C++ RecEvent class.
    """
    def __init__(self, event_type, seq_type, seq_side, priority,
                 nickname):
        self.event_type = event_type
        self.seq_type = seq_type
        self.seq_side = seq_side
        self.priority = priority
        self.realizations = list()
        self.name = ""
        self.nickname = nickname
#        if nickname is not None:
#            self.nickname = nickname
        self.update_name()
        
    def __lt__(self, other):
        return self.priority < other.priority 
#        if ( self.priority < other.priority ):
#            return True
#        elif ( self.priority == other.priority  ):
#            # FIXME: less dependencies should be on top
            



#    def __eq__(self, other):
#        if isinstance(self, other.__class__):
#            return (    (self.event_type   == other.event_type) \
#                    and (self.seq_type     == other.seq_type) \
#                    and (self.seq_side     == other.seq_side) \
#                    and (self.priority     == other.priority) \
#                    and (self.realizations == other.realizations) \
#                    and (self.name         == other.name) \
#                    and (self.nickname     == other.nickname) \
#                    )
#        else:
#            return NotImplemented
#        
#    def __hash__(self):
#        return hash((self.event_type, \
#                     self.seq_type, \
#                     self.seq_side, \
#                     self.priority, \
#                     self.realizations, \
#                     self.name, \
#                     self.nickname))

#    def __str__(self):
#        return self.event_type+";"+self.seq_type+";"+self.seq_side+";"+self.priority+";"+self.nickname
#
#    def __repr__(self):
#        return "Rec_event(" + self.nickname + ")"

    def add_realization(self, realization):
        """Add a realization to the RecEvent realizations list."""
        self.realizations.append(realization)
        self.update_name()

    def update_name(self):
        """Updates the name of the event (will have no effect if the RecEvent
        has not been modified since the last call).

        """
        if self.event_type == "DinucMarkov":
            self.name = self.event_type + "_" + self.seq_type + "_" + \
                        self.seq_side + "_prio" + \
                        str(self.priority) + "_size" + \
                        str(len(self.realizations) ** 2)
        else:
            self.name = self.event_type + "_" + self.seq_type + "_" + \
                        self.seq_side + "_prio" + \
                        str(self.priority) + "_size" + \
                        str(len(self.realizations))

    def get_realization_vector(self):
        """This methods returns the event realizations sorted by the
        realization index as a list.

        """
        if self.event_type == 'GeneChoice':
            tmp = [""] * len(self.realizations)  # empty(, dtype = str)
        else:
            tmp = np.empty(len(self.realizations),
                              dtype=type(self.realizations[0].value))
        # print("Unfinished method get realization vector")
        processed_real_indices = []
        for real in self.realizations:
            if processed_real_indices.count(real.index) == 0:
                if real.name != "":
                    tmp[real.index] = real.name
                else:
                    tmp[real.index] = real.value
                processed_real_indices.append(real.index)
            else:
                print("REALIZATION INDICES ARE DEGENERATE")

        return tmp

    def get_realization_DataFrame(self):
        """ return an Event realizations as a pandas DataFrame to manipulate it.

        """
        return pd.DataFrame.from_records([realiz.to_dict() for realiz in self.realizations], index='index').sort_index()



class Event_realization:
    """A small class storing for each RecEvent realization its name, value and
    corresponding index.

    """
    def __init__(self):
        self.name  = "" #name
        self.value = "" #value
        self.index = "" #index

    def __lt__(self, other):
        return self.index < other.index

    def __gt__(self, other):
        return self.index > other.index

    def __str__(self):
        if self.name == "":
            return self.value+";"+self.index
        else:
            return self.name+";"+self.value+";"+str(self.index)

    def __repr__(self):
        return "Event_realization(" + str(self.index) + ")"

    def to_dict(self):
        return {
            'index' : self.index,
            'value' : self.value,
            'name'  : self.name
            }




### FIXME:
# @recombinationEvent
# $Dim
# #Indices of the realizations of the parent events
# %1d probability array.


class Model_Marginals:
    """
    Class to get a list of Events directly from the *_parms.txt
    """
    def __init__(self, model_marginals_file=None):
#        self.Event_list = list() # list of Rec_event
#        self.Edges      = list()
#        self.error_rate = list()
        self.marginals_dict = {}
        self.network_dict = {}
        #self.G          = nx.DiGraph()
        if model_marginals_file is not None:
            self.read_model_marginals(model_marginals_file)

    #  @d_3_del
    #  $Dim[3,21,21]
    #  #[d_gene,0],[d_5_del,0]
    #  %0,0,0,1.6468e-08,0.00482319,1.08101e-09,0.0195311,0.0210679,0.0359338,0.0328678,2.25686e-05,4.97463e-07,0,9.31048e-08,1.01642e-05,0.000536761,0.0260845,0.0391021,0.319224,0.289631,0.211165
    #  #[d_gene,0],[d_5_del,1]
    #  %0,0,6.86291e-08,2.00464e-09,0.00163832,2.02919e-06,0.0306066,0.0126832,0.000872623,0.016518,0.00495292,0.000776747,4.45576e-05,0.000667902,0.00274004,0.00435049,0.300943,0.182499,0.13817,0.302534,0

    def read_model_marginals(self, filename, dim_names=False):
        """Reads a model marginals file. Returns a tuple containing a dict
        containing the individual events probabilities indexed by the events
        nicknames and a dict containing the list of dimension names/ordering for
        each event.
        """
        with open(filename, "r") as ofile:
            # Model parameters are stored inside a dictionnary of ndarrays
#            marginals_dict = {}
#            network_dict = {}
            element_name = ""
            first = True
            first_dim_line = False
            element_marginal_array = []
            indices_array = []

            for line in ofile:
                strip_line = line.rstrip("\n")  # Remove end of line character
                if strip_line[0] == "@":
                    first_dim_line = True
                    if not first:
                        # Add the previous to the dictionnary
                        self.marginals_dict[element_name] = element_marginal_array
                    else:
                        first = False

                    element_name = strip_line[1:]
                # print element_name

                if strip_line[0] == "$":
                    # define array dimensions
                    coma_index = strip_line.find(",")
                    dimensions = []

                    # Get rid of $Dim[
                    previous_coma_index = 4
                    while coma_index != -1:
                        dimensions.append(
                            int(strip_line[previous_coma_index + 1:coma_index]))
                        previous_coma_index = coma_index
                        coma_index = strip_line.find(",", coma_index + 1)

                    # Add last dimension and get rid of the closing bracket
                    dimensions.append(int(strip_line[previous_coma_index + 1:-1]))

                    element_marginal_array = np.ndarray(shape=dimensions)

                if strip_line[0] == "#":
                    if first_dim_line:
                        dimensions_names = []
                        if len(dimensions) > 1:
                            comma_index = strip_line.find(",")
                            opening_bracket_index = strip_line.find("[")
                            while opening_bracket_index != -1:
                                dimensions_names.append(
                                    strip_line[
                                        opening_bracket_index + 1:comma_index])
                                opening_bracket_index = strip_line.find(
                                    "[", comma_index)
                                comma_index = strip_line.find(
                                    ",", opening_bracket_index)
                        first_dim_line = False
                        dimensions_names.append(element_name)
                        self.network_dict[element_name] = dimensions_names

                    # update indices
                    indices_array = []
                    if len(dimensions) > 1:
                        comma_index = strip_line.find(",")
                        closing_brack_index = strip_line.find("]")
                        while closing_brack_index != -1:
                            indices_array.append(int(
                                strip_line[comma_index + 1:closing_brack_index]))
                            opening_bracket_index = strip_line.find(
                                "[", closing_brack_index)
                            comma_index = strip_line.find(
                                ",", opening_bracket_index)
                            closing_brack_index = strip_line.find(
                                "]", closing_brack_index + 1)

                if strip_line[0] == "%":
                    # read doubles
                    coma_index = strip_line.find(",")
                    marginals_values = []

                    # Get rid of the %
                    previous_coma_index = 0
                    while coma_index != -1:
                        marginals_values.append(
                            float(strip_line[previous_coma_index + 1:coma_index]))
                        previous_coma_index = coma_index
                        coma_index = strip_line.find(",", coma_index + 1)

                    # Add last dimension and get rid of the closing bracket
                    marginals_values.append(
                        float(strip_line[previous_coma_index + 1:]))
                    if len(marginals_values) != dimensions[-1]:
                        print("problem")
                    element_marginal_array[tuple(indices_array)] = marginals_values
            self.marginals_dict[element_name] = element_marginal_array

        #return marginals_dict, network_dict



    def write_model_marginals_from_parms(self, filename, model_parms_file=None, model_marginals_file=None):
        self.marginals_dict = {}
        self.network_dict = {}


#class Event:
#    """Class to manage Event_list in IGoR
#    """
#    def __init__(self):
#        # #GeneChoice;V_gene;Undefined_side;7;v_choice
#        # or
#        # #Deletion;J_gene;Five_prime;5;j_5_del
#        self.eventType     = ""
#        self.targetGene    = ""
#        self.geneSide      = ""
#        self.eventPriority = ""
#        self.eventNickname = ""
#
#        # %U66059|TRBV9*01|Homo sapiens|F|V-REGION|206836..207121|286 nt|1| | | | |286+0=286| | |;GATTCTGGAGTCACACAAACCCCAAAGCACCTGATCACAGCAACTGGACAGCGAGTGA    CGCTGAGATGCTCCCCTAGGTCTGGAGACCTCTCTGTGTACTGGTACCAACAGAGCCTGGACCAGGGCCTCCAGTTCCTCATTCAGTATTATAATGGAGAAGAGAGAGCAAAAGGAAACATTCTTGAACGATTCTCCGCACAACAG    TTCCCTGACTTGCACTCTGAACTAAACCTGAGCTCTCTGGAGCTGGGGGACTCAGCTTTGTATTTCTGTGCCAGCAGCGTAG;88
#        # %18;22
#        self.realization   = ""


#
#    def read_model_parms(self, filename):
#        """Reads a model graph structure from a model params file.
#        Note that for now this method does not read the error rate information.
#
#        """
#        with open(filename, "r") as file:
#            # dictionary containing recarrays?
#            value_str = []
#            index = []
#            real_name = []
#            current_part = ""
#
#            line = file.readline()
#            strip_line = line.rstrip('\n')  # Remove end of line character
#            strip_line = strip_line.rstrip(
#                '\r')  # Remove carriage return character (if needed)
#
#            if strip_line == "@Event_list":
#                line = file.readline()
#                strip_line = line.rstrip('\n')  # Remove end of line character
#                strip_line = strip_line.rstrip(
#                    '\r')  # Remove carriage return character (if needed)
#
#                while strip_line[0] == '#':
#                    semicolon_index = strip_line.find(';')
#                    event_type = strip_line[1:semicolon_index]
#
#                    next_semicolon_index = strip_line.find(';',
#                                                           semicolon_index + 1)
#                    seq_type = strip_line[
#                               semicolon_index + 1:next_semicolon_index]
#
#                    semicolon_index = next_semicolon_index
#                    next_semicolon_index = strip_line.find(';',
#                                                           semicolon_index + 1)
#                    seq_side = strip_line[
#                               semicolon_index + 1:next_semicolon_index]
#
#                    semicolon_index = next_semicolon_index
#                    next_semicolon_index = strip_line.find(';',
#                                                           semicolon_index + 1)
#                    priority = strip_line[
#                               semicolon_index + 1:next_semicolon_index]
#
#                    semicolon_index = next_semicolon_index
#                    next_semicolon_index = strip_line.find(';',
#                                                           semicolon_index + 1)
#                    nickname = strip_line[semicolon_index + 1:]
#
#                    event = Rec_Event(event_type, seq_type, seq_side, priority,
#                                      nickname)
#
#                    line = file.readline()
#                    strip_line = line.rstrip('\n')  # Remove end of line character
#                    strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
#                    while strip_line[0] == '%':
#                        name = ""
#                        semicolon_index = 0
#                        next_semicolon_index = strip_line.find(';')
#                        if event_type == "GeneChoice":
#                            name = strip_line[
#                                   semicolon_index + 1:next_semicolon_index]
#                            semicolon_index = next_semicolon_index
#                            next_semicolon_index = \
#                                strip_line.find(';', semicolon_index + 1)
#                            value = strip_line[
#                                    semicolon_index + 1:next_semicolon_index]
#                        elif event_type == "DinucMarkov":
#                            value = strip_line[
#                                    semicolon_index + 1:next_semicolon_index]
#                        else:
#                            value = int(strip_line[
#                                        semicolon_index + 1:next_semicolon_index])
#
#                        semicolon_index = next_semicolon_index
#                        next_semicolon_index = strip_line\
#                            .find(';', semicolon_index + 1)
#
#                        index = int(strip_line[semicolon_index + 1:])
#
#                        realization = Event_realization(name, value, index)
#
#                        event.add_realization(realization)
#
#                        line = file.readline()
#                        strip_line = line.rstrip('\n')  # Remove end of line character
#                        strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
#                    self.events.append(event)
#
#            if strip_line == "@Edges":
#
#                line = file.readline()
#                strip_line = line.rstrip('\n')  # Remove end of line character
#                strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
#                while strip_line[0] == '%':
#                    semicolon_index = strip_line.find(';')
#                    parent = strip_line[1:semicolon_index]
#                    child = strip_line[semicolon_index + 1:]
#
#                    if not (parent in self.edges):
#                        self.edges[parent] = Adjacency_List()
#
#                    if not (child in self.edges):
#                        self.edges[child] = Adjacency_List()
#
#                    self.edges[parent].children.append(child)
#                    self.edges[child].parents.append(parent)
#
#                    line = file.readline()
#                    strip_line = line.rstrip('\n')  # Remove end of line character
#                    strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
#
#            # FIXME: ErrorRate added
#            if strip_line == "@ErrorRate" :
#                line = file.readline()
#                strip_line = line.rstrip('\n')  # Remove end of line character
#                strip_line = strip_line.rstrip('\r')  # Remove carriage return character (if needed)
#                while strip_line[0] == '#':
#                    if 'SingleErrorRate' == strip_line[1:] :
#                        line = file.readline()
#                        strip_line = line.rstrip('\n').rstrip()
#                        error = strip_line
#                        self.error_rate = {"SingleErrorRate" : error }
#
#    def get_event(self, event_name, by_nickname=False):
#        """Returns the RecEvent with corresponding name or nickname."""
#        if by_nickname:
#            for ev in self.events:
#                if ev.nickname == event_name:
#                    return ev
#            raise Exception(
#                'RecEvent with nickname \"' + event_name + "\" not found.")
#        else:
#            for ev in self.events:
#                if ev.name == event_name:
#                    return ev
#            raise Exception(
#                'RecEvent with name \"' + event_name + "\" not found.")
#
#
#class Rec_Event:
#    """Recombination event class containing event's name, type, realizations,
#    etc... Similar to IGoR's C++ RecEvent class.
#
#    """
#    def __init__(self, event_type, seq_type, seq_side, priority,
#                 nickname=None):
#        self.event_type = event_type
#        self.seq_type = seq_type
#        self.seq_side = seq_side
#        self.priority = priority
#        self.realizations = list()
#        self.name = ""
#        self.nickname = ""
#        if nickname is not None:
#            self.nickname = nickname
#        self.update_name()
#
#    def __str__(self):
#        return self.name
#
#    def __repr__(self):
#        return "Rec_event(" + self.name + ")"
#
#    def add_realization(self, realization):
#        """Add a realization to the RecEvent realizations list."""
#        self.realizations.append(realization)
#        self.update_name()
#
#    def update_name(self):
#        """Updates the name of the event (will have no effect if the RecEvent
#        has not been modified since the last call).
#
#        """
#        if self.event_type == "DinucMarkov":
#            self.name = self.event_type + "_" + self.seq_type + "_" + \
#                        self.seq_side + "_prio" + \
#                        str(self.priority) + "_size" + \
#                        str(len(self.realizations) ** 2)
#        else:
#            self.name = self.event_type + "_" + self.seq_type + "_" + \
#                        self.seq_side + "_prio" + \
#                        str(self.priority) + "_size" + \
#                        str(len(self.realizations))
#
#    def get_realization_vector(self):
#        """This methods returns the event realizations sorted by the
#        realization index as a list.
#
#        """
#        if self.event_type == 'GeneChoice':
#            tmp = [""] * len(self.realizations)  # empty(, dtype = str)
#        else:
#            tmp = numpy.empty(len(self.realizations),
#                              dtype=type(self.realizations[0].value))
#        # print("Unfinished method get realization vector")
#        processed_real_indices = []
#        for real in self.realizations:
#            if processed_real_indices.count(real.index) == 0:
#                if real.name != "":
#                    tmp[real.index] = real.name
#                else:
#                    tmp[real.index] = real.value
#                processed_real_indices.append(real.index)
#            else:
#                print("REALIZATION INDICES ARE DEGENERATE")
#
#        return tmp
#
#
#class Adjacency_List:
#    """Graph adjacency list. Similar to IGoR's C++ Adjacency_List class. Stores
#    adjacent events (parents and children) to a given event.
#
#    """
#    def __init__(self):
#        self.children = list()
#        self.parents = list()
#
#
#class Event_realization:
#    """A small class storing for each RecEvent realization its name, value and
#    corresponding index.
#
#    """
#    def __init__(self, name, value, index):
#        self.name = name
#        self.value = value
#        self.index = index
#
#    def __lt__(self, other):
#        return self.index < other.index
#
#    def __gt__(self, other):
#        return self.index > other.index
#
#
#def compute_average_distribution(event_name, model, averaging_list=None):
#    """Returns the 1D distribution averaged over the event ancestors.
#    Alternatively a list of ancestors can be provided, such that the
#    distribution is only averaged over these ancestors.
#
#    """
#    for event in model.events:
#        if event.name == event_name:
#            event_cluster = list()
#            event_cluster.append(event)
#
#            # Get all events related (even remotely) to this considered event
#            explored = False
#            parents_names = list()
#            parents_names.append(event.name)
#            while not explored:
#                next_parents = list()
#                for parent in parents_names:
#                    if parent not in model.edges:
#                        continue
#                    parent_parents = model.edges[parent].parents
#                    for p in parent_parents:
#                        next_parents.append(p)
#                        if event_cluster.count(model.get_event(p)) == 0:
#                            event_cluster.append(model.get_event(p))
#                parents_names = next_parents
#                if len(parents_names) == 0:
#                    explored = True
#
#            dimension_list = []
#            dimension_names_list = []
#            for ev in event_cluster:
#                dimension_list.append(
#                    model.marginals[0][ev.nickname].shape[-1])
#                dimension_names_list.append(ev.nickname)
#
#            if event.event_type != "DinucMarkov":
#
#                final_array = numpy.ones(tuple(dimension_list))
#                # print(final_array.shape)
#
#                # Reshape all marginals in one big redundant array
#                reshaped_marginals = list()
#
#                # Takes care of model 1(non log part)
#                for ev in event_cluster:
#                    if ev.name not in model.edges:
#                        event_parents = []
#                    else:
#                        event_parents = model.edges[ev.name].parents
#                    # print(event_parents)
#                    dimension_list = []
#                    for e in event_cluster:
#                        if (event_parents.count(e.name) > 0) \
#                                or (ev.name == e.name):
#                            dimension_list.append(
#                                model.marginals[0][e.nickname].shape[-1])  # same as getting the number of realizations
#                        else:
#                            dimension_list.append(1)
#
#                    # print(ev)
#                    # print(dimension_list)
#                    # print(shape(model1.marginals[ev.nickname]))
#                    # print("Dimension list")
#                    # print(dimension_list)
#
#                    swapped_marginal = copy.deepcopy(
#                        model.marginals[0][ev.nickname])
#                    swapped_dimension_names = copy.deepcopy(
#                        model.marginals[1][ev.nickname])
#
#                    is_sorted = False
#                    while not is_sorted:
#                        for i in range(0, len(swapped_marginal.shape)):
#                            if i == len(swapped_marginal.shape) - 1:
#                                is_sorted = True
#                                break
#                            elif (find(numpy.asarray(dimension_names_list) ==
#                                       swapped_dimension_names[i]) > find(
#                                    numpy.asarray(dimension_names_list) ==
#                                    swapped_dimension_names[i + 1])):
#                                swapped_marginal = swapped_marginal.swapaxes(
#                                    i, i + 1)
#                                # Artisanal swap
#                                tmp = swapped_dimension_names[i + 1]
#                                swapped_dimension_names[i + 1] = \
#                                    swapped_dimension_names[i]
#                                swapped_dimension_names[i] = tmp
#                                break
#
#                    reshaped_marginals.append(
#                        swapped_marginal.reshape(tuple(dimension_list)))
#
#                for reshaped_marg in reshaped_marginals:
#                    final_array *= reshaped_marg
#
#                # if no list given average over all dependencies
#                sum_axes = range(0, len(reshaped_marginals[0].shape))
#                sum_axes.remove(find(
#                    numpy.asarray(dimension_names_list) == event.nickname))
#                if averaging_list is not None:
#                    # for dependance_name in averaging_list:
#                    # sum_axes.remove(find(numpy.asarray(dimension_names_list)==dependance_name))
#                    for dimension_name in dimension_names_list:
#                        if dimension_name != event.nickname and all(
#                                numpy.asarray(
#                                        averaging_list) != dimension_name):
#                            sum_axes.remove(find(numpy.asarray(
#                                dimension_names_list) == dimension_name))
#                return final_array.sum(axis=tuple(sum_axes))
#
#
## class Model_Marginals:
## Model marginals are for now quite uneasy to study, should be turned into a
## handy class for each event a small object containing the array and another a
## tuple with the dimension names ordered.
#def read_marginals_txt(filename, dim_names=False):
#    """Reads a model marginals file. Returns a tuple containing a dict
#    containing the individual events probabilities indexed by the events
#    nicknames and a dict containing the list of dimension names/ordering for
#    each event.
#
#    """
#    with open(filename, "r") as f:
#        # Model parameters are stored inside a dictionnary of ndarrays
#        model_dict = {}
#        network_dict = {}
#        element_name = ""
#        first = True
#        first_dim_line = False
#        element_marginal_array = []
#        indices_array = []
#
#        for line in f:
#            strip_line = line.rstrip("\n")  # Remove end of line character
#            if strip_line[0] == "@":
#                first_dim_line = True
#                if not first:
#                    # Add the previous to the dictionnary
#                    model_dict[element_name] = element_marginal_array
#                else:
#                    first = False
#
#                element_name = strip_line[1:]
#            # print element_name
#
#            if strip_line[0] == "$":
#                # define array dimensions
#                coma_index = strip_line.find(",")
#                dimensions = []
#
#                # Get rid of $Dim[
#                previous_coma_index = 4
#                while coma_index != -1:
#                    dimensions.append(
#                        int(strip_line[previous_coma_index + 1:coma_index]))
#                    previous_coma_index = coma_index
#                    coma_index = strip_line.find(",", coma_index + 1)
#
#                # Add last dimension and get rid of the closing bracket
#                dimensions.append(int(strip_line[previous_coma_index + 1:-1]))
#
#                element_marginal_array = numpy.ndarray(shape=dimensions)
#
#            if strip_line[0] == "#":
#                if first_dim_line:
#                    dimensions_names = []
#                    if len(dimensions) > 1:
#                        comma_index = strip_line.find(",")
#                        opening_bracket_index = strip_line.find("[")
#                        while opening_bracket_index != -1:
#                            dimensions_names.append(
#                                strip_line[
#                                    opening_bracket_index + 1:comma_index])
#                            opening_bracket_index = strip_line.find(
#                                "[", comma_index)
#                            comma_index = strip_line.find(
#                                ",", opening_bracket_index)
#                    first_dim_line = False
#                    dimensions_names.append(element_name)
#                    network_dict[element_name] = dimensions_names
#
#                # update indices
#                indices_array = []
#                if len(dimensions) > 1:
#                    comma_index = strip_line.find(",")
#                    closing_brack_index = strip_line.find("]")
#                    while closing_brack_index != -1:
#                        indices_array.append(int(
#                            strip_line[comma_index + 1:closing_brack_index]))
#                        opening_bracket_index = strip_line.find(
#                            "[", closing_brack_index)
#                        comma_index = strip_line.find(
#                            ",", opening_bracket_index)
#                        closing_brack_index = strip_line.find(
#                            "]", closing_brack_index + 1)
#
#            if strip_line[0] == "%":
#                # read doubles
#                coma_index = strip_line.find(",")
#                marginals_values = []
#
#                # Get rid of the %
#                previous_coma_index = 0
#                while coma_index != -1:
#                    marginals_values.append(
#                        float(strip_line[previous_coma_index + 1:coma_index]))
#                    previous_coma_index = coma_index
#                    coma_index = strip_line.find(",", coma_index + 1)
#
#                # Add last dimension and get rid of the closing bracket
#                marginals_values.append(
#                    float(strip_line[previous_coma_index + 1:]))
#                if len(marginals_values) != dimensions[-1]:
#                    print("problem")
#                element_marginal_array[tuple(indices_array)] = marginals_values
#        model_dict[element_name] = element_marginal_array
#
#    return model_dict, network_dict
