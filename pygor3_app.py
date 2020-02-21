#!/home/alfaceor/anaconda3/envs/pygor3/bin/python
import importlib
import pygor3
#importlib.reload(p3)
importlib.reload(pygor3)

igor_specie="human"
igor_chain="tcr_beta"

"""Create a function that receives as input a IgorModel
and plot the dependecies of the events.
"""

#mdl0 = p3.IgorModel.load_default(igor_specie, igor_chain)
mdl0 = pygor3.IgorModel.load_default(igor_specie, igor_chain)


import panel as pn
pn.extension()

button_event = pn.widgets.RadioButtonGroup(
    name='Event nicknames',
    options=[ ev.nickname for ev in mdl0.parms.Event_list],
    button_type='success')

#button_event


def pn_depency_events(strEvent):
    #strEvent = 'd_gene'
    ##strEvent = button_event.value
    event_dependencies = list( mdl0.xdata[strEvent].dims )
    event_dependencies.remove(strEvent)
    event_dependencies
    # create a widget for each dependency
    event_widget_list = list()
    
    for str_ev in event_dependencies:
        labels = (mdl0.xdata[strEvent]['lbl__'+str_ev].values).tolist()
        indexes = (mdl0.xdata[strEvent][str_ev].values).tolist()
        event_options=dict(zip(labels, indexes))
        #print(options, type(options))
        event_widget_list.append( pn.widgets.Select(value=0, options=event_options, name=str_ev) )
    
    #return event_widget_list
    pn_Events = pn.Column(*event_widget_list)
    print("pn_Events: ", pn_Events)
    return pn_Events

    
#layout = pn.Column("Pygor Interface", button_event, update_depency_events) # pn_Events)
layout = pn.Column("Pygor Interface", button_event, pn_depency_events(button_event.value)) # pn_Events)
print(layout)

def update_layout(event):
    layout[2] = pn_depency_events(button_event.value)

button_event.param.watch(update_layout, 'value')
#button_event.param.watch(update_depency_events, 'value')
#pn_Probs = pn.Row(pn_Events, "PLOT")
#layout.servable()
layout.show()
