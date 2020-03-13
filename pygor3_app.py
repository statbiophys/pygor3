#!/home/alfaceor/anaconda3/envs/pygor3/bin/python
import pygor3
import panel as pn
pn.extension()

Titulo = "<h1> Run IGoR<\h1>"
#igor_specie="human"
#igor_chain="tcr_beta"
igor_species_list = ["human", "mouse"]
    
igor_option_path_dict={
    "alpha": "tcr_alpha", 
    "beta" : "tcr_beta", 
    "light": "light", 
    "heavy_naive" : "bcr_heavy", 
    "heavy_memory": "bcr_heavy"
}

wd_igor_specie = pn.widgets.Select(options=igor_species_list, name=str_ev)
wd_igor_chain = pn.widgets.Select(options=igor_option_path_dict, name=str_ev)


"""Create a function that receives as input a IgorModel
and plot the dependecies of the events.
"""

#mdl0 = p3.IgorModel.load_default(igor_specie, igor_chain)
mdl0 = pygor3.IgorModel.load_default(igor_specie, igor_chain)


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
    print("pn_Events: ") #, pn_Events)
    for event_widget in event_widget_list:
        print(event_widget.value)
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
