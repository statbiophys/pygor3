#!/home/alfaceor/anaconda3/envs/pygor3/bin/python
"""Create a function that receives as input a IgorModel
and plot the dependecies of the events.
"""

# panel layout is a Row panel with
# pn_select_model: Select and submit button
# pn_model_events: Show with radiobutton the events of corresponding model and also show plot_Graph
# pn_event_plot: this should be a Row panel with selects and submit at left and plot at rigth

import panel as pn
pn.extension()
pn_description = pn.Row("Explore models")
pn_select_model = pn.Row()
pn_model_events = pn.Row()
pn_event_plot = pn.Row()
pn_layout = pn.Column(pn_description, pn_select_model, pn_model_events, pn_event_plot)

print(pn_layout)

import pygor3
mdl = pygor3.IgorModel()

igor_species_list = ["human", "mouse"]

igor_option_path_dict={
    "alpha": "tcr_alpha", 
    "beta" : "tcr_beta", 
    "light": "light", 
    "heavy_naive" : "bcr_heavy", 
    "heavy_memory": "bcr_heavy"
}

wd_specie = pn.widgets.Select(value="human", options=igor_species_list, name="specie")
wd_chain = pn.widgets.Select(value="tcr_beta", options=igor_option_path_dict, name="chain")
wd_mdl_button = pn.widgets.Button(name="submit", button_type='success')

pn_select_model = pn.Row(wd_specie, wd_chain, wd_mdl_button)

pn_layout[1] = pn_select_model

##############################################

def get_pn_model_events(mdl):
    button_event = pn.widgets.RadioButtonGroup(
        name='Event nicknames',
        options=[ ev.nickname for ev in mdl.parms.Event_list],
        button_type='default')
    return button_event

# change it to select buttono

def update_pn_model_events(event):
    mdl = pygor3.IgorModel.load_default(wd_specie.value, wd_chain.value)
    # Run function to update layout[2] = pn_model_events
    pn_layout[2] = get_pn_model_events(mdl) #(wd_specie.value, wd_chain.value)
    print(pn_layout[2])

wd_mdl_button.on_click(update_pn_model_events) #.param.watch(update_layout_events, 'value')


# Now for the next stages
# 1. checkout the chosen value in pn_layout[2]


pn_layout.show()

"""
def pn_depency_events(strEvent):
    # strEvent = 'd_gene'
    ##strEvent = button_event.value
    event_dependencies = list(mdl0.xdata[strEvent].dims)
    event_dependencies.remove(strEvent)
    event_dependencies
    # create a widget for each dependency
    event_widget_list = list()

    for str_ev in event_dependencies:
        labels = (mdl.xdata[strEvent]['lbl__' + str_ev].values).tolist()
        indexes = (mdl.xdata[strEvent][str_ev].values).tolist()
        event_options = dict(zip(labels, indexes))
        # print(options, type(options))
        event_widget_list.append(pn.widgets.Select(value=0, options=event_options, name=str_ev))

    # return event_widget_list
    pn_Events = pn.Column(*event_widget_list)
    print("pn_Events: ", pn_Events)
    return pn_Events






igor_specie="human"
igor_chain="tcr_beta"


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
    print("pn_Events: ", pn_Events)
    return pn_Events


#layout = pn.Column("Pygor Interface", button_event, update_depency_events) # pn_Events)
pn_plot = pn.Column("EL PLOT")
titulo = "<h1> IGoR model exploration </h1>"
pn_header = pn.Row(titulo, mdl0.parms.plot_Graph())
layout = pn.Column(pn_header, button_event, pn_depency_events(button_event.value), pn_plot) # pn_Events)
print(layout)


def plot_event(strEvent, chosen_event:dict):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    print(chosen_event)
    tipo_evento = mdl0.parms.get_Event(strEvent).event_type
    print(mdl0.xdata[strEvent][chosen_event])
    if tipo_evento == 'GeneChoice':
        print("strEvent : ", strEvent)
        XX = mdl0.xdata[strEvent][strEvent].values
        YY = mdl0.xdata[strEvent][chosen_event].values
        ax.bar(XX, YY)
        ax.set_xticks(mdl0.xdata[strEvent][strEvent].values)
        ax.set_xticklabels(mdl0.xdata[strEvent]['lbl__' + strEvent].values, rotation=90)
        # Plot as Genechoice
        #df = mdl0.xdata[strEvent][chosen_event].to_dataframe()
        #print(df)
        # df.plot.bar(x='lbl_'+strEvent, y='val', rot=90)
        #plot(ax=ax)
        #XX = mdl0.xdata[strEvent][chosen_event]
        #YY = mdl0.xdata[strEvent][chosen_event].values
        #mdl0.xdata[strEvent][chosen_event].plot(ax=ax)
        
    elif tipo_evento == 'Insertion':
        # Plot as Insertion
        print("strEvent : ", strEvent)
        mdl0.xdata[strEvent][chosen_event].plot(ax=ax)
    elif tipo_evento == 'Deletions':
        # Plot as Deletions
        print("strEvent : ", strEvent)
        mdl0.xdata[strEvent][chosen_event].plot(ax=ax)
    elif tipo_evento == 'DinucMarkov':
        # Plot as DinucMarkov
        print("strEvent : ", strEvent)
        mdl0.xdata[strEvent][chosen_event].plot(ax=ax)
    else:
        print("strEvent : ", strEvent)
        ax.plot()
        
    return fig
    


def update_layout_events(event):
    strEvent = button_event.value
    layout[2] = pn_depency_events(strEvent)
    #layout[3] = pn.Column()
    df = mdl0.parms.Event_dict[strEvent]
    #print(df)
    chosen_event = dict()
    for event_widget in layout[2]:
        event_widget.param.watch(update_layout_genes, 'value')

def update_layout_genes(event):
    chosen_event = dict()
    for event_widget in layout[2]:
        chosen_event[event_widget.name] = event_widget.value
    figura  = plot_event(button_event.value, chosen_event)
    layout[3] = pn.Column(figura)
    

button_event.param.watch(update_layout_events, 'value')

layout.show()

"""







#df = mdl0.parms.Event_dict[strEvent]


#event_widget.param.watch(update_plot, 'value')
#button_event.param.watch(update_depency_events, 'value')
#pn_Probs = pn.Row(pn_Events, "PLOT")
#layout.servable()
