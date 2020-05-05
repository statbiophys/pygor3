#!/usr/bin/env python3
import pygor3 as p3


from  optparse import OptionParser
def main():
    #parser = OptionParser()
    #parser.add_option("-s", "--species", dest="species", help='Igor species')
    #parser.add_option("-c", "--chain", dest="chain", help='Igor chain')
    #parser.add_option("-b", "--batch", dest="batch", help='Batchname to identify run. If not set random name is generated')
    #parser.add_option("-o", "--output", dest="output", help='filename of csv file to export data')

    #(options, args) = parser.parse_args()
    #print("hola fasdasdfd 888888")
    species="human"
    chain="tcr_beta"
    mdl = p3.IgorModel.load_default(species, chain)
    # import hvplot
    # hvplot.save(mdl.parms.plot_Graph(), species+"__"+chain+".html")

    # print(mdl.xdata)
    print(mdl.parms.G.edges)
    # print(mdl.parms.G.nodes)
    for node in mdl.parms.G.nodes:
        print(node)

    # TODO: I want the marginals of each event, so for each event
    strEvent = 'v_choice'
    strEvent = 'j_choice'
    strEvent = 'd_gene'
    strEvent = 'v_3_del'
    strEvent = 'd_3_del'
    strEvent = 'd_5_del'
    strEvent = 'j_5_del'
    strEvent = 'vd_ins'
    strEvent = 'dj_ins'
    strEvent = 'vd_dinucl'
    strEvent = 'dj_dinucl'

    # strEvent = 'd_3_del'
    # mdl.get_events_types_list()
    for strEvent in ['d_gene'] : #, 'd_3_del', 'vd_ins', 'dj_dinucl']:
        print(strEvent)
        da = mdl[strEvent]
        pd = mdl.parms.Event_dict[strEvent] #['name'].map(genLabel).values  # FIXME: use the exact name defined in model_parms
        print(pd)
        print(da)

    print(mdl.parms.Event_dict['v_choice'])

    import hvplot
    import hvplot.pandas

    strEvent = 'j_choice'
    print("= " * 20)
    evento = mdl.parms.get_Event(strEvent)
    df = mdl.parms.Event_dict[strEvent]  #mdl.Pmarginal['d_gene'].to_dataframe(name='P')
    print(df)

    print("* " * 20)
    df['name'] = df['name'].apply(p3.genLabel)
    print(df)
    # new_df = df
    evento.update_realizations_from_dataframe(df)
    new_df = mdl.parms.Event_dict[strEvent]
    print(new_df)

    ####################
    print( mdl.Pmarginal['dj_dinucl'] )
    print(mdl.Pmarginal['dj_dinucl'].sum())

    print( mdl.Pmarginal['dj_dinucl'].dims )



    import matplotlib.pyplot as plt
    fig, ax = mdl.plot_Event_Marginal('dj_dinucl')
    plt.show()


    # fig, ax = plt.subplots()
    # da.plot(ax=ax, x='x', y='y')
    # # ax.set_xticks(ticks=mdl[strEvent]["x"].values, labels=mdl[strEvent]["lbl__x"].values)
    # ax.set_xticks(da['x'].values)
    # ax.set_xticklabels(da['lbl__' + 'x'].values)  #, rotation=90)
    # ax.set_yticks(da['y'].values)
    # ax.set_yticklabels(da['lbl__' + 'y'].values)
    # print(da)
    # plt.show()



    # print(df)
    # print( mdl.xdata['d_3_del'] )  # P( D3| D5,D)
    # print(mdl.xdata['d_5_del'])  # P( D3| D5,D)
    # da = mdl.xdata['d_3_del'] * mdl.xdata['d_5_del']
    # print("0000000000"*10)  # P( D3| D5,D)
    #
    # for ii in da[strEvent][strEvent].values:
    #     print("--"*5,ii,da[strEvent]["lbl__"+strEvent].values[ii], "-"*30)
    #     print(da[{strEvent:ii}])

    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots()
    # da[{strEvent: ii}].plot(ax = ax)
    # plt.show()






    """
    directory = "/home/alfaceor/Dropbox/PosDoc/IGoR/Thierry/Lightchain_Thierry/Lightchain/IGL/"
    flnVanchors = directory + "V_gene_CDR3_anchors.csv"
    flnJanchors = directory + "J_gene_CDR3_anchors.csv"
    fln_model_parms = directory + "final_parms.txt"
    fln_model_marginals = directory + "final_marginals.txt"
    mdl = p3.IgorModel(model_parms_file=fln_model_parms, model_marginals_file=fln_model_marginals)
    mdl_anchs = p3.IgorAnchors(flnVanchors, flnJanchors)
    mdl_parms = p3.IgorModel_Parms(model_parms_file=fln_model_parms)
    mdl_marginals = p3.IgorModel_Marginals(model_marginals_file=fln_model_marginals)
    print(mdl_marginals)
    evento = mdl_parms.get_Event('v_choice')  # get_realization_DataFrame()
    print(evento)
    ev_da = evento.get_realization_DataFrame()
    ev_da.head()
    # remove the ones without
    ev_da_drop = ev_da.drop([0, 2, 4])  # .head()
    aaa = ev_da_drop.head()

    import pandas as pd
    df_Vanchors = mdl_anchs.df_Vanchors.rename(columns={'gene': 'name'})
    result = pd.merge(ev_da, df_Vanchors, on='name')  # , left_index=True, right_index=True, how='inner');
    new_df = result.drop(columns=['anchor_index'])
    # set dataframe to envento

    new_df
    """

if __name__ == "__main__":
    main()
