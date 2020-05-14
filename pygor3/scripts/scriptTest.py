#!/usr/bin/env python3
import pygor3 as p3
import numpy as np

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

    igor_genome = p3.IgorGenomics.load_from_path_ref_genome(".")
    df_joined00 = igor_genome.df_V_gene.set_index('name').join(igor_genome.df_V_anchors.set_index('gene'))
    df_joined00 = df_joined00.reset_index()
    print(df_joined00)
    print(type(df_joined00.loc[3]['anchor_index']))
    print(df_joined00.loc[ np.isnan(df_joined00["anchor_index"] ) ] )
    df = df_joined00.loc[df_joined00['anchor_index'] == np.NaN]
    print(df)

    igor_genome.df_V_gene['ID'] = igor_genome.df_V_gene['name'].apply(p3.genLabel)
    igor_genome.df_V_anchors['ID'] = igor_genome.df_V_anchors['gene'].apply(p3.genLabel)


    # print(igor_genome.df_V_gene)
    # print(igor_genome.df_V_anchors)
    df_joined = igor_genome.df_V_gene.set_index('ID').join(igor_genome.df_V_anchors.set_index('ID')).reset_index()
    print(df_joined)
    print(df_joined.loc[np.isnan(df_joined["anchor_index"])])

    print(igor_genome.df_V_gene.loc[igor_genome.df_V_gene['ID']  == 'TRBV17*01'])
    # From dataframe create the two files genomic template and also CDR3 anchors

    # load data into dataframe and then export to fasta and csv. show the excluded



    # print(df_tmp)
    # for gene_id in df_tmp['name']:
    #     print(gene_id.split("*"))



    """
    mdl = p3.IgorModel.load_default(species, chain)
    # import hvplot
    # hvplot.save(mdl.parms.plot_Graph(), species+"__"+chain+".html")

    # genomicDs__imgt.fasta  genomicJs__imgt.fasta  genomicJs__imgt_trim.fasta  genomicVs__imgt.fasta  genomicVs__imgt_trim.fasta

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
    """


if __name__ == "__main__":
    main()
