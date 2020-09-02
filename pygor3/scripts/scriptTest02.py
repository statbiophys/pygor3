#!/usr/bin/env python3
import pygor3 as p3
import argparse

def get_CDR3_data(mdl, scen):
    df_V = mdl.genomic_dataframe_dict['V']
    mdl.parms.get_event_dict()
    print ( str(df_V[df_V['id'] == scen['v_choice']]['gene_name'].values) )

def main():
    parser = argparse.ArgumentParser()
    igor_models = parser.add_argument_group('IGoR default models')

    parser.add_argument("-D", "--database", dest="database", help="Igor database created with database script.",
                        default="Ajam.db")

    args = parser.parse_args()


    db = p3.IgorSqliteDB.create_db(args.database)
    # db.load_Igor
    mdl = db.get_IgorModel()
    # Prob = mdl.xdata['v_choice']*mdl.xdata['j_choice']
    # print(Prob)

    print("="*30)
    #records = db.fetch_IgorGeneAnchors_By_Gene("V")
    #print(records)

    mdl.set_genomic_dataframe_dict(db.get_IgorGenomicDataFrame_dict())

    # Now getting the CDR3 from best scenario
    print( mdl.genomic_dataframe_dict["V"] )

    # Get the CDR3 length from a bs
    seq_index = 5
    scenarios_list = db.get_IgorBestScenarios_By_seq_index(seq_index)
    print(scenarios_list[0].to_dict())

    scenarios_list = db.get_IgorBestScenarios_By_seq_index_IgorModel(seq_index, mdl)

    # Ahora que tengo los scenarios en el formato del mdl
    # quiero obtener el cdr3 de un scenario
    scen = scenarios_list[0]
    print(scen.to_dict())
    print(mdl.parms.Event_list[0].to_dict())
    # mdl.parms.get_event_dict('')

    #lista = db.get_IgorAlignment_data_list_By_seq_index("V", seq_index)
    v_aln_list = db.get_IgorAlignment_data_list_By_seq_index("V", seq_index)
    j_aln_list = db.get_IgorAlignment_data_list_By_seq_index("J", seq_index)
    # str_gene_name = "U66059|TRBV5-1*01|Homo sapiens|F|V-REGION|113806..114091|286 nt|1| | | | |286+0=286| | |"
    gene_templates_dict = db.get_IgorGenomicDataFrame_dict()

    ### JUST ALIGNMENTS
    for aln in v_aln_list:
        # print(gene_templates_dict["V"])
        anchor_index = gene_templates_dict["V"].loc[aln.gene_id].anchor_index
        aln.anchor_in_read = anchor_index + aln.offset
        print(" v: ", aln.gene_id, aln.score, aln.anchor_in_read, anchor_index, aln.offset, aln.strGene_name)

    # print(gene_templates_dict["J"])
    for aln in j_aln_list:
        # print(gene_templates_dict["V"])
        #print(aln.gene_id,)
        anchor_index = gene_templates_dict["J"].loc[aln.gene_id].anchor_index
        aln.anchor_in_read = anchor_index + aln.offset
        print(" j: ", aln.gene_id, aln.score, aln.anchor_in_read, anchor_index, aln.offset, aln.strGene_name)

    print(v_aln_list[0].anchor_in_read)
    print(j_aln_list[0].anchor_in_read)

    # NOW USING THE SCENARIOS
    scenarios_list = db.get_IgorBestScenarios_By_seq_index_IgorModel(seq_index, mdl)
    for scen in scenarios_list:
        # scen = scenarios_list[0]
        print("---> Scenario", scen.scenario_rank)
        print('v_choice', scen['v_choice'])
        print('j_choice', scen['j_choice'])
        # FIXME : I'm assuming the same index for parms and gene templates
        v_anchor_index = gene_templates_dict["V"].loc[scen['v_choice']].anchor_index
        j_anchor_index = gene_templates_dict["J"].loc[scen['j_choice']].anchor_index
        # Now use the alignments to get the list of the alignments with 'v_choice' id a seq_index
        #print(gene_templates_dict["V"].loc[0].anchor_index)
        scen_v_alns = [v_aln.anchor_in_read for v_aln in v_aln_list if v_aln.gene_id == scen['v_choice']]
        scen_j_alns = [j_aln.anchor_in_read for j_aln in j_aln_list if j_aln.gene_id == scen['j_choice']]
        print(scen_v_alns)
        print(scen_j_alns)


    # print(gene_templates_dict["V"])


    # print(list( map(lambda x: (x.strGene_name, x.score, x.offset, x.anchor_in_read), lista)) )

    return 0

    # str_gene_name = "TRBD1*01"
    # print(lista.to_dict())
    lista02 = [aln for aln in lista if aln.strGene_name == str_gene_name]
    print(lista02[0].to_dict())
    # print(mdl.parms.Event_dict['v_choice'] ) #.iloc[0])

    # Make a correspondace table between model gene id and alignment gene id
    # IgorVGeneTemplate -> vgene_id
    # IgorER_v_choice -> id
    """
    SELECT t1.*, t2.* FROM IgorER_v_choice t1 INNER JOIN IgorVGeneTemplate t2 WHERE t1.name==t2.gene_name; 
    """

    # HOW TO GET THE ALIGNMENTS
    """
    SELECT t2.gene_name, t1.*, t2.sequence 
    FROM Igor{upper}Alignments t1  
    INNER JOIN Igor{upper}GeneTemplate t2 
    WHERE t1.seq_index==0 
    AND t1.{lower}gene_id == t2.{lower}gene_id 
    AND t2.gene_name=='TRBD1*01' 
    ORDER BY score DESC
    """

    """
    SELECT t2.gene_name, t1.*, t2.sequence 
    FROM IgorDAlignments t1  INNER JOIN IgorDGeneTemplate t2 
    WHERE t1.seq_index==0 AND t1.dgene_id == t2.dgene_id 
    ORDER BY score DESC
    """

    # TODO: IDEAL
    # scen['v_choice'].get_value(mdl)
    # scen['v_choice'].get_sequence(mdl)
    # scen.get_value('v_choice')
    # scen.get_name('v_choice')
    # scen.get_id('v_choice')


    # get_CDR3_data(mdl, scen)
    # scen = scenarios_list[0]
    # scen_dict = scen.realizations_ids_dict
    #
    # for event in mdl.parms.Event_list:
    #     if not (event.event_type == 'DinucMarkov'):
    #         scen_dict[event.nickname] = scen_dict.pop('id_' + event.nickname)


    # FIXME: DON'T KNOW WHY THE D GENES ARE
    #  A LITTLE BIT WEIRD.
    # strEvent = 'd_gene'
    # strEvent = 'd_3_del'
    # strEvent = 'd_5_del'
    # strEvent = 'vd_dinucl'
    # delta = mdl0.Pmarginal[strEvent] - mdl.Pmarginal[strEvent]
    # print("="*50)
    # print(mdl0.Pmarginal[strEvent])
    # print("-" * 50)
    # print(mdl.Pmarginal[strEvent])
    # print("-" * 50)
    # print(delta)

if __name__ == "__main__":
    main()
