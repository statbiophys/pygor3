#!/usr/bin/env python3
import pygor3 as p3
import argparse

def get_pairwise_prob(mdl, event_nickname1, event_nickname2):
    # create an xarray matrix with event_nickname1 and event_nickname2
    da = mdl.get_zero_xarray_from_list([event_nickname1, event_nickname2])
    # print("da = ", da)
    import xarray as xr

    def tmp_funct(bs:p3.IgorScenario):
        ev_dict_ids = {event_nickname1: bs['id_'+event_nickname1], event_nickname2: bs['id_'+event_nickname2]}
        darr = xr.zeros_like(da)
        darr.loc[ev_dict_ids] = 1
        return darr

    return tmp_funct

def main():
    parser = argparse.ArgumentParser()
    igor_models = parser.add_argument_group('IGoR default models')

    # # Use IGoR default model
    # igor_models.add_argument("-s", "--species", dest="species", help='Igor species', default="human")
    # igor_models.add_argument("-c", "--chain", dest="chain", help='Igor chain', default="tcr_beta")  # , type=str, choices=['TRB', 'TRA'])

    parser.add_argument("-D", "--database", dest="database", help="Igor database created with database script.", default="Ajam.db")
    parser.add_argument("--event_pair", dest="event_pair", help="Events nickname", nargs=2, default=['v_choice', 'j_choice'])

    igor_models = parser.add_argument_group('IGoR envent_pairdefault models')

    args = parser.parse_args()

    print(args)
    print(args.event_pair)
    str_event_nickname1 = args.event_pair[0]
    str_event_nickname2 = args.event_pair[1]

    # Create an IgorModel
    # mdl = p3.IgorModel.load_default(args.species, args.chain)
    db = p3.IgorSqliteDB.create_db(args.database)
    mdl = db.get_IgorModel()
    func_pairwise = get_pairwise_prob(mdl, str_event_nickname1, str_event_nickname2)
    seq_index = 3
    average = db.calc_IgorBestScenarios_average_of(func_pairwise) #, indices_list=[seq_index])
    bs_list = db.get_IgorBestScenarios_By_seq_index_IgorModel(seq_index, mdl)
    bs = bs_list[0]
    bs_realiz = mdl.parms.realiz_dict_from_scenario(bs)
    print(bs_realiz )
    # print(bs_realiz['v_choice'].to_dict())
    print( bs_realiz['mismatches'].value )
    return 0

    print("=" * 50)
    bs_records = db.fetch_IgorBestScenarios_By_seq_index(seq_index)
    # import pandas as pd
    # bs_df = pd.DataFrame.from_records(bs_records)
    bs_df = db.get_IgorBestScenariosDataframe_By_seq_index(seq_index)

    print(bs_df[bs_df.id_v_choice == 59]['scenario_proba_cond_seq'].sum() )
    print(bs_df[bs_df.id_v_choice == 60]['scenario_proba_cond_seq'].sum())

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10, 15))
    aaa = average.plot(ax=ax, cmap='gnuplot2_r') #, vmin=0, vmax=1.0)
    print(type(aaa))

    str_title = "$P($" + str_event_nickname1 + ", " + str_event_nickname2+"$)$"

    ax.set_ylabel(str_event_nickname1.replace("_", " "))
    ax.set_xlabel(str_event_nickname2.replace("_", " "))

    import numpy as np
    v_genlabel = np.vectorize(p3.genLabel)

    XX = average[str_event_nickname2].values
    lbl_XX = average['lbl__' + str_event_nickname2].values

    YY = average[str_event_nickname1].values
    lbl_YY = average['lbl__' + str_event_nickname1].values

    ax.set_xticks(XX)
    ax.set_xticklabels(v_genlabel(lbl_XX), rotation=90)

    ax.set_yticks(YY)
    ax.set_yticklabels(v_genlabel(lbl_YY))

    fig.savefig('figura.pdf')
    # return ax

    plt.show()

    return 0

if __name__ == "__main__":
    main()
