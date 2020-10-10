#!/usr/bin/env python3
import pygor3 as p3
import argparse

def generate_event_delta_Kronecker(event_nickname, event_id):
    def tmp_funct(bs:p3.IgorScenario):
        if bs[event_nickname] == event_id:
            return 1
        else:
            return 0
    return tmp_funct

def get_pairwise_prob(mdl, event_nickname1, event_nickname2):
    # create an xarray matrix with event_nickname1 and event_nickname2
    # da = xr.
    da = mdl.get_zero_xarray_from_list([event_nickname1, event_nickname2])

    def tmp_funct(bs:p3.IgorScenario):
        ev_dict_ids = {event_nickname1: bs['id_'+event_nickname1], event_nickname2: bs['id_'+event_nickname2]}
        da.loc[ev_dict_ids] = 1
        return da

    return tmp_funct


def calc_average(db:p3.IgorSqliteDB, observable_function):
    for sigma in db.fetch_IgorIndexedSeq_indexes():
        print(sigma, len(db.fetch_IgorIndexedSeq_indexes()))

# from  optparse import OptionParser
def main():
    #parser = OptionParser()
    #parser.add_option("-s", "--species", dest="species", help='Igor species')
    #parser.add_option("-c", "--chain", dest="chain", help='Igor chain')
    #parser.add_option("-b", "--batch", dest="batch", help='Batchname to identify run. If not set random name is generated')
    #parser.add_option("-o", "--output", dest="output", help='filename of csv file to export data')

    parser = argparse.ArgumentParser()
    igor_models = parser.add_argument_group('IGoR default models')

    # # Use IGoR default model
    # igor_models.add_argument("-s", "--species", dest="species", help='Igor species', default="human")
    # igor_models.add_argument("-c", "--chain", dest="chain", help='Igor chain', default="tcr_beta")  # , type=str, choices=['TRB', 'TRA'])

    parser.add_argument("-D", "--database", dest="database", help="Igor database created with database script.", default="Ajam.db")
    parser.add_argument("--event_pair", dest="event_pair", help="Events nickname", nargs=2, default=['v_choice', 'j_choice'])

    igor_models = parser.add_argument_group('IGoR envent_pairdefault models')
    mdl = p3.IgorModel.load_default("mouse", "tcr_beta")
    args = parser.parse_args()

    fln_database = args.database
    db = p3.IgorSqliteDB.create_db(fln_database)
    records = db.get_list_of_tables_with_name('Igor%')
    print(records)
    db.delete_IgorModel_Tables()
    records = db.get_list_of_tables_with_name('Igor%')
    print(records)


    return 0

    # Which variable I want to eliminate
    strEvent = 'v_3_del'
    sorted_events = mdl.parms.get_Event_list_sorted()
    sorted_events_to_marginalize = [event for event in sorted_events if not event.event_type == "DinucMarkov"]
    sorted_events_to_marginalize_without_VE = [ event for event in sorted_events_to_marginalize if not event.nickname == strEvent]

    # Start eliminating events
    factors = mdl.VE_get_Pmarginals_initial_factors()
    for event_to_eliminate_VE in sorted_events_to_marginalize_without_VE:
        factors = mdl.VE_get_factors_by_sum_out_variable(event_to_eliminate_VE.nickname, factors)

    Pmarginal = 1
    # Now I want to print each factor
    for factor in factors:
        Pmarginal = Pmarginal * factor
    print(Pmarginal)

    # event_to_eliminate_VE = sorted_events_to_marginalize_without_VE.pop(0)
    # factors = mdl.VE_get_factors_by_sum_out_variable(event_to_eliminate_VE.nickname, factors)
    # print(len(factors))
    #
    # event_to_eliminate_VE = sorted_events_to_marginalize_without_VE.pop(0)
    # factors = mdl.VE_get_factors_by_sum_out_variable(event_to_eliminate_VE.nickname, factors)
    # print(len(factors))
    #
    # event_to_eliminate_VE = sorted_events_to_marginalize_without_VE.pop(0)
    # factors = mdl.VE_get_factors_by_sum_out_variable(event_to_eliminate_VE.nickname, factors)
    # print(len(factors))
    # # print(factors[-1])
    # # print( mdl.xdata['v_choice'] )
    # #
    # # print( mdl.xdata['v_choice']*factors[-1] )
    #
    # event_to_eliminate_VE = sorted_events_to_marginalize_without_VE.pop(0)
    # factors = mdl.VE_get_factors_by_sum_out_variable(event_to_eliminate_VE.nickname, factors)
    # print(len(factors))
    #
    # event_to_eliminate_VE = sorted_events_to_marginalize_without_VE.pop(0)
    # factors = mdl.VE_get_factors_by_sum_out_variable(event_to_eliminate_VE.nickname, factors)
    # print(len(factors))
    #
    # event_to_eliminate_VE = sorted_events_to_marginalize_without_VE.pop(0)
    # factors = mdl.VE_get_factors_by_sum_out_variable(event_to_eliminate_VE.nickname, factors)
    # print(len(factors))
    #
    # event_to_eliminate_VE = sorted_events_to_marginalize_without_VE.pop(0)
    # factors = mdl.VE_get_factors_by_sum_out_variable(event_to_eliminate_VE.nickname, factors)
    # print(len(factors))
    #
    # event_to_eliminate_VE = sorted_events_to_marginalize_without_VE.pop(0)
    # factors = mdl.VE_get_factors_by_sum_out_variable(event_to_eliminate_VE.nickname, factors)
    # print(len(factors))
    #
    # Pmarginal = 1
    # if len(sorted_events_to_marginalize_without_VE) == 0:
    #     # Now I want to print each factor
    #     for factor in factors:
    #         Pmarginal = Pmarginal*factor
    #     print(Pmarginal)
    #




    # print(mdl.Pmarginal['j_choice'])
    # da = mdl.xdata['j_choice']*mdl.xdata['v_choice']
    # print(da.sum('v_choice'))

    return 0
    sorted_events = mdl.parms.get_Event_list_sorted()
    sorted_events_to_marginalize = [event for event in sorted_events if not event.event_type == "DinucMarkov"]
    strEvent = 'j_choice'
    event_to_marginalize = mdl.parms.get_Event(strEvent)
    event_to_eliminate_ordered_list = [event for event in sorted_events_to_marginalize if not
                                       event.nickname == event_to_marginalize.nickname]
    factors = mdl.VE_get_Pmarginals_initial_factors()
    print("event_to_marginalize.nickname : ", event_to_marginalize.nickname)
    for event_to_eliminate in event_to_eliminate_ordered_list:
        print("event_to_eliminate : ", event_to_eliminate.nickname)
        mdl.VE_get_factors_by_sum_out_variable(event_to_eliminate.nickname)
        print(mdl.factors)

    # mdl.Pmarginal[event_to_marginalize.nickname] = mdl.factors
    # print("-------- FACTORS")
    # print(mdl.factors)
    # print("====================")
    # print(mdl.xdata[strEvent])
    # print("--------------------")
    # print( mdl.Pmarginal[strEvent] )
    # print("////////////////////")
    # # At the end only one factor should survive

    return 0

    args = parser.parse_args()
    print(args)
    print(args.event_pair)
    str_event_nickname1 = args.event_pair[0]
    str_event_nickname2 = args.event_pair[1]

    # Create an IgorModel
    #mdl = p3.IgorModel.load_default(args.species, args.chain)
    db = p3.IgorSqliteDB.create_db(args.database)
    mdl = db.get_IgorModel()
    print(mdl['j_choice'])
    # TODO: ahora que tenemos el nombre de
    #  las columnas debemos hacer que esto vaya a un xarray
    #  o a por lo menos un numpy array con los datos
    #  1. Opcion: puede ser pandas luego xarray no se como eso es posible
    # return 0

    fff = get_pairwise_prob(mdl, str_event_nickname1, str_event_nickname2)
    average = db.calc_IgorBestScenarios_average_of(fff)
    print(average)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(10,15))
    average.plot(ax=ax)

    str_title = "$P($" + str_event_nickname1 + ", " + str_event_nickname2+"$)$"
    #ax.set_title(str_title)
    ax.set_ylabel(str_event_nickname1.replace("_", " "))
    ax.set_xlabel(str_event_nickname2.replace("_", " "))

    # plot sutilties
    import numpy as np
    v_genlabel = np.vectorize(p3.genLabel )

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


if __name__ == "__main__":
    main()
