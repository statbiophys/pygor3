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

def calc_average(db:p3.IgorSqliteDB, observable_function):
    for sigma in db.fetch_IgorIndexedSeq_indexes():
        print(sigma, len(db.fetch_IgorIndexedSeq_indexes()))


def observable():
    """
    Get a observable file or database, dataframe, etc ...
    return the result of the observable defined
    """
    obs_task = p3.IgorTask()
    print(obs_task.to_dict())
    obs_task._run_generate(20)
    obs_task.load_IgorModel()
    pass

def main():
    observable()


# from  optparse import OptionParser
def main_old():
    #parser = OptionParser()
    #parser.add_option("-s", "--species", dest="species", help='Igor species')
    #parser.add_option("-c", "--chain", dest="chain", help='Igor chain')
    #parser.add_option("-b", "--batch", dest="batch", help='Batchname to identify run. If not set random name is generated')
    #parser.add_option("-o", "--output", dest="output", help='filename of csv file to export data')

    parser = argparse.ArgumentParser()
    igor_models = parser.add_argument_group('IGoR default models')

    # Use IGoR default model
    igor_models.add_argument("-s", "--species", dest="species", help='Igor species')
    igor_models.add_argument("-c", "--chain", dest="chain", help='Igor chain')  # , type=str, choices=['TRB', 'TRA'])

    parser.add_argument("-D", "--database", dest="database", help="Igor database created with database script.")

    igor_models = parser.add_argument_group('IGoR default models')

    #(options, args) = parser.parse_args()
    # Create an IgorModel
    mdl = p3.IgorModel.load_default("human", "tcr_beta")
    da = mdl.xdata['v_choice']
    import numpy as np

    db = p3.IgorSqliteDB.create_db("Ajam.db")
    db.gen_IgorBestScenarios_cols_list()
    print(db.sql_IgorBestScenarios_cols_name_type_list)
    record = db.fetch_IgorBestScenarios_By_seq_index(0)[0]
    bs = p3.IgorScenario.load_FromSQLRecord(record, db.sql_IgorBestScenarios_cols_name_type_list)
    print(record)
    print(bs)
    print(bs.to_dict())
    func_V_delta = generate_event_delta_Kronecker('id_v_choice', 49)
    func_J_delta = generate_event_delta_Kronecker('id_j_choice', 0)
    # Now I have a kronecker delta function
    f = func_V_delta(bs)*func_J_delta(bs)
    print(f)

    func_V_delta = generate_event_delta_Kronecker('id_v_choice', 49)
    func_J_delta = generate_event_delta_Kronecker('id_j_choice', 0)

    db = p3.IgorSqliteDB.create_db("Ajam.db")
    db.gen_IgorBestScenarios_cols_list()
    indexes_list = db.fetch_IgorIndexedSeq_indexes()
    print("len: ", len(indexes_list))
    function_average = 0
    for indx in indexes_list:
        # aln_data = db.get_IgorAlignment_data_list_By_seq_index('V', indx)
        bs_data_list = db.fetch_IgorBestScenarios_By_seq_index(indx)
        tmp_average = 0
        scenario_norm_factor = 0
        for bs_data in bs_data_list:
            bs = p3.IgorScenario.load_FromSQLRecord(bs_data, db.sql_IgorBestScenarios_cols_name_type_list)
            tmp_average = tmp_average + bs.scenario_proba_cond_seq*func_V_delta(bs)*func_J_delta(bs)
            scenario_norm_factor = scenario_norm_factor + bs.scenario_proba_cond_seq
        tmp_average = tmp_average/scenario_norm_factor
        function_average = function_average + tmp_average

    print( function_average/len(indexes_list) )
    # Insert normalized probability in new table
    return 0




    # Make a query with the variables used for the specification
    # TODO: CORRECT FUNCTION OF CUMULATE PROBABILITY
    observable_vars_value = {
        'd_gene': " TRBD2*01", # id is 1
        'd_5_del': 4 # id is 8
    }
    # FIRST TRY TO WITH NAME IF NAME IS EMPTY THEN change it to value
    # TODO: return dict with the corresponding ids.
    # observable_ids = mdl.get_ids_dict_from_values(observable_vars_value)
    observable_ids = {
        'id_d_gene': [1],  # id is 1
        'id_d_5_del': [4]  # id is 8
    }

    records = db.fetch_IgorBestScenarios_By_events_dict(observable_ids)
    print( db.sql_cols_list )
    # For each record use a mdl to
    for record in records:
        print(record)

    print("*"*50)
    print(db.sql_cols_list)
    "SELECT name, type FROM pragma_table_info('IgorBestScenarios');"
    # ppp = db.execute_select_query("PRAGMA table_info(IgorBestScenarios);")
    ppp = db.execute_select_query("SELECT name, type FROM pragma_table_info('IgorBestScenarios');")
    db.sql_cols_list = db.execute_select_query("SELECT name, type FROM pragma_table_info('IgorBestScenarios');")

    # FIXME: but this is with database ids I want it with mdl ids,
    #  they are the same but is not the same 'v_choice' or 'id_v_choice'
    db.gen_IgorBestScenarios_cols_list()
    print(db.sql_IgorBestScenarios_cols_name_type_list)
    sce = p3.IgorScenario.load_FromSQLRecord(record, db.sql_cols_list)

    for ii, (col_name, tipo) in enumerate(db.sql_cols_list):
        if col_name == 'seq_index' :
            sce.seq_index = int(record[ii])
        elif col_name == 'scenario_rank' :
            sce.scenario_rank = int(record[ii])
        elif col_name == 'scenario_proba_cond_seq' :
            sce.scenario_proba_cond_seq = float(record[ii])
        else :
            if tipo == 'integer':
                sce.realizations_ids_dict[col_name] = int(record[ii])
            else:
                sce.realizations_ids_dict[col_name] = eval(record[ii])
        print(ii, col_name, tipo, record[ii])

if __name__ == "__main__":
    main()
