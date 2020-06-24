#!/usr/bin/env python3
import pygor3 as p3

def Model_Marginals_CreateTable_SQL(event_nickname, list_dependencies:list):
    # from the xarray get the list of events something like:
    lista = [event_nickname]+list_dependencies
    # TODO:

    str_column_ct = "{}_id integer NOT NULL" #.format(event_nickname)
    str_foreign_key_ct = "FOREIGN KEY ({}_id) REFERENCES ER_{} (id)" #.format(event_nickname)

    sqlcmd_table_fields_ct = ",\n".join([str_column_ct.format(evento_nickname) for evento_nickname in lista])
    sqlcmd_foreign_keys_ct = ",\n".join([str_foreign_key_ct.format(evento_nickname, evento_nickname) for evento_nickname in lista])

    sqlcmd_ct_aux = """
            -- MM_XXXXXX table
            CREATE TABLE IF NOT EXISTS MM_{} (
                -- Events id columns
                {},
                P real,
                -- Foreign keys
                {}
            );
            """
    sqlcmd_ct = sqlcmd_ct_aux.format(event_nickname, sqlcmd_table_fields_ct, sqlcmd_foreign_keys_ct)
    return sqlcmd_ct


from  optparse import OptionParser
def main():
    #parser = OptionParser()
    #parser.add_option("-s", "--species", dest="species", help='Igor species')
    #parser.add_option("-c", "--chain", dest="chain", help='Igor chain')
    #parser.add_option("-b", "--batch", dest="batch", help='Batchname to identify run. If not set random name is generated')
    #parser.add_option("-o", "--output", dest="output", help='filename of csv file to export data')

    #(options, args) = parser.parse_args()
    #print("hola fasdasdfd 888888")

    # Create an IgorModel
    mdl_parms = p3.IgorModel_Parms()
    # mdl_parms.read_Event_list()

    # p3.IgorRec_Event(event_type, seq_type, seq_side, priority, nickname)
    print( mdl_parms.Event_list )

    mdl = p3.IgorModel.load_default("human", "tcr_beta")
    # mdl.export_event_to_csv('v_choice', 'putaque')

    db = p3.IgorSqliteDB.create_db("testando.db")
    db.load_IgorModel_Parms(mdl.parms) # TODO: finish edges and ErrorRate
    db.load_IgorModel_Marginals(mdl.xdata)

    # print(mdl.marginals.marginals_dict)

    """
    strEvent = 'vd_dinucl'
    print( mdl.xdata[strEvent] )
    print('='*50)
    da = mdl.xdata[strEvent]
    lista = list(dict(da.coords).keys())
    print(lista)
    for dim in da.dims:
        lista.remove(dim)

    da = da.drop(lista)
    # da = da.drop(['lbl__v_choice', 'lbl__j_choice', 'seq__v_choice', 'seq__j_choice'])
    print(da)
    dimensiones = list(da.dims)
    dimensiones.remove(strEvent)
    dimensiones2 = [strEvent ] +dimensiones
    print( dimensiones2 )

    da = da.transpose(*dimensiones2)
    # print(da)

    #da = da[da.dims]

    df = da.to_dataframe(name='P')
    # print(df) #.head())
    # df = df.set_index('j_choice')
    # df.groupby( by ='j_choice')
    df.reset_index(inplace=True)
    print(df)
    for index, row in df.iterrows():
        print(row.to_dict())
    print(index)
"""


    # db.load_IgorAlignments_FromCSV(strGene,flnAlignments)
    # db.insert_IgorEvent_fromDict(event_dict)
    # for each event in Event_list
    # Add event in table 'MP_Event_list'
    # Add Event realizations in ER_XXXXX

    # change key name
    # dictionary[new_key] = dictionary.pop(old_key)
    # IgorRec_Event
    print( db.flnIgorDB )

    # TODO: Scenarios realizations.
    # mdl.parms.Event_list
    # scenario = p3.IgorScenario()
    # scenario.set_model(mdl)
    #
    # print(scenario.realizations_ids_dict)
    #

    # igor_fln_output_scenarios
    # with open(igor_fln_output_scenarios, "r") as ifile:
    #     inputHeader = ifile.readline()
    #     for line in ifile.readlines():
    #         # print(line.split(strSepChar), len(line.split(strSepChar)) )
    #         line = line.replace("\n", "")
    #         bs = p3.IgorBestScenariosVDJ.load_FromLineBestScenario(line)

    # load Anchors with anchors extract CDR3 and naive alignment

if __name__ == "__main__":
    main()
