#!/usr/bin/env python3
"""
Export Igor Model to csv files.
"""

import pygor3 as p3
import pandas as pd


from optparse import OptionParser

def main():
    task = p3.IgorTask()
    parser = OptionParser()
    parser.add_option("-s", "--species", dest="species", help='Igor species')
    parser.add_option("-c", "--chain", dest="chain", help='Igor chain')
    parser.add_option("-o", "--output", dest="output", help='filename of csv file to export data')

    parser.add_option("-p", "--model_params", dest="model_params", help='IGoR model_params.txt')
    parser.add_option("-m", "--model_marginals", dest="model_marginals", help='IGoR model_marginals.txt')


    (options, args) = parser.parse_args()

    if options.output == None:
        print("output option is mandatory")
        exit()
    else:
        flnPrefix = options.output

    Q_species_chain = ( not (options.species == None) and not(options.chain == None))
    Q_model_files = ( not (options.model_params == None) and not( options.model_marginals == None))
    if ( Q_species_chain or Q_model_files):

        import numpy as np
        import matplotlib.pyplot as plt

        if Q_species_chain :
            task.igor_species = str(options.species)
            task.igor_chain = str(options.chain)
            task.load_IgorModel()

            for event_nickname in task.mdl.Pmarginal.keys():
                fig, ax = plt.subplots()
                # df = task.mdl.Pmarginal[event_nickname].to_dataframe(name=event_nickname)
                task.mdl.plot_Event_Marginal(event_nickname, ax=ax)
                fig.tight_layout()
                flnOutput = flnPrefix+"_"+event_nickname+".pdf"
                fig.savefig(flnOutput)
                print("***** Marginal plot of ", event_nickname, " in ", flnOutput)
                # if task.mdl.parms.get_Event(event_nickname).event_type == 'Dinucl':
                #     df.plot.bar(y=event_nickname, x='lbl__'+event_nickname, ax=ax)
                #     fig.savefig("PM_"+event_nickname+".pdf")


        if Q_model_files:
            task.igor_model_parms_file = options.model_params
            task.igor_model_marginals_file = options.model_marginals
            task.load_IgorModel()

            for event_nickname in task.mdl.Pmarginal.keys():
                fig, ax = plt.subplots()
                # df = task.mdl.Pmarginal[event_nickname].to_dataframe(name=event_nickname)
                task.mdl.plot_Event_Marginal(event_nickname, ax=ax)
                fig.tight_layout()
                flnOutput = flnPrefix+"_"+event_nickname+".pdf"
                fig.savefig(flnOutput)
                print("***** Marginal plot of ", event_nickname, " in ", flnOutput)

            task.mdl.export_csv(options.output)

    else:
        if not Q_species_chain :
            print("species and chain is mandatory")
            exit()

        if not Q_model_files :
            print("models files are mandatory")
            exit()

    # if options.chain == None:
    #     print("species is mandatory")
    #     exit()
    # else:
    #     task.igor_chain = str(options.chain)
    #
    # if options.output == None:
    #     print("output is mandatory")
    #     exit()
    # else:
    #     output = options.output
    #
    # if options.batch == None:
    #     print("Batchname not set. Random batchname will be assing.")

    # if len(args) > 0 and len(args) <2:
    #   filename = str(args[0])
    #   task.igor_read_seqs = filename
    # else:
    #   print("Supply a correct filename")
    #   exit()

    # task.run_evaluate()
    # task.update_batch_filenames()

    #
    # for strEvent in task.mdl.xdata.keys():
    #     da = task.mdl.xdata[strEvent]
    #     dependencias = list(task.mdl.xdata[strEvent].dims)
    #     print("********", dependencias, strEvent)
    #     dependencias.remove(strEvent)
    #     dependencias_dim = [task.mdl.xdata[strEvent][dep].shape[0] for dep in dependencias]
    #     print(strEvent, da.dims, dependencias_dim, dependencias)
    #     if len(dependencias) == 0 :
    #         lbl_file = "P__"+strEvent
    #         print(lbl_file)
    #     elif len(dependencias) == 1:
    #         lbl_file = "P__" + strEvent +"__g__"+dependencias[0]
    #         print(lbl_file)
    #         df = pd.DataFrame(data=da.values, index=da['lbl__' + dependencias[0]].values,
    #                           columns=da['lbl__' + strEvent].values)
    #         #print(df)
    #     elif len(dependencias) == 2 :
    #         lbl_file = "P__" + strEvent + "__g__" + dependencias[0]+"__"+ dependencias[1]
    #         print(lbl_file)
    #     elif len(dependencias) == 3:
    #         lbl_file = "P__" + strEvent + "__g__" + dependencias[0] + "__" + dependencias[1]+ "__" + dependencias[2]
    #         print(lbl_file)
    #     else:
    #         print("OHHHHHHHH")
    #     #df = pd.DataFrame(data=da.values, index=da['lbl__' + 'v_choice'].values, columns=da['lbl__' + 'j_choice'].values)
    #
    # task.run_clean_batch()

if __name__ == "__main__":
    main()
