#!/usr/bin/env python3
import pygor3 as p3
import argparse

def main():
    parser = argparse.ArgumentParser()
    igor_models = parser.add_argument_group('IGoR default models')

    parser.add_argument("-D", "--database", dest="database", help="Igor database created with database script.",
                        default="Ajam.db")

    args = parser.parse_args()


    db = p3.IgorSqliteDB.create_db(args.database)
    mdl = db.get_IgorModel()
    Prob = mdl.xdata['v_choice']*mdl.xdata['j_choice']
    print(Prob)
    # mdl0 = p3.IgorModel.load_default("human", "tcr_beta")

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
