import unittest
import pandas as pd
from pygor3 import IgorModel
from pygor3 import IgorTask
from pygor3 import generate
from pygor3 import get_default_IgorModel
from pygor3.utils import get_dataframe_from_fln_generated_seqs_werr
from pygor3 import IgorRec_Event
from pygor3 import IgorEvent_realization

import os
import subprocess

class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.N_seqs = 10
        self.task = IgorTask.default_model("human", "beta")

    def test__something(self):
        nickname = 'v_3_del'

        real = self.task.mdl.parms.get_Event_realization('v_choice', 0)
        print(self.task.mdl.parms['v_choice'].name.loc[0])
        aaa = self.task.mdl.get_Event_value(nickname, 2)
        print(aaa)



    def test__run_generate(self):
        task = IgorTask.default_model("human", "beta")

        self.assertIsInstance(task, IgorTask)
        fln_generated = task._run_generate(10)
        print(task.igor_fln_generated_seqs_werr)
        print(task.igor_fln_generated_realizations_werr)
        print(task.igor_fln_generation_info)

        df = get_dataframe_from_fln_generated_seqs_werr(task.igor_fln_generated_seqs_werr)
        print("df.columns: ", df.columns)
        # TODO: I WANT TO SHOW THIS AS AN AIRR FORMAT
        df2 = pd.read_csv(task.igor_fln_generated_realizations_werr, sep=';').set_index('seq_index')
        print("df2.columns: ", df2.columns)

        events_name__nickname_dict = task.mdl.parms.get_event_dict('name', 'nickname')
        events_nickname__event_type_dict = task.mdl.parms.get_event_dict('nickname', 'event_type')
        events_nickname__seq_type_dict = task.mdl.parms.get_event_dict('nickname', 'seq_type')

        # 1. Replace the name of the columns to the nickname of the columns
        df2.rename(columns=events_name__nickname_dict, inplace=True)

        # task.mdl.genomic_dataframe_dict[str_gene_type.upper()]['name'].loc[x]

        for column_name in df2.columns:
            try:
                if column_name in events_nickname__event_type_dict.keys():
                    if events_nickname__event_type_dict[column_name] == 'GeneChoice':
                        # remove parenthesis and make it an int column
                        df2[column_name] = df2[column_name].apply(lambda x: int(eval(x)))
                        seq_type = events_nickname__seq_type_dict[column_name]
                        str_gene_type = seq_type[0].lower()
                        gene_call_column_name = (str_gene_type+"_call")
                        # df2[gene_call_column_name] = df2[column_name].apply(lambda x: task.mdl.parms.Event_dict[column_name]['name'].loc[x] )
                        df2[gene_call_column_name] = df2[column_name].apply(lambda x: task.mdl.parms[column_name].name.loc[x])
                        df2[column_name].apply(lambda x: task.mdl.parms[column_name].name.loc[x])


                    elif events_nickname__event_type_dict[column_name] == 'Insertion':
                        # Change to insertions values
                        df2[column_name] = df2[column_name].apply(lambda x: int(eval(x)))
                        # df2["value_"+column_name] = df2[column_name].apply(lambda x: task.mdl[column_name]['name'].loc[x])
                        pass
                    elif events_nickname__event_type_dict[column_name] == 'Deletion':
                        # Change to deletions values
                        df2[column_name] = df2[column_name].apply(lambda x: int(eval(x)))
                        pass
                    elif events_nickname__event_type_dict[column_name] == 'DinucMarkov':
                        # Change to deletions values
                        df2[column_name] = df2[column_name].apply(lambda x: list(eval(x)))
                    else:
                        print("Problem column not found!")
                        pass

                event = task.mdl.parms.get_Event(column_name)
                if event.event_type == 'GeneChoice' :
                    print("event.nickname: ", event.nickname)
            except Exception as e:
                pass
        eee = task.mdl.parms.Event_list[0]
        print("eee.seq_type:", eee.seq_type)

        for event in task.mdl.parms.Event_list:
            if event.seq_type == 'GeneChoice':
                print(event.nickname)
        print(df2.head())

        # print("events_name_nickname_dict", events_name_nickname_dict)
        # bs_events_nickname_list = [events_name_nickname_dict[event_name] for event_name in bs_events_name_list]

        task._run_clean_batch_generate()
        task._run_clean_batch_mdldata()

        # task.run_clean_batch()
        # if return_df:
        #     df = pd.read_csv(self.igor_fln_generated_seqs_werr, delimiter=';').set_index('seq_index')
        #     return df

        # pd_sequences = task.generate(self.N_seqs)
        # self.assertIsInstance(pd_sequences, pd.DataFrame)
        # self.assertTrue(len(pd_sequences) == self.N_seqs)

    def test_pygor_generate(self):
        mdl = get_default_IgorModel("human", "tcr_beta")
        df_seqs = generate(10, mdl)
        print(df_seqs)


    def test_something_human_alpha(self):
        task = IgorTask.default_model("human", "alpha")
        self.assertIsInstance(task, IgorTask)
        pd_sequences = task.generate(self.N_seqs)
        self.assertIsInstance(pd_sequences, pd.DataFrame)
        self.assertTrue(len(pd_sequences) == self.N_seqs)

    def test_generate_keep_directories(self):
        task = IgorTask.default_model("mouse", "beta")
        self.assertIsInstance(task, IgorTask)
        pd_sequences = task.generate(self.N_seqs, igor_wd='here', clean_batch=False)
        self.assertIsInstance(pd_sequences, pd.DataFrame)
        self.assertTrue(len(pd_sequences) == self.N_seqs)

        # Check if files exist
        files_to_check_list = [task.igor_fln_generated_realizations_werr,
        task.igor_fln_generated_realizations_werr,
        task.igor_model_parms_file,
        task.igor_model_marginals_file]

        # FIXME: mld_data directory is created but apparently is not using it
        for fln in files_to_check_list:
            print("fln : ", fln)
            self.assertTrue(os.path.isfile(fln))

        subprocess.call("rm -r here", shell=True)


if __name__ == '__main__':
    unittest.main()
