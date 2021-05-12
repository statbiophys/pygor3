import unittest
from pygor3 import *
import subprocess
import numpy as np
import xarray
import itertools
import matplotlib.pyplot as plt

class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.mdl_hb = IgorModel.load_default("human", "tcr_beta")
        N_seqs = 10
        self.pd_input_sequences = generate(N_seqs, self.mdl_hb, seed=10)

    def test_sequences_generation(self):
        self.assertIsInstance(self.pd_input_sequences, pd.DataFrame)
        task = IgorTask.default_model("human", "beta")
        print("task.b_align: ", task.b_align)
        with tempfile.NamedTemporaryFile(dir='.', prefix=task.igor_batchname, suffix='.csv') as tmp_fln:
            self.pd_input_sequences.to_csv(tmp_fln.name, sep=';')
            # task._run_evaluate(tmp_fln.name, N_scenarios=3)
            # TODO: FOR SOME REASON THE DATABASE IS NOT UPDATING THE FILENAME AND I NEED THAT
            #  TO CREATE A DATABASE AND CALCULATE THE OBSERVABLES WITH THE DATABASE.
            #  ALSO I NEED TO MAKE A NEW MATRIX TO CALCULATE THE OBSERVABLES.
            # df = task.get_dataframe_from_fln_generated_seqs_werr()
            # df = task.evaluate(self.pd_input_sequences, N_scenarios=3

            N_scenarios = 5
            # 3. Write Sequences in file if file not exist
            fln_input_sequences = task.igor_wd + "/" + task.igor_batchname + "input_sequences.csv"
            write_sequences_to_file(self.pd_input_sequences, fln_input_sequences)

            # 4. Export model and ref_genome to model_dir
            path_mdl_data = task.igor_wd + "/" + task.igor_batchname + "_mdldata"
            task.update_model_filenames(igor_model_dir_path=path_mdl_data)
            task.update_ref_genome(igor_model_dir_path=path_mdl_data)
            task.update_batch_filenames()
            task.mdl.write_mdldata_dir(path_mdl_data)

            # 5. Run evaluate model
            task._run_evaluate(igor_read_seqs=fln_input_sequences, N_scenarios=N_scenarios)
            print(task.igor_fln_output_scenarios)


            df_scenarios = task.mdl.get_dataframe_from_fln_generated_realizations_werr(task.igor_fln_output_scenarios)
            print(df_scenarios)


            # 1. Define your observable variable list
            event_tupla_list = [('vd_ins', 'value'), ('v_3_del', 'value')]
            event_tupla_list = [('vd_ins', 'value'), ('v_3_del', 'value')]
            event_lista = [var[0] for var in event_tupla_list]

            # 2. Get the probabilities from scenarios
            da_events_prob = task.mdl.get_probability_matrix_from_event_list_and_scenarios_dataframe(event_lista, df_scenarios)

            # 3. define your observable
            def observable_f0(x, y):
                return x ** 2 + y ** 4

            da_observable = self.mdl_hb.get_observable_xarray_from_function(observable_f0, event_tupla_list)

            # 4. calculate the mean value
            tmp_da = np.multiply(da_observable, da_events_prob)
            observable_value = tmp_da.sum()
            print(observable_value)




            fig, ax = plt.subplots(2,1)
            da_observable.plot(ax=ax[0])
            da_events_prob.plot(ax=ax[1])
            plt.show()

    def test_get_realizations(self):
        df_scenarios = evaluate(self.pd_input_sequences, self.mdl_hb, N_scenarios=3, airr_format=False)
        ps_scenario = df_scenarios.iloc[0]
        print(ps_scenario)
        v_choice = self.mdl_hb.realization(ps_scenario, 'v_choice')
        print(v_choice)
        self.assertIsInstance(v_choice, IgorEvent_realization)
        vd_dinucl = self.mdl_hb.realization(ps_scenario, 'vd_dinucl')
        print(vd_dinucl)
        # self.assertIsInstance(v_choice, IgorEvent_realization)

    def test_observables_realizations(self):
        def observable_n_aa_in_vd(ps_scenario):
            """
            Return the numbers of amino acids in vd insertions
            """
            vd_dinucl = mdl.pd_realization(ps_scenario, 'vd_dinucl')




    def test_observable_02(self):
        # 0. get your sequences
        mdl = IgorModel.load_default("human", "tcr_beta")
        sequences = generate(100, mdl=mdl)
        df_scenarios = evaluate(sequences, mdl, N_scenarios=3, airr_format=False)
        print(df_scenarios)

        # 1. Define your observable variable list
        #event_tupla_list = [('vd_ins', 'value'), ('v_3_del', 'value')]
        # event_tupla_list = [('v_choice', 'id'), ('v_3_del', 'value')]
        event_tupla_list = [('vd_ins', 'value'), ('dj_ins', 'value'), ('v_3_del', 'value')]

        event_lista = [var[0] for var in event_tupla_list]

        # 2. Get the probabilities from scenarios
        da_events_prob = mdl.get_probability_matrix_from_event_list_and_scenarios_dataframe(event_lista,
                                                                                                 df_scenarios)

        # 3. define your observable
        def observable_f0(x, y):
            return x ** 2 + y ** 4

        def observable_f1(x, y):
            if x == 2:
                return y ** 2
            else:
                return 0

        def observable_f2(vd_ins, dj_ins, v_3_del):
            return vd_ins + dj_ins + v_3_del


        da_observable = self.mdl_hb.get_observable_xarray_from_function(observable_f2, event_tupla_list)

        # 4. calculate the mean value
        da_observable_weighted = np.multiply(da_observable, da_events_prob)
        observable_values = da_observable_weighted.sum()
        print(da_events_prob)
        print(da_observable)

        fig, ax = plt.subplots(2, 2)

        da_observable.plot(ax=ax[0][0])
        da_events_prob.plot(ax=ax[1][0])
        da_observable_weighted.plot(ax=ax[0][1])
        # observable_values.plot(ax=ax[1][1])

        ax[0][0].set_title("da_observable")
        ax[1][0].set_title("da_events_prob")
        ax[0][1].set_title("da_observable_weighted")

        plt.show()

    def test_observables_pandas(self):
        # TODO: IN DEV

        # 0. get your sequences
        mdl = IgorModel.load_default("human", "tcr_beta")
        sequences = generate(10, mdl=mdl, seed=10)
        df_scenarios = evaluate(sequences, mdl, N_scenarios=3, airr_format=False)
        print(mdl.V_anchors)
        print("-"*50)
        print("mdl.V_anchor(2) : ", mdl.V_anchor(2) )
        sequences = generate(1000, mdl=mdl, seed=10)

        task = IgorTask(mdl=mdl, igor_wd="delete_me", igor_batchname='ttmmpp')
        task.evaluate(sequences, clean_batch=False)
        task.get_dataframe_scenarios()
        task.get_dataframe_from_fln_generated_realizations_werr()
        print("*"*50)
        print(task.to_dict())


        df_scenarios = evaluate(sequences, mdl, N_scenarios=3, airr_format=False)
        print(df_scenarios.columns)

        df_scenarios['norm_scenario_proba_cond_seq'] = get_df_normalize_prob(df_scenarios)
        # print(df_scenarios)

        ps_scenario = df_scenarios.iloc[2]
        mdl.realization(ps_scenario, 'vd_dinucl')

        # aaa = df_scenarios.apply(lambda x: mdl.parms['v_choice'].loc[x['v_choice']].value )
        def observable_001(id_v_choice):
            v_choice = mdl.parms.get_Event_realization('v_choice', id_v_choice)
            v_anchor = mdl.V_anchors.loc[id_v_choice]
            return len(v_choice.value) - v_anchor

        def observable_V_segment(ps_scenario):
            v_choice = mdl.realization(ps_scenario, 'v_choice')
            v_3_del = mdl.realization(ps_scenario, 'v_3_del')
            try:
                # v_anchor = mdl.V_anchors.loc[v_choice.id]
                v_anchor = mdl.V_anchor(v_choice.id)
            except Exception as e:
                return None # FIXME: what to do in these cases, this will change the normalization.
            return len(v_choice.value) - v_anchor - v_3_del.value


        def observable_CDR3_len(ps_scenario):
            v_choice = mdl.realization(ps_scenario, 'v_choice')
            v_3_del = mdl.realization(ps_scenario, 'v_3_del')
            d_gene = mdl.realization(ps_scenario, 'd_gene')
            vd_ins = mdl.realization(ps_scenario, 'vd_ins')
            d_5_del = mdl.realization(ps_scenario, 'd_5_del')
            d_3_del = mdl.realization(ps_scenario, 'd_3_del')
            dj_ins = mdl.realization(ps_scenario, 'dj_ins')
            j_5_del = mdl.realization(ps_scenario, 'j_5_del')
            j_choice = mdl.realization(ps_scenario, 'j_choice')

            try:
                # v_anchor = mdl.V_anchors.loc[v_choice.id]
                v_anchor = mdl.V_anchor(v_choice.id)
                j_anchor = mdl.V_anchor(j_choice.id)
                cdr3_V_len = len(v_choice.value) - v_anchor - v_3_del.value
                cdr3_VD_len = vd_ins.value
                cdr3_D_len = len(d_gene.value) - d_5_del.value - d_3_del.value
                cdr3_DJ_len = dj_ins.value
                cdr3_J_len = j_anchor - j_5_del.value

                return cdr3_V_len + cdr3_VD_len + cdr3_D_len + cdr3_DJ_len + cdr3_J_len

            except Exception as e:
                return None # FIXME: what to do in these cases, this will change the normalization.

        # vd_dinucl = mdl.realization(ps_scenario, 'vd_dinucl')




        def observable_n_Cys_in_vd(ps_scenario):
            """
            Return the numbers of amino acids in vd insertions
            """
            v_choice = mdl.realization(ps_scenario, 'v_choice')
            v_3_del = mdl.realization(ps_scenario, 'v_3_del')
            vd_ins = mdl.realization(ps_scenario, 'vd_ins')
            v_segment_length = len(v_choice.value) - v_3_del.value
            n_extra_nt = v_segment_length % 3
            n_nt = (3 - n_extra_nt) % 3
            if (vd_ins < n_nt):
                return 0
            else:
                return (vd_ins - n_nt) // 3


        # TODO: make a method from mdl.get_realization(nickname, pd_row)

        # value = observable_V_segment(df_scenarios.loc[0])
        # print("value: ", value)

        df_scenarios['observable'] = df_scenarios.apply(lambda row: observable_V_segment(row), axis=1)
        mean_value = np.multiply(df_scenarios['observable'], df_scenarios['norm_scenario_proba_cond_seq']).sum()

        hist, bin_edges = np.histogram(df_scenarios['observable'])
        bins = 0.5*(bin_edges[:-1] + bin_edges[1:])
        fig, ax = plt.subplots()
        ax.plot(bins, hist)
        plt.show()

        # print(help(np.histogram))
        print(help(IgorModel.load_default))
        print(mean_value)



    def test_observable_mdl(self):
        # TODO: DFAS
        # 1. Define your observable
        def observable_V_segment(ps_scenario):
            v_choice = mdl.get_IgorEvent_realization_for_nickname(ps_scenario, 'v_choice')
            v_3_del = mdl.get_IgorEvent_realization_for_nickname(ps_scenario, 'v_3_del')
            try:
                v_anchor = mdl.V_anchors.loc[v_choice.id]
            except KeyError as e:
                return None
            return len(v_choice.value) - v_anchor - v_3_del.value





    def test_foo(self):
        da = self.mdl_hb.get_ones_xarray_from_list(['vd_ins', 'v_3_del'])
        # print(da[{'vd_ins':7}].coords['vd_ins'])
        dimensions_names = list(da.dims)
        coordinates = list( map(lambda x: da.coords[x].values, dimensions_names) )
        # observable = lambda x:
        # fff = np.meshgrid(*coordinates)
        # print(*fff)

        # this should be my definition of observable

        # 1. What variables (nicknames) I'm gonna use?
        # 2. Calculate the probability matrix to use over the observable

        # 3. From the scenarios IgorScenario define my observable
        def observable_f0(x, y):
            return x**2 + y**4

        print(self.mdl_hb.parms.Event_dict['v_3_del'])
        da_observable = self.mdl_hb.get_observable_xarray_from_function(observable_f0, [('vd_ins', 'value'), ('v_3_del', 'value')])
        print("================= da_observable =================")
        print(da_observable)

        # 4. Get an empty or zero matrix with values

        # scenario['vd_ins'] should give me the id of event.

        # 1. define your number of elements list of variables





        # for elem in itertools.product(*coordinates):
        #     print(elem)
        #     elem_dict = dict(zip(dimensions_names, elem))
        #     da[elem_dict]  = self.mdl_hb.parms
        # print(da.dims)
        #

        # itertools.product()
        # print(da.coords['vd_ins'])

        # for elem in da:
        #     print("--------elem------")
        #     print(elem)


        # array = xr.DataArray([1, 2, 3], coords=[("x", [0.1, 0.2, 0.3])])

        # def observable(a):
        #     func = lambda x: x.coords['vd_ins']
        #     return xr.apply_ufunc(func, a, vectorize=True)
        #
        # aaa = observable(da)
        # print(aaa)

        # event_lista = ['v_choice', 'j_choice', 'v_3_del']
            # event_lista = ['v_choice', 'v_3_del']
            # da_vj_zero = task.mdl.get_zero_xarray_from_list(event_lista)
            # df['norm_scenario_proba_cond_seq'] = get_df_normalize_prob(df)
            # # Now the mean value of
            # print(df[event_lista + ['scenario_proba_cond_seq', 'norm_scenario_proba_cond_seq']])
            # df['norm_scenario_proba_cond_seq']
            # aaa = df.groupby(event_lista)['norm_scenario_proba_cond_seq'].apply(lambda x: x.sum())
            #
            # print(aaa.to_dict())
            # print("aaa.index:", aaa.index)
            # print("aaa.index.names:", aaa.index.names)
            # for iii, value in aaa.iteritems():
            #     coordenadas = dict(zip(aaa.index.names, iii))
            #     print(coordenadas, value)  # da_vj_zero[coordenadas])
            #     da_vj_zero[coordenadas] = value
            #
            # test_value = {'v_choice': 30, 'j_choice': 0, 'v_3_del': 6}  # 0.047007
            # test_value = {'v_choice': 30, 'v_3_del': 6}
            # print(da_vj_zero[test_value])
            # print(da_vj_zero.sum())
            # import matplotlib.pyplot as plt
            # da_vj_zero.plot()
            # plt.show()

            # da_prob_matrix = aaa.to_xarray().fillna(0)

            # bbb = da_prob_matrix.combine_first(da_vj_zero)
            # print(bbb[test_value])
            # print(bbb.sum())

            # ccc = da_vj_zero.update(da_prob_matrix)
            # print(ccc)
            # print(ccc.sum())

            # da_event_prob_matrix = da_vj_zero + da_prob_matrix
            # print(da_event_prob_matrix)
            # print(da_event_prob_matrix.sum())
            # for indx, row in aaa.iterrows():
            #     print('%'*50)
            #     print(indx, row)

            # # print(df.columns)
            # print(df.head())
            # # TODO: SO IN NEED TO for each scenario prob
            # df_ps = df.loc[0]
            # print("-"*50)
            # for indx, row in df_ps.iterrows():
            #     print("."*20)
            #     scen = IgorScenario.load_from_dict(row.to_dict())
            #     scen.seq_index = indx
            # print("scen: ", scen.to_dict())
            #
            # print(row)
            # scenario_normalization = ps['scenario_proba_cond_seq'].sum()
            # print(ps['scenario_proba_cond_seq']/ scenario_normalization)
            # print(task.mdl['v_3_del'])

            # Save evaluations in database
            # task.create_db()
            # task.load_db_from_indexed_sequences()
            # task.load_db_from_indexed_cdr3()
            # task.load_db_from_genomes()
            # task.load_db_from_alignments()
            # task.load_IgorModel()
            # task.load_db_from_models()
            # task.load_db_from_bestscenarios()
            # task.load_db_from_pgen()

    def test_something(self):
        task = IgorTask.default_model("human", "beta")
        task._run_generate(10, seed=10)

        # task.igor_fln_generated_seqs_werr
        # task.igor_fln_generated_realizations_werr
        # task.igor_fln_generation_info
        # cmd = "cat " + task.igor_fln_generated_seqs_werr
        # print(cmd)
        # subprocess.run(cmd, shell=True)
        # cmd = "cat " + task.igor_fln_generated_realizations_werr
        # print(cmd)
        # subprocess.run(cmd, shell=True)

        pd_realizations = task.get_dataframe_from_fln_generated_realizations_werr()
        # TODO: .replace("(", "[").replace(")", "]")
        print(pd_realizations.columns)
        column_name = 'dj_dinucl'
        # pd_realizations[column_name] = pd_realizations[column_name].apply(lambda x: list(eval(x)))
        print(pd_realizations)
        print(pd_realizations.loc[0])

        sce = IgorScenario()
        print(sce.to_dict())
        self.assertEqual(True, True)


    def test_something02(self):
        # 1. First get a model
        mdl = IgorModel.load_default("human", "tcr_beta")
        db = IgorSqliteDB()
        db.calc_IgorBestScenarios_average_of()

        # Make a database file with the evaluation to get the observable.

        # 2. Define a function to act over an IgorScenario variable, given the model

        # 3.
        da = mdl.get_zero_xarray_from_list(['v_choice', 'v_3_del'])
        print("-" * 10)
        print(da)

    def test_function_to_observable(self):
        pass

        """
        def get_pairwise_prob(mdl, event_nickname1, event_nickname2):
            # create an xarray matrix with event_nickname1 and event_nickname2
            da = mdl.get_zero_xarray_from_list([event_nickname1, event_nickname2])
            # print("da = ", da)
            import xarray as xr

            def tmp_funct(bs:IgorScenario):
                ev_dict_ids = {event_nickname1: bs['id_' + event_nickname1],
                               event_nickname2: bs['id_' + event_nickname2]}
                darr = xr.zeros_like(da)
                darr.loc[ev_dict_ids] = 1
                return darr

            return tmp_funct
        """

    def test_IgorModel(self):
        task = IgorTask.default_model("human", "beta", igor_wd='igor_temporal')
        task._run_generate(5, seed=10)
        pd_sequences = task.get_dataframe_from_fln_generated_seqs_werr()
        pd_realizations = task.get_dataframe_from_fln_generated_realizations_werr()

        # pd_realizations = task.mdl.get_dataframe_from_fln_generated_realizations_werr(
        #     task.igor_fln_generated_realizations_werr)
        print(pd_realizations)

        # task.load_db_from_indexed_sequences()
        # task.load_IgorModel()
        # task.load_db_from_genomes()
        # task.load_db_from_models()
        # task.load_db_from_bestscenarios()

        task.evaluate(pd_sequences, N_scenarios=3, clean_batch=False)
        task.igor_fln_output_scenarios

        # task._run_generate(10, seed=20)
        # subprocess.run("rm -r igor_temporal", shell=True)


# 18110919b366adc767a015deb72ab7f8e409476e

if __name__ == '__main__':
    unittest.main()
