import unittest
from pygor3 import *
import subprocess
import numpy as np
import xarray
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
from pygor3.utils import *

class MyTestCase(unittest.TestCase):
    # def setUp(self) -> None:
    #     self.mdl_hb = IgorModel.load_default("human", "tcr_beta")
    #     N_seqs = 10
    #     self.pd_input_sequences = generate(N_seqs, self.mdl_hb, seed=10)

    def test_realiz(self):
        mdl = get_default_IgorModel("human", "tcr_beta")
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_scenarios(fln_scenarios)

        ps_scenario = df_scenarios.iloc[0]

        event_nickname = 'vd_dinucl'
        id_vd_dinucl = ps_scenario[event_nickname]
        realiz_vd_dinucl = mdl.parms.Event_dict[event_nickname].loc[id_vd_dinucl]
        print("realiz_vd_dinucl: ", type(realiz_vd_dinucl))
        print(realiz_vd_dinucl)
        if isinstance(id_vd_dinucl, pd.Series):
            rrr = IgorEvent_realization.from_tuple(id_vd_dinucl.values, realiz_vd_dinucl.value.values, realiz_vd_dinucl.name.values)
        else:
            rrr = IgorEvent_realization.from_tuple(id_vd_dinucl, realiz_vd_dinucl.value.values,
                                                   realiz_vd_dinucl.name.values)
        print("rrr: ", type(rrr))
        ppp = mdl.realization(ps_scenario, event_nickname)


        event_nickname = 'v_choice'
        id_v_choice = ps_scenario[event_nickname]
        realiz_v_choice = mdl.parms.Event_dict[event_nickname].loc[id_v_choice]
        print("realiz_v_choice: ", type(realiz_v_choice))
        print(realiz_v_choice)



    def test_scenario_fasta(self):
        mdl = get_default_IgorModel("human", "tcr_beta")
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_scenarios(fln_scenarios)

        ps_scenario = df_scenarios.iloc[0]
        get_gene_segment


    def test_gene_segments(self):
        # FIXME: IN DEV
        str_seq = "GGTGCTGTCGTCTCTCAACATCCGAGCTGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTGGGGAAGGGCAGCCTGGAAACACCATATATTTTGGAGAGGGAAGTTGGCTCACTGTTGTAG"
        mdl = get_default_IgorModel("human", "tcr_beta")
        df_scenarios = evaluate(str_seq, mdl, N_scenarios=10, airr_format=False)
        ps_scenario = df_scenarios.iloc[0]

        offset = 0
        list_cols_4_alignment = ["segment_description", "gene_description", "offset", "palindrome_5_end", "gene_ini",
                                 "gene_end", "gene_cut", "palindrome_3_end", "gene_segment"]
        df_scenario_aln = pd.DataFrame(columns=list_cols_4_alignment)
        for ii, strGene in enumerate(['V', 'VD', 'D', 'DJ', 'J']):
            ordered_dicto = mdl.get_gene_segment_dict(strGene, ps_scenario)
            dicto = dict(ordered_dicto)
            dicto['segment_description'] = strGene
            dicto['offset'] = offset
            print(dicto)
            df_scenario_aln.loc[ii] = dicto #, ignore_index=True)
            offset = offset + len(ordered_dicto['gene_segment'])

        df_scenario_aln.aln_scenario_len = offset + len(dicto['gene_segment']) # maximun of (offset + len(gene_segment)) - minimun offset

        V_offset = df_scenario_aln.loc[df_scenario_aln['segment_description'] == 'V'].offset.values[0]
        J_offset = df_scenario_aln.loc[df_scenario_aln['segment_description'] == 'J'].offset.values[0]

        df_scenario_aln.aln_pos_V_anchor = V_offset + mdl.V_anchor(ps_scenario[mdl.parms.event_GeneChoice_V.nickname])
        df_scenario_aln.aln_pos_J_anchor = J_offset + mdl.J_anchor(ps_scenario[mdl.parms.event_GeneChoice_J.nickname])

        print(df_scenario_aln)
        print(df_scenario_aln.aln_pos_V_anchor)
        print(df_scenario_aln.aln_pos_J_anchor)
        print(df_scenario_aln.aln_scenario_len)

    def test_get_df_scenario_aln_from_scenario(self):
        str_seq = "GGTGCTGTCGTCTCTCAACATCCGAGCTGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTGGGGAAGGGCAGCCTGGAAACACCATATATTTTGGAGAGGGAAGTTGGCTCACTGTTGTAG"
        mdl = get_default_IgorModel("human", "tcr_beta")
        df_scenarios = evaluate(str_seq, mdl, N_scenarios=10, airr_format=False)
        fln_scenario_fasta = "aver.fasta"
        ps_scenario = df_scenarios.iloc[0]
        mdl.write_df_scenario_aln_FASTA(fln_scenario_fasta, ps_scenario)
        with open(fln_scenario_fasta, 'a') as ofile:
            ofile.write("> input_sequence \n")
            ofile.write(str_seq+"\n")

    def test_numpy_aln_matrix(self):
        # str_seq = "GGTGCTGTCGTCTCTCAACATCCGAGCTGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTGGGGAAGGGCAGCCTGGAAACACCATATATTTTGGAGAGGGAAGTTGGCTCACTGTTGTAG"
        # mdl = get_default_IgorModel("human", "tcr_beta")
        mdl = get_default_IgorModel("human", "tcr_alpha")
        seqs = generate(1, mdl, seed=0)
        str_seq = seqs.loc[0]['nt_sequence']

        df_scenarios = evaluate(str_seq, mdl, N_scenarios=10, airr_format=False)
        ps_scenario = df_scenarios.iloc[0]
        df_scenario_aln = mdl.get_df_scenario_aln_from_scenario(ps_scenario)
        da = from_df_scenario_aln_to_da_scenario_aln(df_scenario_aln)
        da.plot()
        plt.show()



        # nrows = len(df_scenario_aln.index)
        # ncols = df_scenario_aln.aln_scenario_len
        # aln_scenario_np = -np.ones((nrows, ncols))
        # dict_id_2_nt = {-1: '-', 0: 'A', 1: 'C', 2: 'G', 3: 'T'} #mdl.parms['vd_dinucl']['value'].to_dict()
        # dict_nt_2_id = {v: k for k, v in dict_id_2_nt.items()}
        # for ii, row in df_scenario_aln.iterrows():
        #     # add a map with a dictionary defined
        #     np_gene_segment = np.array(list( map(lambda nt: dict_nt_2_id[nt], list(row['gene_segment']) ) ) )
        #     aln_ini_gene_segment = row['offset']
        #     aln_end_gene_segment = row['offset'] + len(row['gene_segment'])
        #     aln_scenario_np[ii, aln_ini_gene_segment:aln_end_gene_segment] = np_gene_segment
        #
        # print(aln_scenario_np)

    def test_plot_scenario_from_da_scenario_aln(self):
        # str_seq = "GGTGCTGTCGTCTCTCAACATCCGAGCTGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTGGGGAAGGGCAGCCTGGAAACACCATATATTTTGGAGAGGGAAGTTGGCTCACTGTTGTAG"
        mdl = get_default_IgorModel("human", "tcr_alpha")
        seqs = generate(1, mdl, seed=0)
        str_seq = seqs.loc[0]['nt_sequence']

        df_scenarios = evaluate(str_seq, mdl, N_scenarios=10, airr_format=False)
        ps_scenario = df_scenarios.iloc[0]
        print("ps_scenario: ", ps_scenario)
        df_scenario_aln = mdl.get_df_scenario_aln_from_scenario(ps_scenario)
        da_scenario_aln = from_df_scenario_aln_to_da_scenario_aln(df_scenario_aln)
        print(da_scenario_aln.sizes)
        print(da_scenario_aln['nucleotide'].size)
        print(da_scenario_aln['nucleotide'].shape)
        fig, ax = plot_scenario_from_da_scenario_aln(da_scenario_aln, nt_lim=(240, 320))
        plt.show()

    def test_IgorModel_plot_scenario(self):
        mdl = get_default_IgorModel("human", "tcr_alpha")
        seqs = generate(1, mdl, seed=0)
        str_seq = seqs.loc[0]['nt_sequence']

        df_scenarios = evaluate(str_seq, mdl, N_scenarios=10, airr_format=False)
        ps_scenario = df_scenarios.iloc[0]
        print("ps_scenario: ", ps_scenario)
        mdl.plot_scenario(ps_scenario)
        df_scenario_aln = mdl.get_df_scenario_aln_from_scenario(ps_scenario)
        da_scenario_aln = from_df_scenario_aln_to_da_scenario_aln(df_scenario_aln)
        print(da_scenario_aln.sizes)
        print(da_scenario_aln['nucleotide'].size)
        print(da_scenario_aln['nucleotide'].shape)
        fig, ax = plot_scenario_from_da_scenario_aln(da_scenario_aln, nt_lim=(240, 320))
        plt.show()



    def test_bahal(self):
        colors = ['white', 'red', 'blue', 'orange', 'green']
        dna_values = ['-'] + list(mdl.parms['vd_dinucl']["value"].values)
        cmap_dna = mpl.colors.ListedColormap(['white', '#fcff92', '#70f970', '#ff99b1', '#4eade1'])
        # cmap_dna

        # fig, ax = plt.subplots(figsize=(40, 40))
        # ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        # ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # ax.set_aspect('equal')

        fig, ax = plt.subplots(figsize=(200, 10))
        ax.imshow(aln_scenario_np, cmap=cmap_dna, vmin=-1.5, vmax=3.5)

        da_scenario_aln.attrs["anchors_CDR3"][1]

        anchors_CDR3_ini = da_scenario_aln.attrs["anchors_CDR3"][0]
        anchors_CDR3_end = da_scenario_aln.attrs["anchors_CDR3"][1]

        ax.axvline(anchors_CDR3_ini-0.5)
        ax.axvline(anchors_CDR3_end-0.5)

        xlim_ini = anchors_CDR3_ini - 3.5
        xlim_end = anchors_CDR3_end + 2.5
        ax.set_xlim(xlim_ini, xlim_end)
        # ax.set_xlim(250, 320)

        # TODO: MODIFY THIS TO GET
        dict_id_2_nt = Igor_dict_id_2_nt
        dict_4_lbls = dict_id_2_nt.copy()
        dict_4_lbls[-1] = ''
        for i in range(aln_scenario_np.shape[0]):
            for j in range(int(xlim_ini+0.5), int(xlim_end+0.5)): #aln_scenario_np.shape[1]):
                text = ax.text(j, i, dict_4_lbls[aln_scenario_np[i, j]], ha="center", va="center", color="gray")

        da_scenario_aln['gene_description'].values
        # ax.set_xticks(np.arange(len(farmers)))
        ax.set_yticks(np.arange(len(da_scenario_aln['gene_description'].values)))
        # ... and label them with the respective list entries
        # ax.set_xticklabels(farmers)
        ax.set_yticklabels(da_scenario_aln['gene_description'].values)



        # da_aln = xr.DataArray(aln_np)
        # da_aln.plot()
        plt.show()


        # mdl.get_gene_segment_dict('D', ps_scenario)['gene_segment']
        # mdl.get_gene_segment_dict('J', ps_scenario)['gene_segment']
        #
        # mdl.get_gene_segment_dict('VD', ps_scenario)['gene_segment']
        # mdl.get_gene_segment_dict('DJ', ps_scenario)['gene_segment']
        # so now I need to use this segments to make a fasta

    def test_mutual_information(self):
        mdl = get_default_IgorModel("human", "tcr_beta")
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_scenarios(fln_scenarios)

        event_nickname_x = 'v_choice'
        event_nickname_y = 'j_choice'
        P_x_y = mdl.get_P_from_scenarios_cols(df_scenarios, [event_nickname_x, event_nickname_y])
        P_x = mdl.get_P_from_scenarios_cols(df_scenarios, [event_nickname_x])
        P_y = mdl.get_P_from_scenarios_cols(df_scenarios, [event_nickname_y])

        I_X_Y = get_D_KL_from_xarray(P_x_y, P_x, P_y)
        print("I_X_Y: ", I_X_Y)

        # print(P_x * P_y.values)
        # print( np.matmul(P_x, P_y) )
        func = lambda xy, x, y: np.log(xy/(x*y))
        log_Pxy_over_Px_Py = xr.apply_ufunc(func, P_x_y, P_x, P_y)
        log_Pxy_over_Px_Py
        print(log_Pxy_over_Px_Py)
        print(np.nan_to_num(log_Pxy_over_Px_Py.values))
        log_Pxy_over_Px_Py.values = np.nan_to_num(log_Pxy_over_Px_Py.values, neginf=0)
        print(log_Pxy_over_Px_Py)
        mutual_information = xr.dot(P_x_y, log_Pxy_over_Px_Py)
        print("mutual_information: ", mutual_information)

    def test_mutual_information_matrix(self):
        # FIXME: DEV
        mdl = get_default_IgorModel("human", "tcr_beta")
        event_nickname1 = 'v_choice'
        event_nickname2 = 'vd_ins'
        da_P_x_y = mdl.get_P_joint([event_nickname1, event_nickname2])
        da_P_x = mdl.Pmarginal[event_nickname1]
        da_P_y = mdl.Pmarginal[event_nickname2]
        # mi = get_D_KL_from_xarray(da_P_x_y, da_P_x, da_P_y)

        # np_P_x_joint_y = np.matmul(da_P_x.values[np.newaxis].T, da_P_y.values[np.newaxis])  # , da_P_x_y.values
        # da_P_x.values, da_P_y.values.T
        # rl = np.outer(np.ones((5,)), np.linspace(-2, 2, 5))
        # rl
        print('-' * 50)
        mi = mdl.get_mutual_information_events(event_nickname1, event_nickname2)
        print(mi)

        event_nickname1 = 'v_choice'
        event_nickname2 = 'j_choice'


        da_P_x_y = mdl.get_P_joint([event_nickname1, event_nickname2])
        da_P_x = mdl.Pmarginal[event_nickname1]
        da_P_y = mdl.Pmarginal[event_nickname2]

        da_P_x_times_P_y = (da_P_x * da_P_y)
        da_P_x_times_P_y

        da_log_P_ratio = xr.zeros_like(da_P_x_y)

        da_log_P_ratio.values = np.nan_to_num(
            np.log2(da_P_x_y / da_P_x_times_P_y), nan=0.0, neginf=0.0
        )

        # da_log_Value.values = np_log_Value
        xr.dot(da_P_x_y, da_log_P_ratio)

        # entropy = mdl.get_entropy_event('v_choice')
        # print("entropy: ", entropy)



        # fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        # df_scenarios = mdl.get_dataframe_scenarios(fln_scenarios)
        #
        #
        # da_mi_scenarios = mdl.get_mutual_information_from_df_scenarios(df_scenarios)
        # mdl.plot_mutual_information(da_mi_scenarios, ax=ax[1][0])
        # # oops =
        #
        # print(da_mi_scenarios)
        # da_oops = xr.zeros_like(da_mi_mdl)
        # print(da_oops)
        # da_oops.values = np.tril(da_mi_mdl) + np.triu(da_mi_scenarios)
        # mdl.plot_mutual_information(da_oops, ax=ax[0][1])
        #
        # plt.show()


    def test_stuffs(self):
        mdl = get_default_IgorModel("human", "tcr_beta")
        mdl.get_events_nicknames_list()
        da_P_joint = mdl.get_P_joint(['v_choice'])
        print(da_P_joint)
        da_P_joint.plot()
        dict_nickname_event_type = self.parms.get_event_dict('nickname', 'event_type')
        dict_events = {key: val for key, val in dict_nickname_event_type.items() if val != 'DinucMarkov'}
        event_lista_nicknames = list(dict_events.keys())
        data_0 = np.zeros((len(event_lista_nicknames), len(event_lista_nicknames)))
        da_mi = xr.DataArray(data_0, dims=('x', 'y'), coords={'x': event_lista_nicknames, 'y': event_lista_nicknames})
        da_mi.name = 'mutual_information'

        mdl.get_mutual_information()
        plt.show()

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
