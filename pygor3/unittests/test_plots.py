import unittest

from pygor3 import *
import itertools
import matplotlib.pyplot as plt
class MyTestCase(unittest.TestCase):
    def test_something(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        fig, ax = mdl.plot_Event('v_choice')
        self.assertIsInstance(fig, plt.Figure)
        fig, ax = mdl.plot_Event('d_gene')
        self.assertIsInstance(fig, plt.Figure)
        fig, ax = mdl.plot_Event('j_choice')
        self.assertIsInstance(fig, plt.Figure)
        fig, ax = mdl.plot_Event('vd_ins')
        self.assertIsInstance(fig, plt.Figure)
        fig, ax = mdl.plot_Event('dj_ins')
        self.assertIsInstance(fig, plt.Figure)
        fig, ax = mdl.plot_Event('vd_dinucl')
        self.assertIsInstance(fig, plt.Figure)
        fig, ax = mdl.plot_Event('dj_dinucl')
        self.assertIsInstance(fig, plt.Figure)

    def test_from_df_scenario_aln_to_da_scenario_aln(self):

        # FIXME: DEV A

        mdl = IgorModel.load_default("human", "tcr_beta")
        mdl.plot_scenario()
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_scenarios(fln_scenarios)

        # print(np.unique(df_scenarios[df_scenarios['v_3_del'] < 4].index ))#.groupby('seq_index'))
        # [ 56  69  87 108 119 138 150 151 157 178 194 197 202 204 208 213 247 252
        #  253 272 283 308 323 352 353 363 371 398 442 471 477 481 482 491 550 561
        #  571 588 590 594 634 638 641 645 706 708 726 731 740 749 754 787 802 809
        #  857 862 864 869 881 883 899 928 958 982 983]

        df_scenario = df_scenarios.loc[87] # 982, 706
        ps_scenario = df_scenario.iloc[0]
        print(ps_scenario)

        dict_scenario = {
            'scenario_rank' : 1,
            'scenario_proba_cond_seq' : 0.9999,
            'v_choice' : 61,
            'j_choice' : 1,
            'd_gene' : 0,
            'v_3_del': 0,
            'd_5_del': 1,
            'j_5_del': 8,
            'd_3_del': 6,
            'vd_ins' : 1,
            'vd_dinucl': [1],
            'dj_ins': 2,
            'dj_dinucl': [0, 3],
            'Mismatches' : [],
            'norm_scenario_proba_cond_seq': 0.000225
        }
        ps_scenario = pd.Series(dict_scenario)

        print("--- >ps_scenario:")
        print(ps_scenario)
        print('v_3_del', mdl.realization(ps_scenario, 'v_3_del').value)
        df_scenario_aln = mdl.get_df_scenario_aln_from_scenario(ps_scenario)
        print("--- >df_scenario_aln:")
        print(df_scenario_aln)
        print("VVVVV: ", df_scenario_aln.loc[0])
        print("len(df_scenario_aln.loc[0]['gene_segment']): ", len(df_scenario_aln.loc[0]['gene_segment']))
        print("len(df_scenario_aln.loc[0]['gene_template']): ", len(df_scenario_aln.loc[0]['gene_template']))
        # da_scenario_aln = from_df_scenario_aln_to_da_scenario_aln(df_scenario_aln)

        nrows = len(df_scenario_aln.index)
        ncols = df_scenario_aln.aln_scenario_len
        aln_scenario_np = -np.ones((nrows, ncols))

        dict_id_2_nt = Igor_dict_id_2_nt  # mdl.parms['vd_dinucl']['value'].to_dict()

        dict_nt_2_id = {v: k for k, v in dict_id_2_nt.items()}
        # FIXME: CHANGE SEGMENT TO TEMPLATE
        for ii, row in df_scenario_aln.iterrows():
            # add a map with a dictionary defined
            np_palindrome_5_end = None
            np_gene_segment = np.array(list(map(lambda nt: dict_nt_2_id[nt], list(row['gene_segment']))))
            print(type(row['gene_segment']))
            np_gene_segment = str_seq_to_np_seq(row['gene_segment'])
            np_palindrome_3_end = None

            np_gene_del_5_end = None
            np_gene_del_3_end = None

            aln_ini_gene_segment = row['offset']
            # aln_end_gene_segment = row['offset'] + len(row['gene_segment'])
            aln_end_gene_segment = row['offset'] + len(row['gene_segment'])

            if row['segment_description'] in ['V', 'D', 'J']:
                np_gene_template = str_seq_to_np_seq(row['gene_template'])
                print('segment_description', row['segment_description'])
                if row['int_gene_5_del'] < 0:
                    # Read palidrome and
                    aln_ini_gene_template = aln_ini_gene_segment # + row['int_gene_5_del'] Already included in segment
                    np_palindrome_5_end = str_seq_to_np_seq(row['palindrome_5_end'])
                    print("np_palindrome_5_end", np_palindrome_5_end)
                else:
                    aln_ini_gene_template = aln_ini_gene_template - row['int_gene_5_del']
                    np_palindrome_5_end = np.array([])

                print("ini: ", aln_ini_gene_segment, aln_ini_gene_template, row['int_gene_5_del'])

                if row['int_gene_3_del'] < 0:
                    # Read palidrome and
                    aln_end_gene_template = aln_end_gene_template # - row['int_gene_3_del']
                    np_palindrome_3_end = str_seq_to_np_seq(row['palindrome_3_end'])
                    print("np_palindrome_3_end", np_palindrome_3_end)
                else:
                    aln_end_gene_template = aln_end_gene_template + row['int_gene_3_del']
                    np_palindrome_3_end = np.array([])
                print("end: ", aln_end_gene_segment, aln_end_gene_template, row['int_gene_3_del'])

                np_gene_template = np.hstack((np_palindrome_5_end, np_gene_template))
                np_gene_template = np.hstack((np_gene_template, np_palindrome_3_end))
                print("len: ", len(np_gene_template), np_gene_template)

                aln_scenario_np[ii, aln_ini_gene_template:aln_end_gene_template] = np_gene_template

            aln_scenario_np[ii, aln_ini_gene_segment:aln_end_gene_segment] = np_gene_segment


        da_scenario_aln = xr.DataArray(aln_scenario_np, dims=('segment', 'nucleotide'),
                                       coords={
                                           "segment_description": (
                                           'segment', df_scenario_aln["segment_description"].values),
                                           "gene_description": ('segment', df_scenario_aln["gene_description"].values),
                                           "gene_template": ('segment', df_scenario_aln["gene_template"].values),
                                           "palindrome_5_end": ('segment', df_scenario_aln["palindrome_5_end"].values),
                                           "gene_ini": ('segment', df_scenario_aln["gene_ini"].values),
                                           "gene_end": ('segment', df_scenario_aln["gene_end"].values),
                                           "gene_cut": ('segment', df_scenario_aln["gene_cut"].values),
                                           "palindrome_3_end": ('segment', df_scenario_aln["palindrome_3_end"].values),
                                           "gene_segment": ('segment', df_scenario_aln["gene_segment"].values),
                                           "offset": ('segment', df_scenario_aln["offset"].values)
                                       },
                                       attrs={
                                           "anchors_CDR3": (
                                           df_scenario_aln.aln_pos_V_anchor, df_scenario_aln.aln_pos_J_anchor),
                                           "dict_id_2_nt": dict_id_2_nt
                                       }
                                       )




        # So gene segment is the one that I use
        # TODO: INSTEAD OF USING GENE_SEGMENT IS BETTER TO USE
        #  palindrome_5_end and assign colors with a pale scale, remember to use len to construct array
        #  gene_cut as usual convetion deletions -2 and color should be dashed
        #

        print("--- >da_scenario_aln:")
        print(da_scenario_aln)
        fig, ax = plot_scenario_from_da_scenario_aln(da_scenario_aln)  # , nt_lim=nt_lim, show_CDR3=show_CDR3, ax=ax)
        x = np.arange(281, 300)-0.5
        y1 = -0.5*np.ones_like(x)
        y2 = y1 + 1
        ax.fill_between(x, y1, y2=y2, hatch='///', alpha=0.0) # For positive deletions
        ax.fill_between(x+5, y1+2, y2=y2+2, hatch='/', alpha=0.1, color='black')
        # da_scenario_aln.loc[{'nucleotide': slice(*da_scenario_aln.attrs["anchors_CDR3"])}].plot()
        plt.show()

    def test_mutual_information_plot(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_scenarios(fln_scenarios)
        # da_mi = mdl.get_mutual_information()
        # print(da_mi)
        print('-'*100)
        mi = mdl.get_mutual_information_events('v_choice', 'vd_ins')
        print(mi)
        # mi = mdl.get_mutual_information_events('v_choice', 'j_choice')
        # print(mi)
        # mi = mdl.get_mutual_information_events('d_gene', 'dj_ins')
        # print(mi)
        # # da_mi_scenarios = mdl.get_mutual_information_from_df_scenarios(df_scenarios)
        # mdl.plot_mutual_information(da_mi)
        # plt.show()


    def test_Pgen_plot(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        df_generated = generate(100, mdl=mdl)
        print( df_generated )
        df_pgen = evaluate_pgen(df_generated, mdl=mdl, airr_format=False)
        print(df_pgen)




    def test_foo(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        dict_nickname_event_type = mdl.parms.get_event_dict('nickname', 'event_type')
        dict_events = {key: val for key, val in dict_nickname_event_type.items() if val != 'DinucMarkov'}
        event_lista_nicknames = list(dict_events.keys())
        for event_nickname_x, event_nickname_y in itertools.combinations_with_replacement(event_lista_nicknames, 2):
            if event_nickname_x != event_nickname_y:
                print("!=", event_nickname_x, event_nickname_y)
            else:
                print(event_nickname_x, event_nickname_y)


    def test_test(self):
        event = self.parms.get_Event(event_nickname, by_nickname=True)
        da = self.Pmarginal[event_nickname]  # real marginal DataArray

        lblEvent = event_nickname.replace("_", " ")
        xEtiqueta = lblEvent
        yEtiqueta = "P"

        if ax is None:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()

        ax.set_xlabel(xEtiqueta)
        ax.set_ylabel(yEtiqueta, rotation=0)

        if event.event_type == 'GeneChoice':
            # Bar plot
            XX = da[event_nickname].values
            YY = da.values

            ax.bar(XX, YY, **kwargs)

            lbl_XX = da['lbl__' + event_nickname].values
            ax.set_xticks(XX)
            ax.set_xticklabels(v_genLabel(lbl_XX), rotation=90)
            # return ax

        elif event.event_type == 'Insertion':
            # Use labels as a coordinate.
            # Insertions are in principle independent,
            # FIXME: but if not what to do.
            XX = da['lbl__' + event_nickname].values
            YY = da.values
            if not 'marker' in kwargs.keys():
                kwargs['marker'] = 'o'
            ax.plot(XX, YY, **kwargs)

        elif event.event_type == 'Deletion':
            # YY = self.xdata[event_nickname].values
            # XX = self.xdata[event_nickname]['lbl__' + event_nickname].values
            # ax.plot(XX, YY)
            XX = da['lbl__' + event_nickname].values
            YY = da.values
            if not 'marker' in kwargs.keys():
                kwargs['marker'] = 's'
            ax.plot(XX, YY, **kwargs)

        elif event.event_type == 'DinucMarkov':
            XX = da['x'].values
            YY = da['y'].values

            lbl__XX = da['lbl__' + 'x'].values
            lbl__YY = da['lbl__' + 'y'].values

            ZZ = da.values

            da.plot(ax=ax, x='x', y='y', vmin=0, vmax=1, cmap='gnuplot2_r', **kwargs)

            ax.set_xlabel('From')
            ax.set_xticks(XX)
            ax.set_xticklabels(lbl__XX, rotation=0)

            ax.set_ylabel('To')
            ax.set_yticks(YY)
            ax.set_yticklabels(lbl__YY, rotation=0)

            ax.set_title(lblEvent)
            ax.set_aspect('equal')
            for i, j in zip(*ZZ.nonzero()):
                ax.text(j, i, ZZ[i, j], color='white', ha='center', va='center')


        else:
            print("Event nickname " + event_nickname + " is not present in this model.")
            print("Accepted Events nicknames are : " + str(self.get_events_nicknames_list()))

        # return self.get_Event_Marginal(nickname)
        # ax.set_title(lblEvent)
        return ax


if __name__ == '__main__':
    unittest.main()
