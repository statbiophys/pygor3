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

    def test_mutual_information_plot(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_scenarios(fln_scenarios)
        da_mi = mdl.get_mutual_information()
        # print(da_mi)
        print('-'*100)
        # mi = mdl.get_mutual_information_events('v_choice', 'vd_ins')
        # print(mi)
        mi = mdl.get_mutual_information_events('d_gene', 'j_choice')
        print(mi)
        mi = mdl.get_mutual_information_events('d_gene', 'dj_ins')
        print(mi)
        # da_mi_scenarios = mdl.get_mutual_information_from_df_scenarios(df_scenarios)
        mdl.plot_mutual_information(da_mi)
        plt.show()


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
