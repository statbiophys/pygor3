import unittest
import numpy as np
from pygor3 import IgorRec_Event, IgorModel, IgorEvent_realization
import Bio.Seq
import pandas as pd
from pygor3.utils import dna_complementary
import matplotlib.pyplot as plt
import numba as nb


class MyTestCase(unittest.TestCase):


    def test_realization(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_scenarios(fln_scenarios)
        v_choice = mdl.realization(df_scenarios, 'v_choice')
        print(v_choice.id)
        print(v_choice.value)
        print(type(v_choice), len(v_choice))
        print(type(v_choice.id), len(v_choice.id))
        print(type(v_choice.value), len(v_choice.value))
        print("*"*50)
        try:
            vd_dinucl = mdl.realization(df_scenarios, 'vd_dinucl')
            print(vd_dinucl.id)
            print(vd_dinucl.value)
            print(type(vd_dinucl))
            print(type(vd_dinucl.id))
            print(type(vd_dinucl.value))
        except Exception as e:
            raise e

    def test_get_df_realizations_dinucl(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_scenarios(fln_scenarios)
        df_realizations = mdl.get_df_realizations_dinucl(df_scenarios, 'vd_dinucl')
        vd_seq = df_realizations.value.apply(lambda x: "".join(x))
        print(vd_seq.value_counts())


    def test_get_CDR3_from_scenarios(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_scenarios(fln_scenarios)

        def observable_get_CDR3(ps_scenario):
            """
            Return the numbers of amino acids in vd insertions
            """
            try:
                v_choice = mdl.realization(ps_scenario, 'v_choice')
                j_choice = mdl.realization(ps_scenario, 'j_choice')
                d_gene = mdl.realization(ps_scenario, 'd_gene')

                v_3_del = mdl.realization(ps_scenario, 'v_3_del')
                d_5_del = mdl.realization(ps_scenario, 'd_5_del')
                d_3_del = mdl.realization(ps_scenario, 'd_3_del')
                j_5_del = mdl.realization(ps_scenario, 'j_5_del')

                # vd_ins = mdl.realization(ps_scenario, 'vd_ins')
                vd_dinucl = mdl.realization(ps_scenario, 'vd_dinucl')

                # dj_ins = mdl.realization(ps_scenario, 'dj_ins')
                dj_dinucl = mdl.realization(ps_scenario, 'dj_dinucl')

                v_anchor = mdl.V_anchor(v_choice.id)
                j_anchor = mdl.J_anchor(j_choice.id)

                # TODO: mdl.get_CDR3_seq(ps_scenario)
                ##### V_Gene
                v_gene_len = len(v_choice.value)
                # mdl.get_CDR3_seq(ps_scenario)
                v_ini = 0
                v_end = v_gene_len
                str_v_3_palidrome = ""
                if v_3_del.value < 0:
                    str_v_3_palidrome = dna_complementary((v_choice.value[v_3_del.value:])[::-1])
                else:
                    v_end = v_end - v_3_del.value

                str_V_segment = v_choice.value[v_ini:v_end] + str_v_3_palidrome

                ##### D_gene
                d_gene_len = len(d_gene.value)
                d_ini = 0
                d_end = d_gene_len
                str_d_5_palidrome = ""
                if d_5_del.value < 0:
                    int_ini = 0
                    str_d_5_palidrome = dna_complementary((d_gene.value[:-d_5_del.value])[::-1])
                else:
                    d_ini = d_5_del.value

                str_d_3_palidrome = ""
                if d_3_del.value < 0:
                    str_d_3_palidrome = dna_complementary((d_gene.value[d_3_del.value:])[::-1])
                else:
                    d_end = d_end - d_3_del.value

                str_D_segment = str_d_5_palidrome + d_gene.value[d_ini:d_end] + str_d_3_palidrome

                ##### J_gene
                j_gene_len = len(j_choice.value)
                j_ini = 0
                j_end = j_gene_len
                str_j_5_palindrome = ""
                if j_5_del.value < 0:
                    j_ini = 0
                    str_j_5_palindrome = dna_complementary((j_choice.value[:-j_5_del.value])[::-1])
                else:
                    j_ini = j_5_del.value

                str_J_segment = str_j_5_palindrome + j_choice.value[j_ini:j_end]

                str_VD_segment = "".join(vd_dinucl.value)
                str_DJ_segment = "".join(dj_dinucl.value[::-1])

                if (v_anchor > v_end) or (j_anchor < j_ini):
                    return np.NaN
                else:
                    str_sequence = str_V_segment[v_anchor:] + str_VD_segment + str_D_segment + str_DJ_segment + str_J_segment[:j_anchor]
                    if len(str_sequence) % 3 == 0:
                        return str(Bio.Seq.Seq(str_sequence).translate())
                    else:
                        return np.NaN
            except Exception as e:
                print(e)
                return None # FIXME: what to do in these cases, this will change the normalization.

        df_scenarios['CDR3'] = mdl.get_observable_from_df_scenarios(observable_get_CDR3, df_scenarios)


        print(df_scenarios['CDR3'])
        aver = df_scenarios.groupby(['v_choice', 'CDR3'])['norm_scenario_proba_cond_seq'].apply(lambda x: x.sum())
        # aver = df_scenarios.groupby(['CDR3'])['norm_scenario_proba_cond_seq'].apply(lambda x: x.sum())
        print(aver)
        da = aver.to_xarray()
        print(da.values)
        da.values = np.nan_to_num(da.values, 0)
        print("da.sum(): ", da.sum())
        print(aver)
        print(da)
        print(da.dims)
        print(da.coords)
        print(da[{'v_choice':1 }])
        print(da.sel({'v_choice': 1}))
        print('*'*50)
        print(da.sel({'CDR3': 'CACIVTRRDKNCFLAV'}))

        mdl.get_probability_matrix_from_event_list_and_scenarios_dataframe()
        mdl.get_P_from_scenarios_cols()
        mdl.get_CDR3_from_scenario()

        # da.sel({'CDR3': 'CACIVTRRDKNCFLAV'}).plot()
        # plt.show()

        # aaa = df_scenarios.groupby('CDR3')['norm_scenario_proba_cond_seq']
        # print(aaa)
        # df_scenarios['norm_scenario_proba_cond_seq'] = get_df_normalize_prob(df_scenarios)
        # mdl.get_scenarios_normalization()




        df_scenarios['CDR3'] = mdl.get_observable_from_df_scenarios(observable_get_CDR3, df_scenarios)


        print(df_scenarios['CDR3'])
        aver = df_scenarios.groupby(['v_choice', 'CDR3'])['norm_scenario_proba_cond_seq'].apply(lambda x: x.sum())
        # aver = df_scenarios.groupby(['CDR3'])['norm_scenario_proba_cond_seq'].apply(lambda x: x.sum())
        print(aver)
        da = aver.to_xarray()
        print(da.values)
        da.values = np.nan_to_num(da.values, 0)
        print("da.sum(): ", da.sum())
        print(aver)
        print(da)
        print(da.dims)
        print(da.coords)
        print(da[{'v_choice':1 }])
        print(da.sel({'v_choice': 1}))
        print('*'*50)
        print(da.sel({'CDR3': 'CACIVTRRDKNCFLAV'}))

    def test_get_CDR3_from_df_scenarios(self):
        def observable_get_CDR3(v_choice, j_choice, d_gene, v_3_del, d_5_del, d_3_del, j_5_del, vd_dinucl, dj_dinucl):
            """
            Return the numbers of amino acids in vd insertions
            """
            try:
                # v_choice = mdl.realization(ps_scenario, 'v_choice')
                # j_choice = mdl.realization(ps_scenario, 'j_choice')
                # d_gene = mdl.realization(ps_scenario, 'd_gene')
                #
                # v_3_del = mdl.realization(ps_scenario, 'v_3_del')
                # d_5_del = mdl.realization(ps_scenario, 'd_5_del')
                # d_3_del = mdl.realization(ps_scenario, 'd_3_del')
                # j_5_del = mdl.realization(ps_scenario, 'j_5_del')
                #
                # # vd_ins = mdl.realization(ps_scenario, 'vd_ins')
                # vd_dinucl = mdl.realization(ps_scenario, 'vd_dinucl')
                #
                # # dj_ins = mdl.realization(ps_scenario, 'dj_ins')
                # dj_dinucl = mdl.realization(ps_scenario, 'dj_dinucl')

                v_anchor = mdl.V_anchor(v_choice.id)
                j_anchor = mdl.J_anchor(j_choice.id)

                # TODO: mdl.get_CDR3_seq(ps_scenario)
                ##### V_Gene
                v_gene_len = len(v_choice.value)
                # mdl.get_CDR3_seq(ps_scenario)
                v_ini = 0
                v_end = v_gene_len
                str_v_3_palidrome = ""
                if v_3_del.value < 0:
                    str_v_3_palidrome = dna_complementary((v_choice.value[v_3_del.value:])[::-1])
                else:
                    v_end = v_end - v_3_del.value

                str_V_segment = v_choice.value[v_ini:v_end] + str_v_3_palidrome

                ##### D_gene
                d_gene_len = len(d_gene.value)
                d_ini = 0
                d_end = d_gene_len
                str_d_5_palidrome = ""
                if d_5_del.value < 0:
                    int_ini = 0
                    str_d_5_palidrome = dna_complementary((d_gene.value[:-d_5_del.value])[::-1])
                else:
                    d_ini = d_5_del.value

                str_d_3_palidrome = ""
                if d_3_del.value < 0:
                    str_d_3_palidrome = dna_complementary((d_gene.value[d_3_del.value:])[::-1])
                else:
                    d_end = d_end - d_3_del.value

                str_D_segment = str_d_5_palidrome + d_gene.value[d_ini:d_end] + str_d_3_palidrome

                ##### J_gene
                j_gene_len = len(j_choice.value)
                j_ini = 0
                j_end = j_gene_len
                str_j_5_palindrome = ""
                if j_5_del.value < 0:
                    j_ini = 0
                    str_j_5_palindrome = dna_complementary((j_choice.value[:-j_5_del.value])[::-1])
                else:
                    j_ini = j_5_del.value

                str_J_segment = str_j_5_palindrome + j_choice.value[j_ini:j_end]

                str_VD_segment = "".join(vd_dinucl.value)
                str_DJ_segment = "".join(dj_dinucl.value[::-1])

                if (v_anchor > v_end) or (j_anchor < j_ini):
                    return np.NaN
                else:
                    str_sequence = str_V_segment[v_anchor:] + str_VD_segment + str_D_segment + str_DJ_segment + str_J_segment[:j_anchor+3]
                    if len(str_sequence % 3 == 0 ):
                        return str(Bio.Seq.Seq(str_sequence).translate())
                    else:
                        return np.NaN
            except Exception as e:
                print(e)
                return None # FIXME: what to do in these cases, this will change the normalization.



    def test_pd_realizations(self):
        # FIXME: IN DEV
        mdl = IgorModel.load_default("human", "tcr_beta")
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_from_fln_generated_realizations_werr(fln_scenarios)

        ps_scenario = df_scenarios.iloc[2]
        print("___vd_dinucl___")
        realiz = mdl.realization(ps_scenario, 'vd_dinucl')
        # aaa = observable_n_Cys_in_vd(ps_scenario)
        print(realiz.to_dict())
        print(realiz.id)
        print("".join(realiz.value))
        print(realiz.name)
        print(type(realiz))

        print("___v_choice___")
        realiz = mdl.realization(ps_scenario, 'v_choice')
        print(realiz.to_dict())
        print(realiz.id)
        print("".join(realiz.value))
        print(realiz.name)
        print(type(realiz))

        # df_scenarios['2ndCys'] = mdl.get_observable_from_scenarios_dataframe(observable_n_Cys_in_vd, df_scenarios)
        # print(df_scenarios)

    def test_observable_something(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_from_fln_generated_realizations_werr(fln_scenarios)

        def observable_anchors_aa(ps_scenario):
            """
            Return the numbers of amino acids in vd insertions
            """
            try:
                v_choice = mdl.realization(ps_scenario, 'v_choice')
                j_choice = mdl.realization(ps_scenario, 'j_choice')
                # v_3_del = mdl.realization(ps_scenario, 'v_3_del')
                # vd_ins = mdl.realization(ps_scenario, 'vd_ins')
                # vd_dinucl = mdl.realization(ps_scenario, 'vd_dinucl')

                v_anchor = mdl.V_anchor(v_choice.id)
                j_anchor = mdl.J_anchor(j_choice.id)

                # anchor_aa = v_choice.value[v_anchor:v_anchor+2]
                str_V_anchor_aa = str(Bio.Seq.Seq(v_choice.value[v_anchor:v_anchor + 3]).translate())
                str_J_anchor_aa = str(Bio.Seq.Seq(j_choice.value[j_anchor:j_anchor + 3]).translate())
                if str_J_anchor_aa == 'F':
                    return np.NaN
                else:
                    return str_J_anchor_aa
                # return str_V_anchor_aa + " " + str_J_anchor_aa
            except Exception as e:
                return None # FIXME: what to do in these cases, this will change the normalization.
            # v_segment_length = len(v_choice.value) - v_3_del.value
            # n_extra_nt = v_segment_length % 3
            # n_nt = (3 - n_extra_nt) % 3
            # if (vd_ins < n_nt):
            #     return 0
            # else:
            #     return (vd_ins - n_nt) // 3


        df_scenarios['tmp_aa'] = mdl.get_observable_from_df_scenarios(observable_anchors_aa, df_scenarios)

        no_F = df_scenarios['tmp_aa']
        df_scenarios.groupby('tmp_aa')
        print(no_F)
        print(len(df_scenarios))

    def test_numpy_recarray(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        np_v_choice = mdl.parms.Event_dict['v_choice'].to_records()
        realiz = np_v_choice[0]
        print(realiz.dtype)
        print(type(np_v_choice[0]))

        np_v_3_del = mdl.parms.Event_dict['v_3_del'].to_records()
        realiz = np_v_3_del[0]
        print(realiz.dtype)
        print(type(np_v_3_del[0]))

    def test_numpy_structured_array(self):
        # dtipo = [('id', '<i8'), ('value', 'O'), ('name', 'O')]
        dtipo = [('id', int), ('value', object), ('name', object)]
        np_realization = np.rec.array(
            [(0, "ACCTAGGATC", "TRAV_GENE_NAME_7"), (1, "CATTATAGGAT", "TRAV_GENE_NAME_5")],
            dtype=[('id', int), ('value', object), ('name', object)]
        )
        np_realization2 = np.rec.array(
            [(0, "ACCTAGGATC", "TRAV_GENE_NAME_7")],
            dtype=[('id', int), ('value', object), ('name', object)]
        )
        print(np_realization2)
        print(len(np_realization2.value))
        zzz = np.array(5)
        print(zzz)
        print(type(zzz))


        mdl = IgorModel.load_default("human", "tcr_beta")
        id = [0, 2, 1, 0, 0]
        # id = 1
        try:
            aver = mdl.parms.Event_dict['vd_dinucl'].loc[id]
            realiz = IgorEvent_realization.from_tuple(np.array(id), np.array(aver.value), np.array(aver.name))
        except Exception as e:
            print(e)

        print("realiz.value: ", realiz.value)
        print("type(realiz.value): ", type(realiz.value))

        print("realiz.id: ", realiz.id)
        print("type(realiz.id): ", type(realiz.id))




    def test_something(self):
        mdl = IgorModel.load_default("human", "tcr_beta")
        # v_choice = mdl.parms.Event_dict['v_choice'].loc[0]
        fln_scenarios = "delete_me/ttmmpp_output/best_scenarios_counts.csv"
        df_scenarios = mdl.get_dataframe_from_fln_generated_realizations_werr(fln_scenarios)
        print(df_scenarios)
        ps_scenario = df_scenarios.loc[0]
        print('^'*50)
        id = 0
        ps_realiz = mdl.parms.Event_dict['v_choice'].loc[id]
        print("ps_realiz: ", type(ps_realiz))
        # np_realiz = np.recarray( dtype=[('id', ''), ('value', ''), ('name', '')])
        v_choice = IgorEvent_realization.from_tuple(id, ps_realiz['value'], ps_realiz['name'])
        print(v_choice)

        print('=' * 50)
        id = [1, 2, 0, 0]
        df_realiz = mdl.parms.Event_dict['vd_dinucl'].loc[id]
        print("df_realiz: ", df_realiz)
        print("type(df_realiz): ", type(df_realiz))




        # vd_dinucl = IgorEvent_realization.from_pandas(df_realiz)
        # print(vd_dinucl)
        # print(vd_dinucl[0].value)

        # TODO if id is a list how can I use the same function



        # mdl.parms.Event_dict['v_choice'].loc[id]


        # aaa = mdl.pd_realization(ps_scenario, 'v_choice')
        # print(aaa.to_dict())

        """
        id = 0
        ps_realiz = mdl.parms.Event_dict['v_choice'].loc[id]
        print("ps_realiz: ", type(ps_realiz))
        v_choice = IgorEvent_realization.from_tuple(id, ps_realiz['value'], ps_realiz['name'])

        print(v_choice.id)
        print(v_choice.value)
        print(v_choice.name)
        
        id = [0, 3, 1, 0, 2]
        ps_realiz = mdl.parms.Event_dict['vd_dinucl'].loc[id]
        ddd = ps_realiz['value'].values
        aaa = ''.join(ddd)
        print(aaa)
        # v_choice = IgorEvent_realization.from_tuple(id, ps_realiz['value'], ps_realiz['name'])

        # evento = IgorRec_Event()
        # self.assertIsInstance(evento, IgorRec_Event)
        # self.assertEqual(True, False
        """


if __name__ == '__main__':
    unittest.main()
