import unittest

from pygor3 import IgorModel
from pygor3 import get_default_IgorModel
from pygor3 import get_df_cross_entropy
from pygor3 import mean_IgorModel
import numpy as np
import pandas as pd
import pandas.testing as pd_testing

class MyTestCase(unittest.TestCase):
    def test_cross_entropy_definition(self):
        mdl_hb = get_default_IgorModel('human', 'beta')
        mdl_0 = IgorModel.make_default_from_Dataframe_dict(mdl_hb.genomic_dataframe_dict)
        P = mdl_hb['v_choice']
        Q = mdl_0['v_choice']
        print(P)
        print(Q)
        H_v_choice = mdl_hb.get_entropy_event('v_choice')
        print("H_v_choice: ", H_v_choice)

        # Cross_entropy H(P,Q) = \sum_x P(x) \log(P(x)/Q(x))
        CE = -( P*np.log2(Q )).sum()
        print("CE: ",CE)

        # DKL(P||Q)
        DKL = CE - H_v_choice
        print("DKL: ", DKL)
        self.assertIsInstance(mdl_hb, IgorModel)
        self.assertIsInstance(mdl_0, IgorModel)

    def test_cross_entropy_definition(self):
        mdl_hb = get_default_IgorModel('human', 'beta')
        mdl_0 = IgorModel.make_default_from_Dataframe_dict(mdl_hb.genomic_dataframe_dict)

        # entropy hb
        df_entropy_hb = mdl_hb.get_df_entropy_decomposition()
        print(df_entropy_hb)
        print(df_entropy_hb['entropy'].sum())

        # entropy 0
        df_entropy_0 = mdl_0.get_df_entropy_decomposition()
        print(df_entropy_0)
        print(df_entropy_0['entropy'].sum())

        # cross entropy (mdl_hb, mdl_hb)
        df_cross_entropy_same_hb = get_df_cross_entropy(mdl_hb, mdl_hb)
        print(df_cross_entropy_same_hb)
        print(df_cross_entropy_same_hb['entropy'].sum())
        self.assertIsInstance(df_cross_entropy_same_hb, pd.DataFrame)

        self.assertTrue(df_entropy_hb.equals(df_cross_entropy_same_hb))

        print("-"*50, "mdl_0, mdl_hb")
        df_cross_entropy_diff_1 = get_df_cross_entropy(mdl_0, mdl_hb)
        print(df_cross_entropy_diff_1)
        print(df_cross_entropy_diff_1['entropy'].sum())
        df_DKL_divergence = df_cross_entropy_diff_1['entropy'] - df_entropy_0['entropy']
        print("DKL: ", df_DKL_divergence)
        print("DKL.sum: ", df_DKL_divergence.sum())

        print("-" * 50, "mdl_hb, mdl_0")
        df_cross_entropy_diff_2 = get_df_cross_entropy(mdl_hb, mdl_0)
        print(df_cross_entropy_diff_2)
        print(df_cross_entropy_diff_2['entropy'].sum())
        df_DKL_divergence = df_cross_entropy_diff_2['entropy'] - df_entropy_hb['entropy']
        print("DKL: ", df_DKL_divergence)
        print("DKL.sum: ", df_DKL_divergence.sum())

    def test_JSD(self):
        # FIXME: IN DEV
        mdl_hb = get_default_IgorModel('human', 'beta')
        mdl_0 = IgorModel.make_default_from_Dataframe_dict(mdl_hb.genomic_dataframe_dict)
        print( 0.5*(mdl_hb['v_choice'] + mdl_0['v_choice']).sum() )

        print(mdl_hb['j_choice'].sum(dim='j_choice'))
        print( 0.5*(mdl_hb['j_choice'] + mdl_0['j_choice']) )
        # print(0.5 * (mdl_hb['j_choice'] + mdl_0['j_choice']).sum())

        # mean_IgorModel(mdl_hb, mdl_0)


if __name__ == '__main__':
    unittest.main()
