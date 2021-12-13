import unittest
from pygor3 import IgorTask
from pygor3 import IgorModel
from pygor3 import generate
from pygor3 import naive_align
import pandas as pd


class MyTestCase(unittest.TestCase):
    def test_something(self):
        mdl_hb = IgorModel.load_default('human', 'tcr_beta')
        df_seqs = generate(3, mdl_hb)
        self.assertIsInstance(df_seqs, pd.DataFrame)
        task = IgorTask(mdl=mdl_hb)
        print(task)
        self.assertIsInstance(task, IgorTask)
        aaa = task.align(df_seqs, mdl=mdl_hb)
        print(aaa['V'])
        print(aaa['V'].columns)
        print('-' * 50)
        bbb = aaa['V'].groupby(by=['seq_index']).apply(lambda g: g.sortby('score', ascending=False))
        print(bbb)

    def test_something02(self):
        mdl_hb = IgorModel.load_default('human', 'tcr_beta')
        df_seqs = generate(3, mdl_hb)
        print(df_seqs)
        self.assertIsInstance(df_seqs, pd.DataFrame)
        df_output = naive_align(df_seqs, mdl=mdl_hb)
        print(df_output[0])
        print(df_output[1])
        print(df_output[1]['v_anchor'])
        print(df_output[1]['seq_index'])

        # self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
