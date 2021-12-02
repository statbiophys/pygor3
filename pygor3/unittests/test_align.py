import unittest
from pygor3 import IgorTask
from pygor3 import IgorModel
from pygor3 import generate
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
        print('-'*50)
        bbb = aaa['V'].groupby(by=['seq_index']).apply(lambda g: g.sortby('score', ascending=False))
        print(bbb)

        # self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
