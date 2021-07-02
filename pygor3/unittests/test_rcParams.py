import unittest
import pygor3 as p3
import pandas as pd

class MyTestCase(unittest.TestCase):
    def test_something(self):

        self.assertEqual(p3.rcParams['paths.igor_exec'], None)

        mdl_hb = p3.IgorModel.load_default('human', 'tcr_beta')
        self.assertIsInstance(mdl_hb, p3.IgorModel)
        p3.rcParams['paths.igor_exec'] = 'mygor'
        df_seqs = p3.generate(10, mdl_hb)

        self.assertIsInstance(df_seqs, pd.DataFrame)




if __name__ == '__main__':
    unittest.main()
