import unittest
from pygor3 import IgorModel
from pygor3 import generate
from pygor3 import infer
from pygor3 import evaluate
import pandas as pd
from pygor3 import v_genLabel



class MyTestCase(unittest.TestCase):
    def test_something(self):
        mdl_hb = IgorModel.load_default('human', 'beta')
        df_seqs = generate(10, mdl_hb)
        df_airr_sequences = evaluate(df_seqs['nt_sequence'], mdl_hb, N_scenarios=5, use_db=True, airr_format=True)
        print(df_airr_sequences)
        self.assertIsInstance(df_airr_sequences, pd.DataFrame)


if __name__ == '__main__':
    unittest.main()
