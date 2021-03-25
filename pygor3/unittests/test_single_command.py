import unittest
import pandas as pd
from pygor3 import IgorModel
from pygor3 import generate
from pygor3 import infer
from pygor3 import evaluate



class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.mdl = IgorModel.load_default("human", "tcr_beta")

    def test_generate(self):
        df = generate(10, mdl=self.mdl)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertTrue(len(df) == 10)

    def test_infer(self):
        self.assertEqual(True, False)


if __name__ == '__main__':
    unittest.main()
