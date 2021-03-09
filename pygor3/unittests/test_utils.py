import unittest
# from pygor3 import IgorTask
from pygor3.utils import *
import subprocess
import os

class MyTestCase(unittest.TestCase):

    def test_run_igor_datadir(self):
        igor_exec_path = run_get_igor_exec_path()
        print(igor_exec_path)
        self.assertTrue(os.path.isfile(igor_exec_path))
        igor_datadir = run_get_igor_datadir()
        print(igor_datadir)
        # self.assertEqual(True, True)
        self.assertTrue(os.path.isdir(igor_datadir))

    """
    def test_something(self):
        dicto = get_default_ref_genomes_species_chain("human", "tcr_beta")
        print(dicto)
        task = IgorTask.default_model("human", "beta")
        sequences = task.run_generate(10)
        self.assertEqual(len(sequences), 10)
        task.run_clean_batch()
        sequences_np = sequences['nt_sequence'].values
        # print(sequences)

        write_sequences_to_file(sequences_np, 'sequences.csv')
        # subprocess.call("cat sequences.csv", shell=True)
        task.run_infer(igor_read_seqs='sequences.csv')
        self.assertEqual(True, True)
    """




if __name__ == '__main__':
    unittest.main()
