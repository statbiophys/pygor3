import unittest
from pygor3 import *
import subprocess


class MyTestCase(unittest.TestCase):
    def test_something(self):
        task = IgorTask.default_model("human", "beta")
        task._run_generate(10)

        # task.igor_fln_generated_seqs_werr
        # task.igor_fln_generated_realizations_werr
        # task.igor_fln_generation_info
        # cmd = "cat " + task.igor_fln_generated_seqs_werr
        # print(cmd)
        # subprocess.run(cmd, shell=True)
        # cmd = "cat " + task.igor_fln_generated_realizations_werr
        # print(cmd)
        # subprocess.run(cmd, shell=True)

        pd_realizations = get_dataframe_from_fln_generated_seqs_werr(task.igor_fln_generated_realizations_werr)
        print(pd_realizations)

        sce = IgorScenario()
        print(sce.to_dict())
        self.assertEqual(True, True)

    def test_IgorModel(self):
        task = IgorTask.default_model("human", "beta", igor_wd='igor_temporal')
        task._run_generate(5, seed=10)
        pd_realizations = task.mdl.get_dataframe_from_fln_generated_realizations_werr(
            task.igor_fln_generated_realizations_werr)
        print(pd_realizations)

        # task.load_db_from_indexed_sequences()
        # task.load_IgorModel()
        # task.load_db_from_genomes()
        # task.load_db_from_models()
        # task.load_db_from_bestscenarios()

        task.evaluate(pd_realizations, N_scenarios=3, clean_batch=False)
        # task._run_generate(10, seed=20)
        pd_realizations2 = task.mdl.get_dataframe_from_fln_generated_realizations_werr(
            task.igor_fln_generated_realizations_werr)
        print(pd_realizations2)
        # subprocess.run("rm -r igor_temporal", shell=True)





if __name__ == '__main__':
    unittest.main()
