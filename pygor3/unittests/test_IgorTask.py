import unittest
from pygor3 import *
import subprocess

class MyTestCase(unittest.TestCase):
    null_task = IgorTask()
    null_task.update_ref_genome()
    null_task.update_model_filenames()
    null_task.update_batch_filenames()
    null_task.update_batchname()

    null_genomes = IgorRefGenome()
    null_genomes.update_fln_names()

    null_task.igor_model_dir_path

    def test_generate_sequences(self):
        default_task = IgorTask.default_model("human", "beta")
        sequences_pd = default_task.run_generate(10)
        sequences_np = sequences_pd["nt_sequence"].values
        default_task.run_evaluate()
        fln_txt = "random.txt"
        np.savetxt(fln_txt, sequences_np, fmt="%s")
        fln_csv = "random.csv"
        sequences_pd.to_csv(fln_csv, sep=';')
        self.assertEqual(True, True)

    # def test_load_default(self):
    #     """test default task"""
    #     default_task = IgorTask.default_model("human", "beta")
    #     self.assertEqual(type(self.null_task), type(default_task))
    #
    #     # Load from default location explicitaly
    #     path_genomes = IgorRefGenome.load_from_path(default_task.igor_path_ref_genome)
    #     tmp_flag = default_task.genomes.df_genomicVs.equals(path_genomes.df_genomicVs)
    #     self.assertTrue(tmp_flag)
    #     # Now write it in a database
    #     default_task.create_db("hb_genome.db")
    #     print(default_task.igor_fln_db)
    #     default_task.load_db_from_genomes()
    #     default_task.db_export_IgorGenomes("export_genome_from_db")
    #
    #     fromdb_genomes = IgorRefGenome.load_from_path("export_genome_from_db")
    #     tmp_flag = default_task.genomes.df_genomicVs.equals(fromdb_genomes.df_genomicVs)
    #     self.assertTrue(tmp_flag)
    #
    #
    #     subprocess.call("rm -r hb_genome.db", shell=True)
    #     subprocess.call("rm -r export_genome_from_db", shell=True)````


if __name__ == '__main__':
    unittest.main()
