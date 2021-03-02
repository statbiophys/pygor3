import unittest
from pygor3 import IgorTask
from pygor3 import IgorModel
from pygor3 import IgorRefGenome
from pygor3 import generate
from pygor3 import infer
import pandas as pd
from pygor3 import v_genLabel

class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.species = "human"
        self.chain = "tcr_beta"
        self.mdl = IgorModel.load_default(self.species, self.chain)
        # mdl = IgorTask.run_generate(return_df=True)
        self.pd_sequences = generate(10, self.mdl)


    def test_something(self):
        # print(self.pd_sequences)
        self.assertIsInstance(self.pd_sequences, pd.DataFrame)
        # ref_genome_default = IgorRefGenome.load_default(self.species, self.chain)
        ref_genome_default = self.mdl.parms.get_IgorRefGenome()
        ref_genome_default.clean_empty_anchors()
        # From IgorRefGenome to IgorModel
        print(ref_genome_default.V)
        print(ref_genome_default.D)
        print(ref_genome_default.J)
        # ref_genome_default.df_V_ref_genome['name'] = v_genLabel(ref_genome_default.df_V_ref_genome['name'])
        # ref_genome_default.df_genomicDs['name'] = v_genLabel(ref_genome_default.df_genomicDs['name'])
        # ref_genome_default.df_J_ref_genome['name'] = v_genLabel(ref_genome_default.df_J_ref_genome['name'])
        # ref_genome_default.clean_empty_anchors()
        # ref_genome.clean_empty_anchors()

        """
        mdl_ini = IgorModel.make_default_model_from_genomes(ref_genome_default)

        # mdl_ini = IgorModel.make_default_VDJ(ref_genome_default.df_V_ref_genome,
        #                                      ref_genome_default.df_genomicDs,
        #                                      ref_genome_default.df_J_ref_genome)

        # print(self.mdl['d_gene'])
        mdl_infer = infer(self.pd_sequences, mdl=mdl_ini)
        print("="*30, "inferred")
        print(mdl_infer)
        # self.assertIsInstance(mdl_infer)
        """


if __name__ == '__main__':
    unittest.main()
