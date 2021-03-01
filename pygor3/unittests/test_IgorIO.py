#!/usr/bin/env python3
import unittest

from pygor3 import *
from pygor3.imgt import *
import os
import subprocess

class TestPygor3(unittest.TestCase):
    null_genome = IgorRefGenome()
    null_mdl = IgorModel()
    
    # def test_1(self):
    #     mdl = IgorModel()
    #     mdl = self.null_mdl
    #
    #     self.assertEqual(self.null_mdl, mdl)
    #
    # def test_IgorRefGenome(self):
    #     #hs_genomes = IgorRefGenome.load_VDJ_from_IMGT_website("Homo+sapiens", "TRB")
    #     self.assertEqual(self.null_genome.df_genomicVs, None)
    #     self.assertEqual(self.null_genome.df_genomicDs, None)
    #     self.assertEqual(self.null_genome.df_genomicJs, None)
    #     # self.assertEqual(type(hs_genomes), type(null_genomes))
    #
    # def test_IgorRefGenome_download_imgt_default_path(self):
    #     # 1 Download using imgt
    #     df_genomics_dict = download_ref_genome_VDJ("Homo+sapiens", "TRB")
    #     self.assertTrue(os.path.exists("models"))
    #     genome = IgorRefGenome.load_from_path("models/Homo+sapiens/TRB/ref_genome")
    #     self.assertEqual(type(self.null_genome), type(genome))
    #     subprocess.call("rm -r models", shell=True)
    #     self.assertTrue(not os.path.exists("models"))
    #
    # def test_pygor_imgt_ref_genome(self):
    #     my_ref_genome = "my_ref_genome"
    #     df_genomics_dict = download_ref_genome_VDJ("Homo+sapiens", "TRB", ref_genes_path=my_ref_genome)
    #     self.assertTrue(os.path.exists(my_ref_genome))
    #     self.assertTrue(os.path.exists(my_ref_genome + "/genomicVs.fasta"))
    #     self.assertTrue(os.path.exists(my_ref_genome + "/genomicDs.fasta"))
    #     self.assertTrue(os.path.exists(my_ref_genome + "/genomicJs.fasta"))
    #     self.assertTrue(os.path.exists(my_ref_genome + "/J_gene_CDR3_anchors.csv"))
    #     self.assertTrue(os.path.exists(my_ref_genome + "/V_gene_CDR3_anchors.csv"))
    #
    #     # self.assertTrue(os.path.exists("my_ref_genome/genomicJs.fasta"))
    #     subprocess.call("rm -r " + my_ref_genome, shell=True)
    #     self.assertTrue(not os.path.exists(my_ref_genome))
    #
    # def test_pygor_imgt_modelspath(self):
    #     modelspath = "my_models"
    #     imgt_species = "Homo+sapiens"
    #     imgt_chain = "TRB"
    #     df_genomics_dict = download_ref_genome_VDJ(imgt_species, imgt_chain, modelspath=modelspath)
    #     self.assertTrue(os.path.exists(modelspath))
    #     self.assertEqual(type(df_genomics_dict['V']), type(pd.DataFrame()) )
    #     # self.assertTrue(os.path.exists(modelspath + "/genomicVs.fasta"))
    #     # self.assertTrue(os.path.exists(modelspath + "/genomicDs.fasta"))
    #     # self.assertTrue(os.path.exists(modelspath + "/genomicJs.fasta"))
    #     # self.assertTrue(os.path.exists(modelspath + "/J_gene_CDR3_anchors.csv"))
    #     # self.assertTrue(os.path.exists(modelspath + "/V_gene_CDR3_anchors.csv"))
    #
    #     # self.assertTrue(os.path.exists("my_ref_genome/genomicJs.fasta"))
    #     subprocess.call("rm -r " + modelspath, shell=True)
    #     self.assertTrue(not os.path.exists(modelspath))


    def test_IgorRefGenome_download_imgt_in_path(self):
        # 1 Download using imgt
        modelspath = "my_models"
        imgt_species = "Homo+sapiens"
        imgt_chain = "TRB"
        df_genomics_dict = download_ref_genome_VDJ(imgt_species, imgt_chain, modelspath=modelspath)
        self.assertTrue(os.path.exists(modelspath))
        genome_from_df = IgorRefGenome.load_from_dataframe_genomics_dict(df_genomics_dict)

        # print(df_genomics_dict['V'])

        df_genomics_dict['V'].to_csv("tmp_V.csv", sep=';')
        print(df_genomics_dict['V'].columns)

        print(genome_from_df.df_V_ref_genome.columns)
        task = IgorTask.default_model("mouse", "beta")

        task.genomes.df_V_ref_genome.columns
        task.genomes.df_genomicVs
        task.igor_path_ref_genome
        task.load_IgorRefGenome()

        IgorTask.load_from_batchname()


        genome_from_path = IgorRefGenome.load_from_path(modelspath+"/"+imgt_species+"/"+imgt_chain+"/ref_genome")
        self.assertEqual(len(genome_from_df.df_V_ref_genome), len(genome_from_path.df_V_ref_genome))
        """
        print(genome_from_path.df_V_ref_genome)

        subprocess.call("rm -r " + modelspath, shell=True)
        self.assertTrue(not os.path.exists(modelspath))
        # genome = IgorRefGenome.load_VDJ_from_IMGT_website()

        # Now I want to create a temporary directory

        # self.assertEqual(type(df_genomics_dict), type(dict()))
        """




if __name__ == '__main__':
    unittest.main()

