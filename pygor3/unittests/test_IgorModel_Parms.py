import unittest
import tempfile
from pygor3.utils import *
from pygor3 import IgorModel_Parms
from pygor3 import IgorRec_Event
from pygor3 import IgorRefGenome
import os

import time


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_dir = tempfile.TemporaryDirectory(dir=".", prefix="models")
        self.fln_dict = get_default_fln_names_for_model_dir(self.tmp_dir.name)

    def test_IgorModel_Parms(self):
        species = "human"
        chain = "tcr_beta"

        mdl_parms = IgorModel_Parms()

        fln_model_parms, fln_model_marginals = get_default_models_paths_species_chain(species, chain)
        self.assertTrue(os.path.isfile(fln_model_parms))
        mdl_parms.read_model_parms(fln_model_parms)

        self.assertFalse('anchor_index' in mdl_parms.df_V_ref_genome.columns.to_list())
        self.assertFalse('anchor_index' in mdl_parms.df_J_ref_genome.columns.to_list())

        fln_dict = get_default_fln_dict_ref_genomes_species_chain(species, chain)
        self.assertTrue(os.path.isfile(fln_dict["fln_V_gene_CDR3_anchors"]))
        self.assertTrue(os.path.isfile(fln_dict["fln_J_gene_CDR3_anchors"]))

        mdl_parms.attach_V_anchors_from_file(fln_dict["fln_V_gene_CDR3_anchors"])
        mdl_parms.attach_J_anchors_from_file(fln_dict["fln_J_gene_CDR3_anchors"])

        self.assertTrue('anchor_index' in mdl_parms.df_V_ref_genome.columns.to_list())
        self.assertTrue('anchor_index' in mdl_parms.df_J_ref_genome.columns.to_list())

    def test_IgorModel_Parms_with_anchors(self):
        species = "human"
        chain = "tcr_beta"
        fln_model_parms, fln_model_marginals = get_default_models_paths_species_chain(species, chain)
        fln_dict = get_default_fln_dict_ref_genomes_species_chain(species, chain)
        print("fln_dict: ", fln_dict)

        mdl_parms = IgorModel_Parms(model_parms_file=fln_model_parms,
                                    fln_V_gene_CDR3_anchors=fln_dict["fln_V_gene_CDR3_anchors"],
                                    fln_J_gene_CDR3_anchors=fln_dict["fln_J_gene_CDR3_anchors"])
        # no function in anchors file
        self.assertTrue('anchor_index' in mdl_parms.df_V_ref_genome.columns.to_list())
        self.assertTrue('anchor_index' in mdl_parms.df_J_ref_genome.columns.to_list())


    def test_IgorModel_Parms_from_IgorRefGenome(self):
        ref_genome = IgorRefGenome.load_default("human", "tcr_alpha")
        self.assertIsInstance(ref_genome, IgorRefGenome)
        ref_genome_dict = ref_genome.to_dict()
        print(ref_genome.df_genomicVs)
        print(ref_genome.df_V_anchors)
        print(IgorRefGenome.V)
        mdl_parms = IgorModel_Parms.make_default_VDJ_from_IgorRefGenome(ref_genome)
        # mdl_parms = IgorModel_Parms.make_default_VDJ(ref_genome.df_genomicVs, ref_genome.df_genomicDs, ref_genome.df_genomicJs)
        # mdl_parms.attach_anchors_from_files()
        # IgorModel_Parms.make_default_VDJ_from_IgorRefGenome()
        # print(mdl_parms.df_V_anchors)
        # print(mdl_parms.df_J_anchors)

    def test_IgorModel_Parms_get_IgorRefGenome_VDJ(self):
        species = "human"
        chain = "tcr_beta"
        fln_model_parms, fln_model_marginals = get_default_models_paths_species_chain(species, chain)
        fln_dict = get_default_fln_dict_ref_genomes_species_chain(species, chain)

        mdl_parms = IgorModel_Parms(model_parms_file=fln_model_parms,
                                    fln_V_gene_CDR3_anchors=fln_dict["fln_V_gene_CDR3_anchors"],
                                    fln_J_gene_CDR3_anchors=fln_dict["fln_J_gene_CDR3_anchors"])

        self.assertFalse(mdl_parms.df_V_anchors is None)
        self.assertFalse(mdl_parms.df_J_anchors is None)
        self.assertFalse(mdl_parms.event_GeneChoice_D is None)
        self.assertFalse(mdl_parms.df_D_ref_genome is None)
        ref_genome = mdl_parms.get_IgorRefGenome()
        mdl_parms.gen_EventDict_DataFrame()

        tmp_ref_genome_dir = tempfile.TemporaryDirectory(dir='.', prefix="ref_genome")
        ref_genome.write_ref_genome_dir(tmp_ref_genome_dir.name)
        print(ref_genome.to_dict())

        fln_tmp_dict = dict()
        fln_tmp_dict['fln_genomicVs'] = tmp_ref_genome_dir.name + "/" + "genomicVs.fasta"
        fln_tmp_dict['fln_genomicDs'] = tmp_ref_genome_dir.name + "/" + "genomicDs.fasta"
        fln_tmp_dict['fln_genomicJs'] = tmp_ref_genome_dir.name + "/" + "genomicJs.fasta"
        fln_tmp_dict['fln_V_gene_CDR3_anchors'] = tmp_ref_genome_dir.name + "/" + "V_gene_CDR3_anchors.csv"
        fln_tmp_dict['fln_J_gene_CDR3_anchors'] = tmp_ref_genome_dir.name + "/" + "J_gene_CDR3_anchors.csv"
        # time.sleep(20)

        self.assertTrue(os.path.isfile(fln_tmp_dict["fln_genomicVs"]))
        self.assertTrue(os.path.isfile(fln_tmp_dict["fln_genomicDs"]))
        self.assertTrue(os.path.isfile(fln_tmp_dict["fln_genomicJs"]))
        self.assertTrue(os.path.isfile(fln_tmp_dict["fln_V_gene_CDR3_anchors"]))
        self.assertTrue(os.path.isfile(fln_tmp_dict["fln_J_gene_CDR3_anchors"]))

        tmp_ref_genome_dir.cleanup()

    def test_IgorModel_Parms_get_IgorRefGenome_VJ(self):
        """
        Return an IgorRefGenome object generated from GeneChoice events
        and use it to write a ref_genome directory that will be use to run IGoR.
        """
        species = "human"
        chain = "tcr_alpha"
        fln_model_parms, fln_model_marginals = get_default_models_paths_species_chain(species, chain)
        fln_dict = get_default_fln_dict_ref_genomes_species_chain(species, chain)

        # 1. Make an IgorModel_Parms from scratch
        mdl_parms = IgorModel_Parms(model_parms_file=fln_model_parms,
                                    fln_V_gene_CDR3_anchors=fln_dict["fln_V_gene_CDR3_anchors"],
                                    fln_J_gene_CDR3_anchors=fln_dict["fln_J_gene_CDR3_anchors"])

        self.assertTrue(mdl_parms.event_GeneChoice_D is None)
        self.assertTrue(mdl_parms.df_D_ref_genome is None)

        # 2. Get IgorRefGenome from events
        ref_genome = mdl_parms.get_IgorRefGenome()
        # 3. Write a ref_genome_dir
        tmp_ref_genome_dir = tempfile.TemporaryDirectory(dir='.', prefix="ref_genome")
        ref_genome.write_ref_genome_dir(tmp_ref_genome_dir.name)

        fln_tmp_dict = dict()
        fln_tmp_dict['fln_genomicVs'] = tmp_ref_genome_dir.name + "/" + "genomicVs.fasta"
        fln_tmp_dict['fln_genomicJs'] = tmp_ref_genome_dir.name + "/" + "genomicJs.fasta"
        fln_tmp_dict['fln_V_gene_CDR3_anchors'] = tmp_ref_genome_dir.name + "/" + "V_gene_CDR3_anchors.csv"
        fln_tmp_dict['fln_J_gene_CDR3_anchors'] = tmp_ref_genome_dir.name + "/" + "J_gene_CDR3_anchors.csv"

        self.assertTrue(os.path.isfile(fln_tmp_dict["fln_genomicVs"]))
        self.assertTrue(os.path.isfile(fln_tmp_dict["fln_genomicJs"]))
        self.assertTrue(os.path.isfile(fln_tmp_dict["fln_V_gene_CDR3_anchors"]))
        self.assertTrue(os.path.isfile(fln_tmp_dict["fln_J_gene_CDR3_anchors"]))

        self.assertFalse(os.path.isfile(tmp_ref_genome_dir.name + "/" + "genomicDs.fasta"))

        tmp_ref_genome_dir.cleanup()

    def tearDown(self) -> None:
        self.tmp_dir.cleanup()

if __name__ == '__main__':
    unittest.main()
