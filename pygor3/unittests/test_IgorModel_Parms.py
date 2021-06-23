import unittest
import tempfile
import matplotlib.pyplot as plt
from pygor3.utils import *
from pygor3 import IgorModel_Parms
from pygor3 import IgorRec_Event
from pygor3 import IgorRefGenome
from pygor3 import *
import os
import io

import time


str_mock_VDJ_fln_genomicVs = \
""">TRBV1*01
GATACTGGAATTACCCAGACACCAAAATACCTGGTCACAGCAATGGGGAGTAAAAGGACA
ATGAAACGTGAGCATCTGGGACATGATTCTATGTATTGGTACAGACAGAAAGCTAAGAAA
TCCCTGGAGTTCATGTTTTACTACAACTGTAAGGAATTCATTGAAAACAAGACTGTGCCA
AATCACTTCACACCTGAATGCCCTGACAGCTCTCGCTTATACCTTCATGTGGTCGCACTG
CAGCAAGAAGACTCAGCTGCGTATCTCTGCACCAGCAGCCAAGA
>TRBV2*01
GAACCTGAAGTCACCCAGACTCCCAGCCATCAGGTCACACAGATGGGACAGGAAGTGATC
TTGCGCTGTGTCCCCATCTCTAATCACTTATACTTCTATTGGTACAGACAAATCTTGGGG
CAGAAAGTCGAGTTTCTGGTTTCCTTTTATAATAATGAAATCTCAGAGAAGTCTGAAATA
TTCGATGATCAATTCTCAGTTGAAAGGCCTGATGGATCAAATTTCACTCTGAAGATCCGG
TCCACAAAGCTGGAGGACTCAGCCATGTACTTCTGTGCCAGCAGTGAAGC
>TRBV2*02
GAACCTGAAGTCACCCAGACTCCCAGCCATCAGGTCACACAGATGGGACAGGAAGTGATC
TTGCACTGTGTCCCCATCTCTAATCACTTATACTTCTATTGGTACAGACAAATCTTGGGG
CAGAAAGTCGAGTTTCTGGTTTCCTTTTATAATAATGAAATCTCAGAGAAGTCTGAAATA
TTCGATGATCAATTCTCAGTTGAAAGGCCTGATGGATCAAATTTCACTCTGAAGATCCGG
TCCACAAAGCTGGAGGACTCAGCCATGTACTTCTGTGCCAGCAGT
"""
str_mock_VDJ_fln_genomicDs = \
""">TRBD1*01
GGGACAGGGGGC
>TRBD2*01
GGGACTAGCGGGGGGG
>TRBD2*02
GGGACTAGCGGGAGGG
"""
str_mock_VDJ_fln_genomicJs = \
""">TRBJ1-1*01
TGAACACTGAAGCTTTCTTTGGACAAGGCACCAGACTCACAGTTGTAG
>TRBJ1-2*01
CTAACTATGGCTACACCTTCGGTTCGGGGACCAGGTTAACCGTTGTAG
>TRBJ1-3*01
CTCTGGAAACACCATATATTTTGGAGAGGGAAGTTGGCTCACTGTTGTAG
>TRBJ1-4*01
CAACTAATGAAAAACTGTTTTTTGGCAGTGGAACCCAGCTCTCTGTCTTGG
>TRBJ1-5*01
TAGCAATCAGCCCCAGCATTTTGGTGATGGGACTCGACTCTCCATCCTAG
"""

str_mock_VDJ_fln_V_gene_CDR3_anchors = \
"""gene;anchor_index;gfunction
TRBV1*01;267;P
TRBV2*01;273;F
TRBV2*02;273;(F)
TRBV2*03;273;(F)
TRBV3-1*01;270;F
"""

str_mock_VDJ_fln_J_gene_CDR3_anchors = \
"""gene;anchor_index;function
TRBJ1-1*01;17;F
TRBJ1-2*01;17;F
TRBJ1-3*01;19;F
TRBJ1-4*01;20;F
"""

class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp_dir = tempfile.TemporaryDirectory(dir=".", prefix="models")
        self.fln_dict = get_default_fln_names_for_model_dir(self.tmp_dir.name)

    def test_IgorModel_Parms_from_dataframe(self):
        mdl_hb = get_default_IgorModel("human", "tcr_beta")
        import copy
        genomic_dict = copy.deepcopy(mdl_hb.genomic_dataframe_dict)

        genomic_dict['V']['name'] = v_genLabel(genomic_dict['V']['name'])
        genomic_dict['J']['name'] = v_genLabel(genomic_dict['J']['name'])

        new_V_gene_dict = {
            'name': 'my_pseudo_TRBV',
            'value': 'AAACCCTTTGGGACCCAGAGCCCAAGACACAAGATCACAGAGACAGGAAGGCAGGTGACCTTGGCGTGTCACCAGACTTGGAACCACAACAATATGTTCTGGTATCGACAAGACCTGGGACATGGGCTGAGGCTGATCCATTACTCATATGGTGTTCACGACACTAACAAAGGAGAAGTCTCAGATGGCTACAGTGTCTCTAGATCAAACACAGAGGACCTCCCCCTCACTCTGTAGTCTGCTGCCTCCTCCCAGACATCTGTATATTTCTGCGCCAGCAGTGAGTC',
            'anchor_index': 270
        }
        df_V = genomic_dict['V'].loc[10:15]
        df_V = df_V.append(new_V_gene_dict, ignore_index=True)
        df_V.index.name = 'id'
        df_V


        mdl_parms_0 = IgorModel_Parms.make_default_VDJ(df_V, genomic_dict['D'], genomic_dict['J'])
        mdl_marginals_0 = IgorModel_Marginals.make_uniform_from_parms(mdl_parms_0)
        print("parms.Edges: ")
        print(mdl_hb.parms.Edges)
        print(mdl_parms_0.Edges)

        print("parms.Edges_dict:")
        print(mdl_hb.parms.Edges_dict)
        print(mdl_parms_0.Edges_dict)
        print("marginals.network_dict")
        print(mdl_hb.marginals.network_dict)
        print(mdl_marginals_0.network_dict)

        self.assertIsInstance(mdl_parms_0, IgorModel_Parms)

    def test_get_df_ref_genome_from_files(self):
        ofile_mock_VDJ_fln_genomicVs = io.StringIO(str_mock_VDJ_fln_genomicVs)
        ofile_mock_VDJ_fln_V_gene_CDR3_anchors = io.StringIO(str_mock_VDJ_fln_V_gene_CDR3_anchors)

        ofile_mock_VDJ_fln_genomicDs = io.StringIO(str_mock_VDJ_fln_genomicDs)

        ofile_mock_VDJ_fln_genomicJs = io.StringIO(str_mock_VDJ_fln_genomicJs)
        ofile_mock_VDJ_fln_J_gene_CDR3_anchors = io.StringIO(str_mock_VDJ_fln_J_gene_CDR3_anchors)

        df_V_ref_genome = get_dataframe_from_fasta_and_csv_anchors(ofile_mock_VDJ_fln_genomicVs,
                                                                   ofile_mock_VDJ_fln_V_gene_CDR3_anchors)
        df_D_ref_genome = get_dataframe_from_fasta_and_csv_anchors(ofile_mock_VDJ_fln_genomicDs)

        df_J_ref_genome = get_dataframe_from_fasta_and_csv_anchors(ofile_mock_VDJ_fln_genomicJs,
                                                                   ofile_mock_VDJ_fln_J_gene_CDR3_anchors)
        self.assertIsInstance(df_V_ref_genome, pd.DataFrame)
        self.assertIsInstance(df_D_ref_genome, pd.DataFrame)
        self.assertIsInstance(df_J_ref_genome, pd.DataFrame)

        mdl_parms = IgorModel_Parms.make_default_VDJ(df_V_ref_genome, df_D_ref_genome, df_J_ref_genome)
        self.assertIsInstance(mdl_parms, IgorModel_Parms)

    def test_get_df_ref_genome_from_files_VJ(self):
        ofile_mock_VJ_fln_genomicVs = io.StringIO(str_mock_VDJ_fln_genomicVs)
        ofile_mock_VJ_fln_V_gene_CDR3_anchors = io.StringIO(str_mock_VDJ_fln_V_gene_CDR3_anchors)

        ofile_mock_VJ_fln_genomicDs = io.StringIO(str_mock_VDJ_fln_genomicDs)

        ofile_mock_VJ_fln_genomicJs = io.StringIO(str_mock_VDJ_fln_genomicJs)
        ofile_mock_VJ_fln_J_gene_CDR3_anchors = io.StringIO(str_mock_VDJ_fln_J_gene_CDR3_anchors)

        df_V_ref_genome = get_dataframe_from_fasta_and_csv_anchors(ofile_mock_VJ_fln_genomicVs,
                                                                   ofile_mock_VJ_fln_V_gene_CDR3_anchors)

        df_J_ref_genome = get_dataframe_from_fasta_and_csv_anchors(ofile_mock_VJ_fln_genomicJs,
                                                                   ofile_mock_VJ_fln_J_gene_CDR3_anchors)
        self.assertIsInstance(df_V_ref_genome, pd.DataFrame)
        self.assertIsInstance(df_J_ref_genome, pd.DataFrame)

        mdl_parms = IgorModel_Parms.make_default_VJ(df_V_ref_genome, df_J_ref_genome)
        self.assertIsInstance(mdl_parms, IgorModel_Parms)

        # mdl_parms.plot_Graph() # DON'T KNOW WHY the plot is a mess.
        # plt.show()

        # get_dataframe_from_fasta_and_csv_anchors()
        # mdl_parms.make_default_VDJ()
        # print(Igor_VDJ_default_nickname_list)
        # print(IgorRec_Event_default_dict)
        # GeneChoice;V_gene;Undefined_side;7;v_choice

        # event_type ='GeneChoice'
        # seq_type = '' seq_side, priority,
        # nickname
        # event = IgorRec_Event(event_type, seq_type, seq_side, priority,
        #          nickname)
        # mdl_parms.add_Event()

    def test_set_event_from_dataframe(self):
        mdl_copy = IgorModel.load_default("human", "tcr_beta")
        mdl_copy.genomic_dataframe_dict['V']
        mdl_copy.parms['v_choice']['name'] = v_genLabel(mdl_copy.parms['v_choice']['name'])
        mdl_copy.parms['j_choice']['name'] = v_genLabel(mdl_copy.parms['j_choice']['name'])
        mdl_copy.parms.Event_dict['v_choice']
        help(mdl_copy.parms.set_event_realizations_from_DataFrame)
        mdl_copy.parms.set_event_realizations_from_DataFrame('v_choice', mdl_copy.parms['v_choice'])
        mdl_copy.parms.set_event_realizations_from_DataFrame('j_choice', mdl_copy.parms['j_choice'])
        mdl_copy.generate_xdata()

        # TODO: HOW TO EDIT MODEL MARGINALS.
        mdl_hb.genomic_dataframe_dict['V'].loc[0].value
        mdl_hb.genomic_dataframe_dict['V'].loc[0].anchor_index
        mdl_hb.set_genomic_dataframe_dict(dictio)
        mdl_hb.generate_xdata()




        # mdl_copy.marginals.marginals_dict['v_choice']

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
