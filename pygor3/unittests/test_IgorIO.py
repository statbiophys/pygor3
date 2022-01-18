#!/usr/bin/env python3
import unittest

from pygor3 import *
from pygor3.imgt import *
import os
import pandas as pd
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

    def test_AIRR_from_scenario(self):
        species = "human"
        chain = "tcr_beta"

        str_sequence = "CAAGACCCAGGACTGGGCCTACGGTTGATCTATTACTCCTTTGATGTCAAAGATATAAACAAAGGAGAGATCTCTGATGGATACAGTGTCTCTCGACAGGCACAGGCTAAATTCTCCCTGTCCCTAGAGTCTGCCATCCCCAACCAGACAGCTCTTTACTTCTGTGCCACTCCCCCGGTGGCTGGCTACACCTTCGGTTCGGGGACCAGGTTAACCGTTGTAG"

        airr_fields = AIRR_VDJ_rearrangement.list_of_fields()
        mdl = get_default_IgorModel(species, chain)
        df_scenarios = evaluate(str_sequence, mdl=mdl, N_scenarios=20)
        self.assertIsInstance(df_scenarios, pd.DataFrame)
        print(df_scenarios)
        ps_scenario = df_scenarios.iloc[0]
        aaa = mdl.get_AIRR_from_ps_scenario(ps_scenario, v_offset=0)
        print(aaa)
        # df_ps_scenario = mdl.get_df_scenario_aln_from_scenario(ps_scenario)
        # print(df_ps_scenario)
        # ps_scenario
        # v_offset
        # airr_rearrangement_dict = mdl.get_AIRR_VDJ_rearragement_dict_from_scenario(ps_scenario, str_sequence, v_offset=0)
        # print(airr_rearrangement_dict)

    def test_AIRR_format(self):
        import airr

        # Cigar alignment from scenario and V_offset
        directory_name = os.path.dirname(flnAIRR_arrangement)
        pathlib.Path(directory_name).mkdir(parents=True, exist_ok=True)

        # IF VDJ THEN:
        b_D_gene = (len([event.event_type for event in mdl.parms.Event_list if event.nickname == 'd_gene']) > 0)

        if b_D_gene:
            airr_fields = AIRR_VDJ_rearrangement.list_of_fields()
            airr_rearrangement_writer = airr.create_rearrangement(flnAIRR_arrangement, fields=airr_fields)
            for seq_index, sequence in self.fetch_IgorIndexedSeq_records():
                scenarios_list = self.get_IgorBestScenarios_By_seq_index_IgorModel(seq_index, mdl)
                for scenario in scenarios_list:
                    v_best_aln = self.get_best_IgorAlignment_data_By_seq_index("V", seq_index)

                    airr_rearrangement_dict = mdl.get_AIRR_VDJ_rearragement_dict_from_scenario(scenario, sequence,
                                                                                               v_offset=v_best_aln.offset)
                    airr_rearrangement_dict['scenario_rank'] = scenario.scenario_rank
                    airr_rearrangement_dict['scenario_proba_cond_seq'] = scenario.scenario_proba_cond_seq
                    airr_rearrangement_dict['pgen'] = self.fetch_IgorPgen_By_seq_index(seq_index)[1]
                    cdr3_record = self.fetch_IgorIndexedCDR3_By_seq_index(seq_index)
                    airr_rearrangement_dict['junction'] = cdr3_record[3]
                    airr_rearrangement_dict['junction_aa'] = cdr3_record[4]
                    airr_rearrangement_writer.write(airr_rearrangement_dict)

            airr_rearrangement_writer.close()

        else:
            # FIXME: IF VJ THEN:
            airr_fields = AIRR_VDJ_rearrangement.list_of_fields()
            airr_rearrangement_writer = airr.create_rearrangement(flnAIRR_arrangement, fields=airr_fields)




if __name__ == '__main__':
    unittest.main()

