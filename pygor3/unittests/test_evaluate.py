import unittest

import pandas as pd

from pygor3 import IgorTask
from pygor3 import IgorModel
from pygor3 import IgorRefGenome
from pygor3 import generate
from pygor3 import evaluate
from pygor3 import evaluate_pgen
from pygor3 import infer
from pygor3 import get_default_IgorModel
from pygor3 import get_IgorRefGenome_VDJ_from_IMGT
from pygor3 import rcParams
from pygor3.utils import dna_translate

from pygor3 import from_df_scenario_aln_to_da_scenario_aln, plot_scenario_from_da_scenario_aln
import matplotlib.pyplot as plt
import matplotlib as mpl



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
>TRBV2*03
GAACCTGAAGTCACCCAGACTCCCAGCCATCAGGTCACACAGATGGGACAGGAAGTGATC
TTGCGCTGTGTCCCCATCTCTAATCACTTATACTTCTATTGGTACAGACAAATCTTGGGG
CAGAAAGTCGAGTTTCTGGTTTCCTTTTATAATAATGAAATCTCAGAGAAGTCTGAAATA
TTCGATGATCAATTCTCAGTTGAGAGGCCTGATGGATCAAATTTCACTCTGAAGATCCGG
TCCACAAAGCTGGAGGACTCAGCCATGTACTTCTGTGCCAGCAGTGAA
>TRBV3-1*01
GACACAGCTGTTTCCCAGACTCCAAAATACCTGGTCACACAGATGGGAAACGACAAGTCC
ATTAAATGTGAACAAAATCTGGGCCATGATACTATGTATTGGTATAAACAGGACTCTAAG
AAATTTCTGAAGATAATGTTTAGCTACAATAATAAGGAGCTCATTATAAATGAAACAGTT
CCAAATCGCTTCTCACCTAAATCTCCAGACAAAGCTCACTTAAATCTTCACATCAATTCC
CTGGAGCTTGGTGACTCTGCTGTGTATTTCTGTGCCAGCAGCCAAGA
>TRBV3-1*02
GACACAGCTGTTTCCCAGACTCCAAAATACCTGGTCACACAGATGGGAAACGACAAGTCC
ATTAAATGTGAACAAAATCTGGGCCATGATACTATGTATTGGTATAAACAGGACTCTAAG
AAATTTCTGAAGATAATGTTTAGCTACAATAACAAGGAGATCATTATAAATGAAACAGTT
CCAAATCGATTCTCACCTAAATCTCCAGACAAAGCTAAATTAAATCTTCACATCAATTCC
CTGGAGCTTGGTGACTCTGCTGTGTATTTCTGTGCCAGC
>TRBV3-2*01
GACACAGCCGTTTCCCAGACTCCAAAATACCTGGTCACACAGATGGGAAAAAAGGAGTCT
CTTAAATGAGAACAAAATCTGGGCCATAATGCTATGTATTGGTATAAACAGGACTCTAAG
AAATTTCTGAAGACAATGTTTATCTACAGTAACAAGGAGCCAATTTTAAATGAAACAGTT
CCAAATCGCTTCTCACCTGACTCTCCAGACAAAGCTCATTTAAATCTTCACATCAATTCC
CTGGAGCTTGGTGACTCTGCTGTGTATTTCTGTGCCAGCAGCCAAGA
>TRBV3-2*02
GACACAGCCGTTTCCCAGACTCCAAAATACCTGGTCACACAGATGGGAAAAAAGGAGTCT
CTTAAATGAGAACAAAATCTGGGCCATAATGCTATGTATTGGTATAAACAGGACTCTAAG
AAATTTCTGAAGACAATGTTTATCTACAGTAACAAGGAGCCAATTTTAAATGAAACAGTT
CCAAATCGCTTCTCACCTGACTCTCCAGACAAAGTTCATTTAAATCTTCACATCAATTCC
CTGGAGCTTGGTGACTCTGCTGTGTATTTCTGTGCCAGCAGCCAAGA
>TRBV3-2*03
GACACAGCCGTTTCCCAGACTCCAAAATACCTGGTCACACAGACGGGAAAAAAGGAGTCT
CTTAAATGAGAACAAAATCTGGGCCATAATGCTATGTATTGGTATAAACAGGACTCTAAG
AAATTTCTGAAGACAATGTTTATCTACAGTAACAAGGAGCCAATTTTAAATGAAACAGTT
CCAAATCGCTTCTCACCTGACTCTCCAGACAAAGTTCATTTAAATCTTCACATCAATTCC
CTGGAGCTTGGTGACTCTGCTGTGTATTTCTGTGCCAGCAGCCAAG
>TRBV4-1*02
CACCTGGTCATGGGAATGACAAATAAGAAGTCTTTGAAATGTGAACAACATATGGGGCAC
AGGGCAATGTATTGGTACAAGCAGAAAGCTAAGAAGCCACCGGAGCTCATGTTTGTCTAC
AGCTATGAGAAACTCTCTATAAATGAAAGTGTGCCAAGTCGCTTCTCACCTGAATGCCCC
AACAGCTCTCTCTTAAACCTTCACCTACACGCCCTGCAGCCAGAAGACTCAGCCCTGTAT
CTCTGCGCCAGCAGCCAAG
>TRBV4-2*01
GAAACGGGAGTTACGCAGACACCAAGACACCTGGTCATGGGAATGACAAATAAGAAGTCT
TTGAAATGTGAACAACATCTGGGGCATAACGCTATGTATTGGTACAAGCAAAGTGCTAAG
AAGCCACTGGAGCTCATGTTTGTCTACAACTTTAAAGAACAGACTGAAAACAACAGTGTG
CCAAGTCGCTTCTCACCTGAATGCCCCAACAGCTCTCACTTATTCCTTCACCTACACACC
CTGCAGCCAGAAGACTCGGCCCTGTATCTCTGTGCCAGCAGCCAAGA
>TRBV4-2*02
GAAACGGGAGTTACGCAGACACCAAGACACCTGGTCATGGGAATGACAAATAAGAAGTCT
TTGAAATGTGAACAACATCTGGGGCATAACGCTATGTATTGGTACAAGCAAAGTGCTAAG
AAGCCACTGGAGCTCATGTTTGTCTACAACTTTAAAGAACAGACTGAAAACAACAGTGTG
CCAAGTCGCTTCTCACCTGAATGCCCCAACAGCTCTCACTTATGCCTTCACCTACACACC
CTGCAGCCAGAAGACTCGGCCCTGTATCTCTGTGCCAGCACC
>TRBV4-3*01
GAAACGGGAGTTACGCAGACACCAAGACACCTGGTCATGGGAATGACAAATAAGAAGTCT
TTGAAATGTGAACAACATCTGGGTCATAACGCTATGTATTGGTACAAGCAAAGTGCTAAG
AAGCCACTGGAGCTCATGTTTGTCTACAGTCTTGAAGAACGGGTTGAAAACAACAGTGTG
CCAAGTCGCTTCTCACCTGAATGCCCCAACAGCTCTCACTTATTCCTTCACCTACACACC
CTGCAGCCAGAAGACTCGGCCCTGTATCTCTGCGCCAGCAGCCAAGA
>TRBV4-3*02
GAAACGGGAGTTACGCAGACACCAAGACACCTGGTCATGGGAATGACAAATAAGAAGTCT
TTGAAATGTGAACAACATCTGGGTCATAACGCTATGTATTGGTACAAGCAAAGTGCTAAG
AAGCCACTGGAGCTCATGTTTGTCTACAGTCTTGAAGAACGGGTTGAAAACAACAGTGTG
CCAAGTCGCTTCTCACCTGAATGCCCCAACAGCTCTCACTTATCCCTTCACCTACACACC
CTGCAGCCAGAAGACTCGGCCCTGTATCTCTGCGCCAGCAGC
>TRBV4-3*03
GAAACGGGAGTTACGCAGACACCAAGACACCTGGTCATGGGAATGACAAATAAGAAGTCT
TTGAAATGTGAACAACATCTGGGTCATAACGCTATGTATTGGTACAAGCAAAGTGCTAAG
AAGCCACTGGAGCTCATGTTTGTCTACAGTCTTGAAGAACGTGTTGAAAACAACAGTGTG
CCAAGTCGCTTCTCACCTGAATGCCCCAACAGCTCTCACTTATTCCTTCACCTACACACC
CTGCAGCCAGAAGACTCGGCCCTGTATCTCTGCGCCAGCAGC
>TRBV4-3*04
AAGAAGTCTTTGAAATGTGAACAACATCTGGGGCATAACGCTATGTATTGGTACAAGCAA
AGTGCTAAGAAGCCACTGGAGCTCATGTTTGTCTACAGTCTTGAAGAACGGGTTGAAAAC
AACAGTGTGCCAAGTCGCTTCTCACCTGAATGCCCCAACAGCTCTCACTTATTCCTTCAC
CTACACACCCTGCAGCCAGAAGACTCGGCCCTGTATCTCTGCGCCAGCAGC
>TRBV5-1*01
AAGGCTGGAGTCACTCAAACTCCAAGATATCTGATCAAAACGAGAGGACAGCAAGTGACA
CTGAGCTGCTCCCCTATCTCTGGGCATAGGAGTGTATCCTGGTACCAACAGACCCCAGGA
CAGGGCCTTCAGTTCCTCTTTGAATACTTCAGTGAGACACAGAGAAACAAAGGAAACTTC
CCTGGTCGATTCTCAGGGCGCCAGTTCTCTAACTCTCGCTCTGAGATGAATGTGAGCACC
TTGGAGCTGGGGGACTCGGCCCTTTATCTTTGCGCCAGCAGCTTGG
>TRBV5-1*02
AGGGCTGGGGTCACTCAAACTCCAAGACATCTGATCAAAACGAGAGGACAGCAAGTGACA
CTGGGCTGCTCCCCTATCTCTGGGCATAGGAGTGTATCCTGGTACCAACAGACCCTAGGA
CAGGGCCTTCAGTTCCTCTTTGAATACTTCAGTGAGACACAGAGAAACAAAGGAAACTTC
CTTGGTCGATTCTCAGGGCGCCAGTTCTCTAACTCTCGCTCTGAGATGAATGTGAGCACC
TTGGAGCTGGGGGACTCGGCCCTTTATCTTTGCGCCAGC
>TRBV5-2*01
GAGGCTGGAATCACCCAAGCTCCAAGACACCTGATCAAAACAAGAGACCAGCAAGTGACA
CTGAGATGCTCCCCTGCCTCTGGGCATAACTGTGTGTCCTGGTACCTACGAACTCCAAGT
CAGCCCCTCTAGTTATTGTTACAATATTGTAATAGGTTACAAAGAGCAAAAGGAAACTTG
CCTAATTGATTCTCAGCTCACCACGTCCATAACTATTACTGAGTCAAACACGGAGCTAGG
GGACTCAGCCCTGTATCTCTGTGCCAGCAACTTGATG
>TRBV5-3*01
GAGGCTGGAGTCACCCAAAGTCCCACACACCTGATCAAAACGAGAGGACAGCAAGTGACT
CTGAGATGCTCTCCTATCTCTGGGCACAGCAGTGTGTCCTGGTACCAACAGGCCCCGGGT
CAGGGGCCCCAGTTTATCTTTGAATATGCTAATGAGTTAAGGAGATCAGAAGGAAACTTC
CCTAATCGATTCTCAGGGCGCCAGTTCCATGACTGTTGCTCTGAGATGAATGTGAGTGCC
TTGGAGCTGGGGGACTCGGCCCTGTATCTCTGTGCCAGAAGCTTGG
>TRBV5-3*02
GAGGCTGGAGTCACCCAAAGTCCCACACACCTGATCAAAACGAGAGGACAGCAAGTGACT
CTGAGATGCTCTCCTATCTCTGGGCACAGCAGTGTGTCCTGGTACCAACAGGCCCCGGGT
CAGGGGCCCCAGTTTATCTTTGAATATGCTAATGAGTTAAGGAGATCAGAAGGAAACTTC
CCTAATCGATTCTCAGGGCGCCAGTTCCATGACTATTGCTCTGAGATGAATGTGAGTGCC
TTGGAGCTGGGGGACTCGGCCCTGTATCTCTGTGCCAGAAGCTTGG
>TRBV5-4*01
GAGACTGGAGTCACCCAAAGTCCCACACACCTGATCAAAACGAGAGGACAGCAAGTGACT
CTGAGATGCTCTTCTCAGTCTGGGCACAACACTGTGTCCTGGTACCAACAGGCCCTGGGT
CAGGGGCCCCAGTTTATCTTTCAGTATTATAGGGAGGAAGAGAATGGCAGAGGAAACTTC
CCTCCTAGATTCTCAGGTCTCCAGTTCCCTAATTATAGCTCTGAGCTGAATGTGAACGCC
TTGGAGCTGGACGACTCGGCCCTGTATCTCTGTGCCAGCAGCTTGG
>TRBV5-4*02
GAGACTGGAGTCACCCAAAGTCCCACACACCTGATCAAAACGAGAGGACAGCAAGTGACT
CTGAGATGCTCTTCTCAGTCTGGGCACAACACTGTGTCCTGGTACCAACAGGCCCTGGGT
CAGGGGCCCCAGTTTATCTTTCAGTATTATAGGGAGGAAGAGAATGGCAGAGGAAACTTC
CCTCCTAGATTCTCAGGTCTCCAGTTCCCTAATTATAACTCTGAGCTGAATGTGAACGCC
TTGGAGCTGGACGACTCGGCCCTGTATCTCTGTGCCAGCAGC
>TRBV5-4*03
CAGCAAGTGACACTGAGATGCTCTTCTCAGTCTGGGCACAACACTGTGTCCTGGTACCAA
CAGGCCCTGGGTCAGGGGCCCCAGTTTATCTTTCAGTATTATAGGGAGGAAGAGAATGGC
AGAGGAAACTTCCCTCCTAGATTCTCAGGTCTCCAGTTCCCTAATTATAGCTCTGAGCTG
AATGTGAACGCCTTGGAGCTGGACGACTCGGCCCTGTATCTCTGTGCCAGCAGC
>TRBV5-4*04
ACTGTGTCCTGGTACCAACAGGCCCTGGGTCAGGGGCCCCAGTTTATCTTTCAGTATTAT
AGGGAGGAAGAGAATGGCAGAGGAAACTCCCCTCCTAGATTCTCAGGTCTCCAGTTCCCT
AATTATAGCTCTGAGCTGAATGTGAACGCCTTGGAGCTGGACGACTCGGCCCTGTATCTC
TGTGCCAGCAGC
>TRBV5-5*01
GACGCTGGAGTCACCCAAAGTCCCACACACCTGATCAAAACGAGAGGACAGCAAGTGACT
CTGAGATGCTCTCCTATCTCTGGGCACAAGAGTGTGTCCTGGTACCAACAGGTCCTGGGT
CAGGGGCCCCAGTTTATCTTTCAGTATTATGAGAAAGAAGAGAGAGGAAGAGGAAACTTC
CCTGATCGATTCTCAGCTCGCCAGTTCCCTAACTATAGCTCTGAGCTGAATGTGAACGCC
TTGTTGCTGGGGGACTCGGCCCTGTATCTCTGTGCCAGCAGCTTGG
>TRBV5-5*02
GACGCTGGAGTCACCCAAAGTCCCACACACCTGATCAAAACGAGAGGACAGCACGTGACT
CTGAGATGCTCTCCTATCTCTGGGCACAAGAGTGTGTCCTGGTACCAACAGGTCCTGGGT
CAGGGGCCCCAGTTTATCTTTCAGTATTATGAGAAAGAAGAGAGAGGAAGAGGAAACTTC
CCTGATCGATTCTCAGCTCGCCAGTTCCCTAACTATAGCTCTGAGCTGAATGTGAACGCC
TTGTTGCTGGGGGACTCGGCCCTGTATCTCTGTGCCAGCAGC
>TRBV5-5*03
GACGCTGGAGTCACCCAAAGTCCCACACACCTGATCAAAACGAGAGGACAGCAAGTGACT
CTGAGATGCTCTCCTATCTCTGAGCACAAGAGTGTGTCCTGGTACCAACAGGTCCTGGGT
CAGGGGCCCCAGTTTATCTTTCAGTATTATGAGAAAGAAGAGAGAGGAAGAGGAAACTTC
CCTGATCGATTCTCAGCTCGCCAGTTCCCTAACTATAGCTCTGAGCTGAATGTGAACGCC
TTGTTGCTGGGGGACTCGGCCCTGTATCTCTGTGCCAGCAGC
>TRBV5-6*01
GACGCTGGAGTCACCCAAAGTCCCACACACCTGATCAAAACGAGAGGACAGCAAGTGACT
CTGAGATGCTCTCCTAAGTCTGGGCATGACACTGTGTCCTGGTACCAACAGGCCCTGGGT
CAGGGGCCCCAGTTTATCTTTCAGTATTATGAGGAGGAAGAGAGACAGAGAGGCAACTTC
CCTGATCGATTCTCAGGTCACCAGTTCCCTAACTATAGCTCTGAGCTGAATGTGAACGCC
TTGTTGCTGGGGGACTCGGCCCTCTATCTCTGTGCCAGCAGCTTGG
>TRBV5-7*01
GACGCTGGAGTCACCCAAAGTCCCACACACCTGATCAAAACGAGAGGACAGCACGTGACT
CTGAGATGCTCTCCTATCTCTGGGCACACCAGTGTGTCCTCGTACCAACAGGCCCTGGGT
CAGGGGCCCCAGTTTATCTTTCAGTATTATGAGAAAGAAGAGAGAGGAAGAGGAAACTTC
CCTGATCAATTCTCAGGTCACCAGTTCCCTAACTATAGCTCTGAGCTGAATGTGAACGCC
TTGTTGCTAGGGGACTCGGCCCTCTATCTCTGTGCCAGCAGCTTGG
>TRBV5-8*01
GAGGCTGGAGTCACACAAAGTCCCACACACCTGATCAAAACGAGAGGACAGCAAGCGACT
CTGAGATGCTCTCCTATCTCTGGGCACACCAGTGTGTACTGGTACCAACAGGCCCTGGGT
CTGGGCCTCCAGTTCCTCCTTTGGTATGACGAGGGTGAAGAGAGAAACAGAGGAAACTTC
CCTCCTAGATTTTCAGGTCGCCAGTTCCCTAATTATAGCTCTGAGCTGAATGTGAACGCC
TTGGAGCTGGAGGACTCGGCCCTGTATCTCTGTGCCAGCAGCTTGG
>TRBV5-8*02
AGGACAGCAAGCGACTCTGAGATGCTCTCCTATCTCTGGGCACACCAGTGTGTACTGGTA
CCAACAGGCCCTGGGTCTGGGCCTCCAGCTCCTCCTTTGGTATGACGAGGGTGAAGAGAG
AAACAGAGGAAACTTCCCTCCTAGATTTTCAGGTCGCCAGTTCCCTAATTATAGCTCTGA
GCTGAATGTGAACGCCTTGGAGCTGGAGGACTCGGCCCTGTATCTCTGTGCCAGCAGC
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
>TRBJ1-6*01
CTCCTATAATTCACCCCTCCACTTTGGGAATGGGACCAGGCTCACTGTGACAG
>TRBJ1-6*02
CTCCTATAATTCACCCCTCCACTTTGGGAACGGGACCAGGCTCACTGTGACAG
>TRBJ2-1*01
CTCCTACAATGAGCAGTTCTTCGGGCCAGGGACACGGCTCACCGTGCTAG
>TRBJ2-2*01
CGAACACCGGGGAGCTGTTTTTTGGAGAAGGCTCTAGGCTGACCGTACTGG
>TRBJ2-2P*01
CTGAGAGGCGCTGCTGGGCGTCTGGGCGGAGGACTCCTGGTTCTGG
>TRBJ2-3*01
AGCACAGATACGCAGTATTTTGGCCCAGGCACCCGGCTGACAGTGCTCG
>TRBJ2-4*01
AGCCAAAAACATTCAGTACTTCGGCGCCGGGACCCGGCTCTCAGTGCTGG
>TRBJ2-5*01
ACCAAGAGACCCAGTACTTCGGGCCAGGCACGCGGCTCCTGGTGCTCG
>TRBJ2-6*01
CTCTGGGGCCAACGTCCTGACTTTCGGGGCCGGCAGCAGGCTGACCGTGCTGG
>TRBJ2-7*01
CTCCTACGAGCAGTACTTCGGGCCGGGCACCAGGCTCACGGTCACAG
>TRBJ2-7*02
CTCCTACGAGCAGTACGTCGGGCCGGGCACCAGGCTCACGGTCACAG
"""
str_mock_VDJ_fln_V_gene_CDR3_anchors = \
"""gene;anchor_index;gfunction
TRBV1*01;267;P
TRBV2*01;273;F
TRBV2*02;273;(F)
TRBV2*03;273;(F)
TRBV3-1*01;270;F
TRBV3-1*02;270;(F)
TRBV3-2*01;270;P
TRBV3-2*02;270;P
TRBV3-2*03;270;(P)
TRBV4-1*01;270;F
TRBV4-1*02;243;(F)
TRBV4-2*01;270;F
TRBV4-2*02;270;(F)
TRBV4-3*01;270;F
TRBV4-3*02;270;(F)
TRBV4-3*03;270;(F)
TRBV4-3*04;219;(F)
TRBV5-1*01;270;F
TRBV5-1*02;270;(F)
TRBV5-2*01;259;P
TRBV5-3*01;270;ORF
TRBV5-3*02;270;ORF
TRBV5-4*01;270;F
TRBV5-4*02;270;(F)
TRBV5-4*03;222;(F)
TRBV5-4*04;180;(F)
TRBV5-5*01;270;F
TRBV5-5*02;270;(F)
TRBV5-5*03;270;(F)
TRBV5-6*01;270;F
TRBV5-7*01;270;ORF
TRBV5-8*01;270;F
TRBV5-8*02;226;(F)
"""
str_mock_VDJ_fln_J_gene_CDR3_anchors = \
"""gene;anchor_index;function
TRBJ1-1*01;17;F
TRBJ1-2*01;17;F
TRBJ1-3*01;19;F
TRBJ1-4*01;20;F
TRBJ1-5*01;19;F
TRBJ1-6*01;22;F
TRBJ1-6*02;22;F
TRBJ2-1*01;19;F
TRBJ2-2*01;20;F
TRBJ2-3*01;18;F
TRBJ2-4*01;19;F
TRBJ2-5*01;17;F
TRBJ2-6*01;22;F
TRBJ2-7*01;16;F
"""

str_mock_VDJ_fln_dict = dict()
str_mock_VDJ_fln_dict['fln_genomicVs'] = str_mock_VDJ_fln_genomicVs
str_mock_VDJ_fln_dict['fln_genomicDs'] = str_mock_VDJ_fln_genomicDs
str_mock_VDJ_fln_dict['fln_genomicJs'] = str_mock_VDJ_fln_genomicJs
str_mock_VDJ_fln_dict['fln_V_gene_CDR3_anchors'] = str_mock_VDJ_fln_V_gene_CDR3_anchors
str_mock_VDJ_fln_dict['fln_J_gene_CDR3_anchors'] = str_mock_VDJ_fln_J_gene_CDR3_anchors


class MyTestCase(unittest.TestCase):
    def setUp(self) -> None:
        self.species = "human"
        self.chain = "tcr_beta"
        self.mdl = IgorModel.load_default(self.species, self.chain)
        rcParams["paths.igor_exec"] = "/home/olivares/GitHub/statbiophys/IGoR/digor"
        # rcParams["paths.igor_exec"] = None
        # mdl = IgorTask.run_generate(return_df=True)
        self.pd_sequences = generate(10, self.mdl)

    def test_evaluate(self):
        # 0. Get your input sequences, in this case generated sequences
        print(self.pd_sequences)
        self.assertIsInstance(self.pd_sequences, pd.DataFrame)

        # 2. Use a Model to evaluate sequences
        hb_mdl = IgorModel.load_default("human", "tcr_beta")
        self.assertIsInstance(hb_mdl, IgorModel)

        # 3. infer a new model using the initial  model.
        evaluated_seqs = evaluate_pgen(self.pd_sequences, hb_mdl)
        print(evaluated_seqs)

    def test_from_generation_to_evaluation(self):
        hb_mdl = get_default_IgorModel("human", "tcr_beta")
        self.assertIsInstance(hb_mdl, IgorModel)
        sequences = generate(10, hb_mdl)
        self.assertIsInstance(sequences, pd.DataFrame)

        ref_genome = get_IgorRefGenome_VDJ_from_IMGT("Homo+sapiens", "TRB")
        self.assertIsInstance(ref_genome, IgorRefGenome)
        ref_genome.clean_empty_anchors()
        self.assertIsInstance(ref_genome, IgorRefGenome)

        mdl_ini = IgorModel.make_default_model_from_IgorRefGenome(ref_genome)
        self.assertIsInstance(mdl_ini, IgorModel)

        mdl_new, df_likelihood = infer(sequences, mdl_ini)
        self.assertIsInstance(mdl_new, IgorModel)
        self.assertIsInstance(df_likelihood, pd.DataFrame)

    def test_evaluate_sequences(self):
        hb_mdl = get_default_IgorModel("human", "tcr_beta")
        sequences = generate(10, hb_mdl)
        df = evaluate(sequences, hb_mdl)
        print(df)


    def test_one_sequence_evaluate(self):
        hb_mdl = get_default_IgorModel("human", "tcr_beta")
        self.assertIsInstance(hb_mdl, IgorModel)
        sequences = generate(10, hb_mdl)
        self.assertIsInstance(sequences, pd.DataFrame)
        ref_genome = get_IgorRefGenome_VDJ_from_IMGT("Homo+sapiens", "TRB")
        self.assertIsInstance(ref_genome, IgorRefGenome)
        ref_genome.clean_empty_anchors()
        self.assertIsInstance(ref_genome, IgorRefGenome)
        mdl_ini = IgorModel.make_default_model_from_IgorRefGenome(ref_genome)
        self.assertIsInstance(mdl_ini, IgorModel)
        mdl_new, df_likelihood = infer(sequences, mdl_ini)
        one_sequence = 'GGTGCTGTCGTCTCTCAACATCCGAGCTGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTCAGTTCGCGGGAATTAGGAACACTGAAGCTTTCTTTGGACAAGGCACCAGACTCACAGTTGTAG'
        one_pgen = evaluate(one_sequence, mdl_new, N_scenarios=5)
        print(one_pgen)

    def test_one_sequence_list_evaluate(self):
        hb_mdl = get_default_IgorModel("human", "tcr_beta")
        one_sequence = 'GGTGCTGTCGTCTCTCAACATCCGAGCTGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTCAGTTCGCGGGAATTAGGAACACTGAAGCTTTCTTTGGACAAGGCACCAGACTCACAGTTGTAG'
        tuplita = (7, one_sequence)
        one_pgen = evaluate(tuplita , hb_mdl, N_scenarios=5)
        print(one_pgen)

    def test_one_sequence_scenarios_output(self):
        hb_mdl = get_default_IgorModel("human", "tcr_beta")
        one_sequence = 'GGTGCTGTCGTCTCTCAACATCCGAGCTGGGTTATCTGTAAGAGTGGAACCTCTGTGAAGATCGAGTGCCGTTCCCTGGACTTTCAGGCCACAACTATGTTTTGGTATCGTCAGTTCCCGAAACAGAGTCTCATGCTGATGGCAACTTCCAATGAGGGCTCCAAGGCCACATACGAGCAAGGCGTCGAGAAGGACAAGTTTCTCATCAACCATGCAAGCCTGACCTTGTCCACTCTGACAGTGACCAGTGCCCATCCTGAAGACAGCAGCTTCTACATCTGCAGTGCTCAGTTCGCGGGAATTAGGAACACTGAAGCTTTCTTTGGACAAGGCACCAGACTCACAGTTGTAG'
        tuplita = (7, one_sequence)
        fln_output = 'df_scenarios.csv'
        df_scenarios = evaluate(tuplita, hb_mdl, N_scenarios=5, airr_format=False, fln_output=fln_output)
        self.assertIsInstance(df_scenarios, pd.DataFrame)

    def test_one_sequence_plot(self):
        hb_mdl = get_default_IgorModel("human", "tcr_beta")
        fln_output = 'df_scenarios.csv'
        df_scenarios = hb_mdl.get_dataframe_scenarios(fln_output)
        ps_scenario = df_scenarios.iloc[0]
        print(ps_scenario)

        hb_mdl.plot_scenario(ps_scenario)
        plt.show()

        ###########################################
        # hb_mdl.plot_scenario(ps_scenario)
        df_scenario_aln = hb_mdl.get_df_scenario_aln_from_scenario(ps_scenario)
        da_scenario_aln = from_df_scenario_aln_to_da_scenario_aln(df_scenario_aln)
        # fig, ax = plot_scenario_from_da_scenario_aln(da_scenario_aln)
        aln_scenario_np = da_scenario_aln.values

        cmap_dna = mpl.colors.ListedColormap(['white', '#fcff92', '#70f970', '#ff99b1', '#4eade1'])

        # fig, ax = plt.subplots(figsize=(40, 40))

        # ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        # ax.set_aspect('equal')

        fig, ax = plt.subplots(figsize=(20, 10))
        ax.imshow(aln_scenario_np, cmap=cmap_dna, vmin=-1.5, vmax=3.5)
        ###########################################
        plt.show()

        # hb_mdl.plot_scenario(ps_scenario)
        # plt.show()

        # fig, ax = plt.subplots(figsize=(20, 10))
        # ax.imshow(aln_scenario_np, cmap=cmap_dna, vmin=-1.5, vmax=3.5)
        # ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    def test_IgorTask_evaluate(self):
        task = IgorTask(mdl=self.mdl)
        pd_rearrangement = task.evaluate(self.pd_sequences, igor_wd='here', clean_batch=False)
        self.assertIsInstance(pd_rearrangement, pd.DataFrame)
        print(pd_rearrangement)

    def test_CDR3(self):
        hb_mdl = get_default_IgorModel("human", "tcr_beta")
        df_sequences = generate(10, hb_mdl, seed=0)
        df_scenarios, df_offsets = evaluate(df_sequences, hb_mdl, airr_format=False, b_V_offset=True)

        seq_index = 5
        ps_scenario = df_scenarios.loc[seq_index].iloc[0]

        V_offset = int(df_offsets.loc[seq_index].loc[ps_scenario[hb_mdl.event_GeneChoice_V_nickname]]) # None
        print("V_offset: ", V_offset)



        if V_offset is None:
            V_offset = 0
        len_scenario = 0
        len_GeneChoice = 0
        len_Insertion = 0
        len_Deletion = 0

        print("ps_scenario: ", ps_scenario)
        for event_nickname in hb_mdl.event_GeneChoice_nickname_list:
            len_GeneChoice += len(hb_mdl.realization(ps_scenario, event_nickname).value)
            print(event_nickname, len_GeneChoice)
        for event_nickname in hb_mdl.event_Insertion_nickname_list:
            len_Insertion += hb_mdl.realization(ps_scenario, event_nickname).value
        for event_nickname in hb_mdl.event_Deletion_nickname_list:
            len_Deletion += -hb_mdl.realization(ps_scenario, event_nickname).value

        len_scenario = len_GeneChoice + len_Insertion + len_Deletion
        print('*'*80)
        print("len_scenario: ", len_scenario)

        V_nickname = hb_mdl.event_GeneChoice_V_nickname
        J_nickname = hb_mdl.event_GeneChoice_J_nickname

        V_choice_realization = hb_mdl.realization(ps_scenario, V_nickname)
        J_choice_realization = hb_mdl.realization(ps_scenario, J_nickname)

        V_anchor_in_seq = V_offset + int( hb_mdl.V_anchor(V_choice_realization.id) )
        J_anchor_in_seq = V_offset + int( hb_mdl.J_anchor(J_choice_realization.id) ) - len(J_choice_realization.value) + len_scenario + 3

        V_anchor_in_seq, J_anchor_in_seq

        print("V_anchor_in_seq: ", V_anchor_in_seq)
        print('-'*40)
        print("J_anchor_in_seq: ", J_anchor_in_seq)
        print(df_sequences)
        print('-'*50)
        str_seq = df_sequences['nt_sequence'].loc[seq_index]
        print("str_seq: ", len(str_seq), str_seq)

        str_CDR3_nt = str_seq[V_anchor_in_seq:J_anchor_in_seq]
        print("str_CDR3_nt: ", str_CDR3_nt)
        print( dna_translate(str_CDR3_nt) )
        self.assertIsInstance(hb_mdl, IgorModel)
        # df_sequences = generate(10, hb_mdl)
        # self.assertIsInstance(df_sequences, pd.DataFrame)
        # df_scenarios, df_V_offsets = evaluate(df_sequences,hb_mdl,5, b_V_offset=True, airr_format=False)
        # self.assertIsInstance(df_scenarios, pd.DataFrame)
        # # self.assertIsInstance(df_V_offsets, pd.DataFrame)
        # aaa = hb_mdl.get_VDJ_CDR3_from_df_scenario(df_scenarios)
        # print(aaa)
        hb_mdl.plot_scenario(ps_scenario)
        plt.show()





if __name__ == '__main__':
    unittest.main()
