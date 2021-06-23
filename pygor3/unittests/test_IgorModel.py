import unittest
import tempfile
from pygor3 import IgorTask, IgorModel, IgorModel_Parms, IgorEvent_realization, IgorRefGenome
from pygor3 import *
import time
import subprocess
import copy

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
        self.tmp_dir = tempfile.TemporaryDirectory(dir=".", prefix="models")
        self.fln_dict = get_default_fln_names_for_model_dir(self.tmp_dir.name)
        # print(self.fln_dict)

        # Write a ref genome to use in other tests.
        self.ref_genome_path_dir = self.tmp_dir.name + "/ref_genome"
        os.makedirs(self.ref_genome_path_dir)

        self.model_path_dir = self.tmp_dir.name + "/models"
        os.makedirs(self.model_path_dir)

        for fln_key in self.fln_dict.keys():
            if fln_key in str_mock_VDJ_fln_dict.keys():
                with open(self.fln_dict[fln_key], mode='w') as ofile:
                    ofile.write(str_mock_VDJ_fln_dict[fln_key])


    def test_check_files(self):
        for fln_key in str_mock_VDJ_fln_dict.keys():
            self.assertTrue(os.path.isfile(self.fln_dict[fln_key]))


    def test_IgorModel_default(self):
        species = "human"
        chain = "tcr_beta"
        mdl = IgorModel.load_default(species, chain)
        self.assertIsInstance(mdl, IgorModel)
        self.assertIsInstance(mdl.parms, IgorModel_Parms)
        self.assertIsInstance(mdl.marginals, IgorModel_Marginals)

    def test_IgorModel_anchors(self):
        # TODO: ADD ANCHORS TO IGORMODEL
        species = "human"
        chain  = "tcr_beta"
        mdl = IgorModel.load_default(species, chain)
        print("+" * 40)
        print("mdl.genomic_dataframe_dict['V']: ", mdl.genomic_dataframe_dict['V'])
        print("-"*40)
        print("mdl.genomic_dataframe_dict['D']: ", mdl.genomic_dataframe_dict['D'])
        print("-" * 40)
        print("mdl.genomic_dataframe_dict['J']: ", mdl.genomic_dataframe_dict['J'])
        print("-" * 40)
        print("mdl.anchors_CDR3_V: ", mdl.anchors_CDR3_V)
        print("-" * 40)
        print("mdl.V_anchors: ", mdl.V_anchors)

    def test_IgorModel_anchors_from_genomic_dict(self):
        species = "human"
        chain = "tcr_beta"
        mdl = IgorModel.load_default(species, chain)
        mdl.parms.ErrorRate_dict
        # modify genomic_dataframe_dict
        print(mdl.genomic_dataframe_dict['V'])
        dictio = copy.deepcopy(mdl.genomic_dataframe_dict)
        mdl.set_genomic_dataframe_dict(dictio)
        mdl.generate_xdata()
        print("9"*10)
        print(mdl.genomic_dataframe_dict['V'])

    def test_IgorModel_make_default_VDJ(self):
        pass

    def test_IgorModel_edit_model_parms(self):

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

        mdl_hb.generate_xdata()
        mdl_hb.generate_Pmarginals()


        mdl_0 = IgorModel.make_default_VDJ(df_V, genomic_dict['D'], genomic_dict['J'])
        mdl_0.export_csv()
        mdl_hb.export_plot_events(fln_output_prefix + "_CP")

        self.assertIsInstance(mdl_0, IgorModel)

    def test_abalabdada(self):
        pass
        IgorModel.make_default_VDJ
        import copy
        mdl_copy = copy.deepcopy(mdl_hb)
        mdl_copy.genomic_dataframe_dict['V']['name'] = p3.v_genLabel(mdl_copy.genomic_dataframe_dict['V']['name'])
        mdl_copy.genomic_dataframe_dict['J']['name'] = p3.v_genLabel(mdl_copy.genomic_dataframe_dict['J']['name'])
        mdl_copy.genomic_dataframe_dict
        mdl_hb.write_model('model_parms.txt', 'model_marginals.txt',
                           fln_V_gene_CDR3_anchors="V_gene_CDR3_anchors.csv",
                           fln_J_gene_CDR3_anchors="J_gene_CDR3_anchors.csv")
        dictio = mdl_hb.genomic_dataframe_dict
        p3.v_genLabel(mdl_copy.parms['v_choice']['name'])
        mdl_copy.generate_xdata()
        mdl_copy['d_gene']
        mdl_copy.parms['j_choice']['name'] = p3.v_genLabel(mdl_copy.parms['j_choice']['name'])
        mdl_copy.parms['v_choice'], mdl_copy.parms['j_choice']
        mdl_copy.parms.df_V_anchors
        mdl_copy.genomic_dataframe_dict['V'][['name', 'anchor_index']].rename(columns={'name': 'gene'}).set_index(
            'gene').dropna()
        mdl_copy.genomic_dataframe_dict['V']
        mdl_copy.parms['v_choice']['name'] = p3.v_genLabel(mdl_copy.parms['v_choice']['name'])
        mdl_copy.parms['j_choice']['name'] = p3.v_genLabel(mdl_copy.parms['j_choice']['name'])
        help(mdl_copy.parms.set_event_realizations_from_DataFrame)
        # mdl_copy.marginals.marginals_dict['v_choice']
        mdl_copy.parms.set_event_realizations_from_DataFrame('v_choice', mdl_copy.parms['v_choice'])
        mdl_copy.parms.set_event_realizations_from_DataFrame('j_choice', mdl_copy.parms['j_choice'])

        help(mdl_copy.get_df_realizations)






        # mdl.parms.get_IgorRefGenome()
        # mdl.generate_xdata()
        # mdl.anchors_CDR3_V
        #
        # # To access the anchors use mdl.V_anchors
        # # this method should be use when mdl.V_anchors are
        # print("mdl.V_anchors: ", mdl.V_anchors)
        # # Make a copy of anchors in functions.
        # # and accessed by mdl.parms.df_V_anchors
        # print("mdl.parms.df_V_anchors: ", mdl.parms.df_V_anchors)
        #
        # print("mdl.anchors_CDR3_V: ", mdl.anchors_CDR3_V)

        # how to populate genomic_dataframe_dict

        # 1. Use parms.df_V_anchors as a copy to mdl.anchors_CDR3_V

        # print("mdl.parms.df_V_ref_genome: ", mdl.parms.df_V_ref_genome)
        # TODO: ADD genomic_dataframe_dict whenv generate_xdata is call
        # print("mdl.V_anchors: ", mdl.V_anchors)

        # 2. to attach mdl.anchors_CDR3_V to make a copy


        # print(mdl.V_anchors)

    def test_IgorModel_from_dataframes(self):
        ref_genome = IgorRefGenome.load_from_path(self.ref_genome_path_dir)
        # Because the model depends has VDJ genes
        self.assertIsInstance(ref_genome, IgorRefGenome)

        mdl_from_ref_genome = IgorModel.make_default_model_from_IgorRefGenome(ref_genome)
        self.assertIsInstance(mdl_from_ref_genome, IgorModel)

        print("mdl_from_ref_genome.V_anchors: ", mdl_from_ref_genome.V_anchors)
        print("mdl_from_ref_genome.J_anchors: ", mdl_from_ref_genome.J_anchors)

        path_mdl_data = self.tmp_dir.name + "/batch_mdldata"
        aaa = path_mdl_data + "/ref_genome"
        mdl_from_ref_genome
        mdl_from_ref_genome.write_mdldata_dir()

        ref_genome_again = mdl_from_ref_genome.parms.get_IgorRefGenome()
        self.assertIsInstance(ref_genome_again, IgorRefGenome)

        os.system("mkdir -p " + aaa)
        ref_genome_again.write_ref_genome_dir(aaa)
        fln_dict_tmp = get_default_ref_genome_fln_paths(ref_genome_path=aaa)
        print(fln_dict_tmp, str_mock_VDJ_fln_dict.keys())
        for fln_key in str_mock_VDJ_fln_dict.keys():
            print("fln_dict_tmp[" + fln_key + "]:", fln_dict_tmp[fln_key])
            self.assertTrue(os.path.isfile(fln_dict_tmp[fln_key]))

        """
        No V genes event found!
        No J genes event found!
        ==================================================
        []
        No J genes event found!
        None
        """

    def test_IgorModel_write_model(self):
        mdl_hb = IgorModel.load_default("human", "tcr_beta")
        fln_model_parms = 'model_parms.txt'
        fln_model_marginals = 'model_marginals.txt'
        fln_V_gene_CDR3_anchors = 'V_gene_CDR3_anchors.csv'
        fln_J_gene_CDR3_anchors = 'J_gene_CDR3_anchors.csv'

        ## TODO: ADD anchors
        mdl_hb.write_model(fln_model_parms, fln_model_marginals, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors)
        mdl_hb_2 = IgorModel(fln_model_parms, fln_model_marginals,
                             fln_V_gene_CDR3_anchors=fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors=fln_J_gene_CDR3_anchors)

        # TODO: CHANGE ANCHORS WITH DATAFRAME.
        mdl_hb_2

        # FIXME: SOLVE THE REFERENCES TO ANCHORS AND SEQUENCES NAMES PROBLEM
        # mdl_hb.parms.df_V_ref_genome
        print("mdl_hb_2.parms.dictNameNickname: ", mdl_hb_2.parms.dictNameNickname)
        df = mdl_hb_2.get_event_realizations_DataFrame('j_choice')

        print("mdl_hb_2.genomic_dataframe_dict: ", mdl_hb_2.genomic_dataframe_dict)
        new_df = df[:4]
        mdl_hb_2.parms.gen_NameNickname_dict()

        mdl_hb_2.set_realization_event_from_DataFrame('j_choice', new_df)
        mdl_hb_2.set_event_realizations_from_DataFrame('j_choice', new_df)
        print("mdl_hb_2.parms.dictNameNickname: ", mdl_hb_2.parms.dictNameNickname)
        mdl_copy.set_genomic_dataframe_dict()

        mdl_parms = IgorModel_Parms()
        mdl_parms.Event_list  # add_event() #add_Event()

        # import copy
        # mdl_copy = copy.deepcopy(mdl_hb)
        # mdl_copy

        # mdl_hb.parms.df_V_ref_genome
        # v_genLabel(mdl_hb.parms.df_V_ref_genome['name'])

        # print(mdl_hb_2.parms.df_V_ref_genome)
        # mdl_hb_2.genomic_dataframe_dict


        self.assertIsInstance(mdl_hb_2, IgorModel)
        self.assertTrue(os.path.isfile(fln_model_parms))
        self.assertTrue(os.path.isfile(fln_model_marginals))
        self.assertTrue(os.path.isfile(fln_V_gene_CDR3_anchors))
        self.assertTrue(os.path.isfile(fln_J_gene_CDR3_anchors))
        cmd = "rm {} {} {} {}".format(fln_model_parms, fln_model_marginals, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors)
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True)



    def test_IgorModel_name_change(self):
        pass

    """
    def test_IgorModel_default_from_system(self):
        species = "human"
        chain = "tcr_beta"
        mdl_hb = IgorModel.load_default(species, chain)
        mdl_hb.parms.Event_list[2].to_dict()
        mdl_hb.parms.df_V_ref_genome
        mdl_hb.genomic_dataframe_dict
        self.assertEqual(True, True)

    def test_IgorModel_default_VDJ(self):
        IgorModel.make_default_VDJ()
        IgorModel.make_model_default_VDJ_from_dataframes()

        species = "human"
        chain = "tcr_beta"
        mdl_hb = IgorModel.load_default(species, chain)

        mdl_hb.parms.Event_dict['d_gene']
        mdl_hb.parms.df_V_ref_genome
        mdl_hb.parms.df_D_ref_genome
        mdl_hb.parms.df_J_ref_genome
        mdl_hb.parms.Event_dict['v_choice']

        write_genetemplate_dataframe_to_fasta("v_genes.fasta", mdl_hb.parms.df_V_ref_genome)
        write_geneanchors_dataframe_to_csv("v_anchors.csv", mdl_hb.parms.df_V_ref_genome)
        task = IgorTask.default_model(species, "beta")
        task.fln_V_gene_CDR3_anchors
        task.fln_J_gene_CDR3_anchors
        mdl_hb.parms.attach_anchors(task.fln_V_gene_CDR3_anchors, task.fln_J_gene_CDR3_anchors)
        genomes = IgorRefGenome.load_default(species, chain)
        genomes.dict_genomicVs
        genomes.df_V_ref_genome
        aaa = genomes.get_anchors_dict()
        aaa
        df_genetemplates = mdl_hb.parms.get_Event('v_choice').get_realization_DataFrame()
        df_anchors = pd.read_csv(genomes.fln_V_gene_CDR3_anchors, sep=';')
        dddd = df_anchors.set_index('gene')
        dddd.set_index('gene')
        df_V_ref_genome = join_genomics_anchors_dataframes(df_genetemplates, df_anchors)
        df_V_ref_genome['anchor_index'].dropna()
        mdl_hb.parms.Event_dict['d_gene']

        df_fff
        df_fff['anchor_index']

        GeneChoice_list = [event for event in mdl_hb.parms.Event_list if event.event_type == 'GeneChoice']
        try:
            event_D = [event for event in GeneChoice_list if event.seq_type == 'D_gene'][0]
        except IndexError:
            print("WARNING: No D GeneChoice event found.")
            pass
        except Exception as e:
            raise e



        mdl_hb.parms.Event_dict['v_choice']
        # Gene templates from parms


        df_genetemplates = df_genetemplates[['name', 'value']]
        df_genetemplates['id'] = df_genetemplates.index.get_level_values('id')
        df_genetemplates.set_index('name', inplace=True)


        df_anchors

        df_all = df_genetemplates.join(df_anchors)
        df_all['name'] = df_all.index.get_level_values('name')
        df_all.set_index('id', inplace=True)
        df_all

        columnas = df_all.columns.to_list()
        ini_cols = ['name', 'value', 'anchor_index']
        other_cols = list()
        for col in columnas:
            if not col in ini_cols:
                other_cols.append(col)
        new_order = ini_cols + other_cols
        new_order
        df_all = df_all[new_order]
        type(df_all['anchor_index'])
        anchor_dict = df_all['anchor_index'].dropna().to_dict()
        anchor_dict[2]

        print( df_genetemplates.columns, len(df_genetemplates))

        dffff = df_genetemplates.set_index('name')
        df2222 = dffff.drop(['TRBV1*01'])
        df2222
        df2222.reset_index()



        print(df_anchors.columns, len(df_anchors))

        pd.set_option('display.max_columns', None)
        df_all = df_genetemplates.join(df_anchors)
        print(df_all.columns, len(df_all))
        print(df_all.head())
        df_all.set_index('id')

        df_anchors.loc['TRBV17*01']
        genomes.df_V_ref_genome.columns



        df_genetemplates['anchor_index'] = np.nan


        df_genetemplates.columns

        mdl_hb.parms.df_V_ref_genome #get_Event('v_choice').get_realization_DataFrame()

        mdl_hb.parms.df_V_ref_genome




        genomes.df_V_ref_genome['anchor_index']

        mdl_hb.anchors_CDR3_J[0]
        mdl_hb.parms.Event_list[0].to_dict()
        mdl_hb.parms.df_V_ref_genome

        df_genetemplates['anchor_index'] =
        event.to_dict()
        event.export_realizations_to_fasta()


        # Create a temporary directory with the genomic references tu run IGoR
        # directly from a model

        try:
            genomes.df_V_ref_genome.to_csv('joder.csv', index=False, columns=['name', 'anchor_index', 'function'], sep=';')
        except KeyError:
            print("KeyError")
            genomes.df_V_ref_genome.to_csv('joder.csv', index=False, columns=['name', 'anchor_index'],
                                           sep=';')
            pass
        except Exception as e:
            print(e)
            pass
        tmp_dir = tempfile.TemporaryDirectory(prefix='genomic', dir='.')
        genomes.update_fln_names(tmp_dir.name)
        genomes.write_ref_genome()
        tmp_dir.cleanup()

        mdl = IgorModel.make_model_default_VDJ_from_dataframes()
        dicto = genomes.get_anchors_dict()
        dicto['V']
        da = mdl.xdata['v_3_del']
        mdl.genomic_dataframe_dict
        print(da)
        # Make a temporal file
        tmp_file = tempfile.NamedTemporaryFile(prefix='genomic', dir='.', delete=True)
        # write IgorRefGenome in IGoR's files
        ref = IgorRefGenome()
        # write_dataframe_to_fasta(fln_fasta, df_genomic)
        ref.write_ref_genome()
        parms = IgorModel_Parms()
        mdl = IgorModel()
        # Use previous method to make a IgorRefGenome from IgorModel_Parms


        self.assertEqual(True, True)
        # lo que quiero es que
        # mdl.get_realizations_dict_from_scenario_dict()
        # realiz = IgorEvent_realization()
        # len(realiz.value)
        # mdl.v_anchor[realiz.id]
        # task = IgorTask.default_model("human", "beta")
        mdl = IgorModel_Parms.get_Event_list_sorted()

    def test_something(self):
        # TODO: MOVE THIS TO test IgorModel
        # 1. Download genomic templates
        # 2. Remove genes without anchors
        # 3. Create a default model with the templates
        # 4. Use some sequences to infer a model and return a new model
        #   4.1 Create a temporal directory to infer a model
        #   4.2 Run IGoR to infer
        #   4.3 Return mdl and also the likelihood log.
        #   4.4 Remove directory used to infer
        # 5. Then evaluate model

        # 3.
        tmp_dir = tempfile.TemporaryDirectory(dir=".", prefix="models")
        ref_genome_pd_dict = imgt.download_ref_genome_VDJ("Homo+sapiens", "TRB", modelspath=tmp_dir.name)
        # ref_genome_pd_dict = imgt.download_ref_genome_VJ("Homo+sapiens", "TRA", modelspath=tmp_dir.name)


        ref_genome = IgorRefGenome.load_from_dataframe_genomics_dict(ref_genome_pd_dict)
        ref_genome.clean_empty_anchors()
        ref_genome_dict = ref_genome.to_dict()
        mdl = IgorModel.make_default_from_Dataframe_dict(ref_genome_dict)
        mdl.parms.Event_dict
        mdl.V_anchors
        mdl.J_anchors
        mdl.genomic_dataframe_dict
        mdl.genomic_dataframe_dict['V']
        mdl.parms.Event_dict['v_choice']
        tmp_dir_Model_dir = tempfile.TemporaryDirectory(dir=".", prefix="igor_dir_model")
        tmp_dir_Model_dir.name
        # dictionary of default_paths given


        ref_genome_dir = "ref_genome/"
        fln_genomicVs = ref_genome_dir + "genomicVs.fasta"
        fln_genomicDs = ref_genome_dir + "genomicDs.fasta"
        fln_genomicJs = ref_genome_dir + "genomicJs.fasta"

        fln_V_gene_CDR3_anchors = ref_genome_dir + "V_gene_CDR3_anchors.csv"
        fln_J_gene_CDR3_anchors = ref_genome_dir + "J_gene_CDR3_anchors.csv"
        mdl.write_ref_genome(fln_genomicVs, fln_genomicDs, fln_genomicJs, fln_V_gene_CDR3_anchors, fln_J_gene_CDR3_anchors)

        ref_genome_files_dict = get_default_ref_genomes_species_chain("human", "tcr_alpha", ref_genome_path=".")
        ref_genome_files_dict
        mdl.genomic_dataframe_dict
        mdl.write_ref_genome(**ref_genome_files_dict)
        mdl.write_model()
        mdl.write_model()

        tmp_dir_Model_dir.cleanup()

        mdl = IgorModel.load_default("human", "tcr_beta")
        mdl.anchors_CDR3_J

        # FIXME: WHEN DEFAULT IS CALL mdl.parms.df_V_ref_genome should be fill if anchors are present
        #  also anchors can be add it manually.
        #

        mdl.parms.df_V_ref_genome
        self.assertEqual(True, False)

    """

    def tearDown(self) -> None:
        self.tmp_dir.cleanup()

if __name__ == '__main__':
    unittest.main()
