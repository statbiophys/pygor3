import unittest
from pygor3 import imgt
from pygor3 import IgorRefGenome
import tempfile

class MyTestCase(unittest.TestCase):
    def test_get_species_list(self):
        imgt_species_list = imgt.get_species_list()
        self.assertIsInstance(imgt_species_list, list)

    def test_download_ref_genome_VDJ(self):
        species = "Homo+sapiens"
        chain = "TRB"
        try:
            tmp_dir = tempfile.TemporaryDirectory(dir=".", prefix="models")
            type(tmp_dir.name)
            ref_genome_pd_dict = imgt.download_ref_genome_VDJ(species, chain, modelspath=tmp_dir.name)
            self.assertTrue('V' in ref_genome_pd_dict)
            ref_genome_pd_dict.pop('V', None)
            self.assertTrue('D' in ref_genome_pd_dict)
            ref_genome_pd_dict.pop('D', None)
            self.assertTrue('J' in ref_genome_pd_dict)
            ref_genome_pd_dict.pop('J', None)
        except Exception as e:
            raise e
        finally:
            tmp_dir.cleanup()

    def test_download_ref_genome_VJ(self):
        species = "Homo+sapiens"
        chain = "TRA"
        try:
            tmp_dir = tempfile.TemporaryDirectory(dir=".", prefix="models")
            type(tmp_dir.name)
            ref_genome_pd_dict = imgt.download_ref_genome_VJ(species, chain, modelspath=tmp_dir.name)
            self.assertTrue('V' in ref_genome_pd_dict)
            ref_genome_pd_dict.pop('V', None)
            self.assertTrue('J' in ref_genome_pd_dict)
            ref_genome_pd_dict.pop('J', None)
        except Exception as e:
            raise e
        finally:
            tmp_dir.cleanup()


if __name__ == '__main__':
    unittest.main()
