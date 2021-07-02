import unittest
from click.testing import CliRunner
from pygor3.scripts.cli import *
class MyTestCase(unittest.TestCase):

    def test_imgt_get_genomes_info(self):
        opciones = ["--info"]
        runner = CliRunner()
        result = runner.invoke(get_ref_genome, opciones)
        print(result.exit_code)
        # print(result.output)
        # print(result.exception)

    def test_imgt_get_genomes_VJ(self):
        runner = CliRunner()
        result = runner.invoke(get_ref_genome, ["--imgt-species", "Homo+sapiens", "--imgt-chain", "IGL", "-t", "VJ"])
        # print(result.exit_code)
        # print(result.output)
        # print(result.exception)

    def test_imgt_get_genomes_VDJ(self):
        opciones = ["--imgt-species", "Homo+sapiens", "--imgt-chain", "IGK", "-t", "VDJ"]
        runner = CliRunner()
        result = runner.invoke(get_ref_genome, opciones)
        # print(result.exit_code)
        # print(result.output)
        # print(result.exception)

    def test_model_create_from_fasta_genes(self):
        opciones = ["--Vgene", "models/Homo+sapiens/IGL/ref_genome/genomicVs.fasta",
                    "--Jgene", "models/Homo+sapiens/IGL/ref_genome/genomicJs.fasta",
                    "-t", "VJ",
                    "-o", "nuevo"]
        runner = CliRunner()
        result = runner.invoke(model_create, opciones)
        # print(result.exit_code)
        # print(result.output)
        # print(result.exception)

    def test_model_create(self):
        opciones = ["-M", "models/Homo+sapiens/IGL/", "-t", "VJ"]
        runner = CliRunner()
        result = runner.invoke(model_create, opciones)
        print(result.exit_code)
        print(result.output)
        print(result.exception)

    def test_something(self):
        runner = CliRunner()
        result = runner.invoke(run_generate, ["-w", "human", "-c", "tcr_beta", "-N", "10"])
        # run_generate(10, "human", "tcr_beta")
        print(result.exit_code)
        print(result.output)
        print(result.exception)
        # print(result.stderr)
        self.assertEqual(result.exit_code, 0)


if __name__ == '__main__':
    unittest.main()
