import unittest
from click.testing import CliRunner
from pygor3.scripts.cli import *
class MyTestCase(unittest.TestCase):

    def test_cli_demo_get_data(self):
        opciones = ["demo-get-data"]
        runner = CliRunner()
        result = runner.invoke(get_ref_genome, opciones)
        print(result.exit_code)
        # print(result.output)
        # print(result.exception)

    def test_cli_imgt_get_genomes_info(self):
        opciones = ["--info"]
        runner = CliRunner()
        result = runner.invoke(get_ref_genome, opciones)
        print(result.exit_code)
        # print(result.output)
        # print(result.exception)

    def test_cli_imgt_get_genomes_VJ(self):
        runner = CliRunner()
        result = runner.invoke(get_ref_genome, ["--imgt-species", "Homo+sapiens", "--imgt-chain", "IGL", "-t", "VJ"])
        # print(result.exit_code)
        # print(result.output)
        # print(result.exception)

    def test_cli_imgt_get_genomes_VDJ(self):
        opciones = ["--imgt-species", "Homo+sapiens", "--imgt-chain", "IGK", "-t", "VDJ"]
        runner = CliRunner()
        result = runner.invoke(get_ref_genome, opciones)
        # print(result.exit_code)
        # print(result.output)
        # print(result.exception)

    def test_cli_model_create_from_fasta_genes(self):
        opciones = ["--Vgene", "models/Homo+sapiens/IGL/ref_genome/genomicVs.fasta",
                    "--Jgene", "models/Homo+sapiens/IGL/ref_genome/genomicJs.fasta",
                    "-t", "VJ",
                    "-o", "nuevo"]
        runner = CliRunner()
        result = runner.invoke(model_create, opciones)
        # print(result.exit_code)
        # print(result.output)
        # print(result.exception)

    def test_cli_model_create(self):
        opciones = ["-M", "models/Homo+sapiens/IGL/", "-t", "VJ"]
        runner = CliRunner()
        result = runner.invoke(model_create, opciones)
        print(result.exit_code)
        print(result.output)
        print(result.exception)

    def test_cli_run_generate(self):
        "pygor igor-generate -s human -c beta -N 10"
        opciones = ["-s", "human", "-c", "beta", "-N", "10"]
        runner = CliRunner()
        result = runner.invoke(run_generate, opciones)
        print(result.exit_code)
        print(result.output)
        print(result.exception)

    def test_cli_run_infer(self):
        "pygor igor-infer -M models/Homo+sapiens/IGL/ -i demo/data/IgL/IgL_seqs_naive_Nofunctional.txt -o IgL_new"
        opciones = ["-M", "models/Homo+sapiens/IGL/", "-i", "demo/data/IgL/IgL_seqs_naive_Nofunctional.txt", "-o", "IgL_new"]
        runner = CliRunner()
        result = runner.invoke(run_infer, opciones)
        print(result.exit_code)
        print(result.output)
        print(result.exception)

    def test_cli_run_evaluate_m(self):
        "pygor igor-evaluate -m IgL_new_parms.txt IgL_new_marginals.txt --Vanchors IgL_new_V_gene_CDR3_anchors.csv --Janchors IgL_new_J_gene_CDR3_anchors.csv -i IgL_seq_memory_Functional_sample.txt -o IgL_evaluation"
        opciones = ["-m", "IgL_new_parms.txt", "IgL_new_marginals.txt",
                    "--Vanchors", "IgL_new_V_gene_CDR3_anchors.csv",
                    "--Janchors", "IgL_new_J_gene_CDR3_anchors.csv",
                    "-i", "IgL_seq_memory_Functional_sample.txt",
                    "-o", "IgL_evaluation"]
        runner = CliRunner()
        result = runner.invoke(run_evaluate, opciones)
        print(result.exit_code)
        print(result.output)
        print(result.exception)

    def test_cli_run_evaluate(self):
        "pygor igor-evaluate -s human -c beta -i joder_sequences.csv -o evaluados"
        opciones = ["-s", "human", "-c", "beta", "-i", "joder_sequences.csv", "-o", "evaluados"]
        runner = CliRunner()
        result = runner.invoke(run_evaluate, opciones)
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
