import unittest
from click.testing import CliRunner
from pygor3.scripts.cli import *
class MyTestCase(unittest.TestCase):
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
