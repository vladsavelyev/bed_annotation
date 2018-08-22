import os
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath
from ngs_utils.testing import BaseTestCase, swap_output, check_call
from ngs_utils.file_utils import add_suffix


class AnnotateBedTests(BaseTestCase):
    script = 'bed_annotation'

    data_dir = join(dirname(__file__), BaseTestCase.data_dir)
    results_dir = join(dirname(__file__), BaseTestCase.results_dir)
    gold_standard_dir = join(join(dirname(__file__), BaseTestCase.gold_standard_dir))

    def test_hg19(self):
        self._test('hg19', 'hg19')

    def test_grch37_ext(self):
        self._test('grch37_ext', 'GRCh37', '--extended')

    def test_hg38_features_ext(self):
        self._test('hg38_features_ext', 'hg38', '--output-features --extended')

    def test_many_options(self):
        self._test('many_opts', 'hg38', '--extended --output-features --high-confidence --canonical --coding-only')

    def test_amb_best_one(self):
        self._test('amb_best_one', 'GRCh37', '--ambiguities best_one')

    def test_amb_best_all(self):
        self._test('amb_best_all', 'GRCh37', '--ambiguities best_all')

    def test_amb_all(self):
        self._test('amb_all', 'GRCh37', '--ambiguities all')

    def test_hg38_short(self):
        self._test('hg38_short', 'hg38', '--short')

    def test_mm10(self):
        self._test('mm10', 'mm10', '--extended')

    def _test(self, name, genome, opts=''):
        os.chdir(self.results_dir)
        input_fname = genome + '.bed'
        output_fname = add_suffix(input_fname, f'anno_{name}')
        input_fpath = join(self.data_dir, input_fname)
        output_fpath = join(self.results_dir, output_fname)

        cmdl = f'{self.script} {input_fpath} -o {output_fpath} {opts} -g {genome}'
        self._run_cmd(cmdl, [input_fpath], output_fpath)
        self._check_file_throws(output_fpath)

