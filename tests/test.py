import os
from collections import namedtuple
from datetime import datetime
from genericpath import getmtime
from nose import SkipTest
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath
from ngs_utils.testing import BaseTestCase, swap_output, check_call
from ngs_utils.file_utils import add_suffix


class AnnotateBedTests(BaseTestCase):
    script = 'annotate_bed.py'

    data_dir = join(dirname(__file__), BaseTestCase.data_dir)
    results_dir = join(dirname(__file__), BaseTestCase.results_dir)
    gold_standard_dir = join(join(dirname(__file__), BaseTestCase.gold_standard_dir))

    def test_hg19(self):
        self._test('hg19', 'hg19')

    def test_hg38(self):
        self._test('hg38', 'hg38')

    def test_grch37(self):
        self._test('GRCh37', 'GRCh37')

    def test_hg38_ext(self):
        self._test('hg38_ext', 'hg38', '--extended')

    def test_hg38_extended_features(self):
        self._test('hg38_ext_features', 'hg38', '--extended --output-features')

    def test_hg38_many_options(self):
        self._test('hg38_ext_features_high_conf', 'hg38', '--extended --output-features --high-confidence --canonical --coding-only')

    def test_hg38_ambiguities_best_one(self):
        self._test('hg38_short', 'GRCh37', '--ambiguities best_one')

    def test_hg38_ambiguities_best_all(self):
        self._test('hg38_short', 'GRCh37', '--ambiguities best_all')

    def test_hg38_ambiguities_all(self):
        self._test('hg38_short', 'GRCh37', '--ambiguities all')

    def test_hg38_short(self):
        self._test('hg38_short', 'hg38', '--short')

    def test_mm10(self):
        self._test('mm10', 'mm10', '--extended')

    def _test(self, name, genome, opts=None):
        os.chdir(self.results_dir)
        input_fname = genome + '.bed'
        output_fname = add_suffix(input_fname, 'anno')
        input_fpath = join(self.data_dir, input_fname)
        output_fpath = join(self.results_dir, output_fname)

        cmdl = f'{self.script} {input_fpath} -o {output_fpath} {opts} -g {genome}'

        swap_output(output_fpath)

        print('-' * 100)
        check_call(cmdl, shell=True)
        print('-' * 100)
        print('')

        self._check_file(output_fpath)

