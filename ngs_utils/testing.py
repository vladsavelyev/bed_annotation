import subprocess
import unittest
import os
from os.path import dirname, join, exists, isfile, splitext, basename, isdir, relpath, getctime, getsize, realpath, \
    islink
import sys
from datetime import datetime
from glob import glob

from ngs_utils.file_utils import verify_dir, verify_file, safe_mkdir
from ngs_utils.utils import is_travis


def echo(msg=''):
    sys.stderr.write(msg + '\n')

def call(cmdl, suppress_output=False):
    echo(cmdl if isinstance(cmdl, str) else subprocess.list2cmdline(cmdl))
    if isinstance(cmdl, str):
        return subprocess.call(cmdl, shell=True, executable='/bin/bash',
                               stdout=subprocess.DEVNULL if suppress_output else None)
    else:
        return subprocess.call(cmdl,
                               stdout=subprocess.DEVNULL if suppress_output else None)

def check_call(cmdl):
    echo(cmdl if isinstance(cmdl, str) else subprocess.list2cmdline(cmdl))
    if isinstance(cmdl, str):
        subprocess.check_call(cmdl, shell=True, executable='/bin/bash')
    else:
        subprocess.check_call(cmdl)

def check_output(cmdl):
    echo(cmdl if isinstance(cmdl, str) else subprocess.list2cmdline(cmdl))
    if isinstance(cmdl, str):
        return subprocess.check_output(cmdl, shell=True, executable='/bin/bash', stderr=subprocess.STDOUT)
    else:
        return subprocess.check_output(cmdl, stderr=subprocess.STDOUT)

def swap_output(output_path):
    if not exists(output_path):
        return

    last_changed = datetime.fromtimestamp(getctime(output_path))
    prev_output_path = output_path + '_' + last_changed.strftime('%Y_%m_%d__%H_%M_%S')
    os.rename(output_path, prev_output_path)

    # Adding _prev symlink
    swap_prev_symlink(output_path, prev_output_path)

    return prev_output_path

def swap_prev_symlink(output_path, prev_output_path=None):
    """
    :param output_path:
    :param prev_output_path: the older "previous" path that will be symlinked to <output_path>_prev
    """
    prev_output_link = get_prev(output_path)
    if exists(prev_output_link) or islink(prev_output_link):
        os.remove(prev_output_link)
    assert not exists(prev_output_link), f'{prev_output_link} still exists and cannot be removed'
    if prev_output_path:
        os.symlink(prev_output_path, prev_output_link)

def get_prev(fpath):
    return fpath + '_prev'


class BaseTestCase(unittest.TestCase):
    # Run on top of existing latest results. Also controlled with TEST_REUSE
    reuse = False
    # Do not run, just diff the latest results against the gold standard. Also controlled with TEST_ONLY_DIFF
    only_diff = False

    script = None

    data_dir = 'data'
    results_dir = 'results'
    gold_standard_dir = 'gold_standard'

    remove_work_dir_on_success = False

    def setUp(self):
        if not isdir(self.data_dir):
            os.makedirs(self.data_dir)
        if not exists(self.results_dir):
            os.makedirs(self.results_dir)

    def _run_cmd(self, cmdl, input_paths, output_path, before_run_fn=None):
        only_diff  = BaseTestCase.only_diff or any('TEST' in e.upper() and 'DIFF'  in e.upper() for e in os.environ)
        reuse      = BaseTestCase.reuse     or any('TEST' in e.upper() and 'REUSE' in e.upper() for e in os.environ)
        if only_diff:
            echo('TESRT_DIFF_ONLY set: not actually running the program, only checking diffs with the previous results')
        if reuse:
            echo('TEST_REUSE set: running on top of the previous results')
        tools_opts = next((os.environ[e] for e in os.environ if 'TEST' in e.upper() and 'OPTS' in e.upper()), '')

        if not only_diff:
            input_paths = [input_paths] if isinstance(input_paths, str) else input_paths
            for ip in input_paths:
                assert exists(ip), f'Data {ip} does not exist.'

            if not reuse:
                swap_output(output_path)
            safe_mkdir(dirname(output_path))

            if before_run_fn:
                before_run_fn()

            echo('-' * 100)
            check_call(cmdl + ' ' + tools_opts)
            echo('-' * 100)
            echo('')

    def _check_file_throws(self, wc_fpath, ignore_matching_lines=None, wrapper=None, cmp_line_number_only=True,
                           check_diff=True):
        found = glob(wc_fpath)
        assert found, 'file is not found ' + wc_fpath
        assert len(found) == 1, 'more than 1 file is found as ' + wc_fpath
        fpath = found[0]
        assert isfile(fpath), 'is not a file: ' + fpath
        assert getsize(fpath) > 0, 'file is empty: ' + fpath

        if isdir(self.gold_standard_dir):
            cmp_wc_path = join(self.gold_standard_dir, relpath(wc_fpath, self.results_dir))
        else:
            cmp_wc_path = get_prev(wc_fpath)

        cmp_found = glob(cmp_wc_path)
        assert cmp_found, 'cmp_file not found (no gold standard dir or *_prev file): ' + cmp_wc_path
        assert len(cmp_found) == 1, 'more than 1 cmp_file is found as ' + cmp_wc_path
        cmp_fpath = cmp_found[0]
        assert isfile(cmp_fpath), 'cmp_file is not a file: ' + cmp_fpath

        if check_diff:
            cmdl = 'diff'
            if ignore_matching_lines:
                if isinstance(ignore_matching_lines, str):
                    ignore_matching_lines = [ignore_matching_lines]
                for r in ignore_matching_lines:
                    cmdl += ' -I ' + subprocess.list2cmdline([r])
            if wrapper:
                if isinstance(wrapper, list):
                    wrapper = subprocess.list2cmdline(wrapper)
                if not fpath.endswith('.gz'):
                    fpath = '<(cat ' + fpath + ' | ' + wrapper + ')'
                    cmp_fpath = '<(cat ' + cmp_fpath + ' | ' + wrapper + ')'
                else:
                    fpath = '<(gunzip -c ' + fpath + ' | ' + wrapper + ')'
                    cmp_fpath = '<(gunzip -c ' + cmp_fpath + ' | ' + wrapper + ')'
            elif fpath.endswith('.gz'):
                fpath = '<(gunzip -c ' + fpath + ')'
                cmp_fpath = '<(gunzip -c ' + cmp_fpath + ')'
            cmdl += ' ' + fpath + ' ' + cmp_fpath
            ret_code = call(cmdl, suppress_output=not is_travis())
            assert ret_code == 0, 'diff returned non-zero: ' + fpath

    @staticmethod
    def _check_dir_not_empty(dirpath, description=None):
        assert verify_dir(dirpath, description=description), dirpath
        contents = [join(dirpath, fname) for fname in os.listdir(dirpath)
                    if not fname.startswith('.')]
        assert len(contents) >= 1, dirpath + ': ' + str(contents)
        assert all(verify_file(realpath(fpath), is_critical=True)
                   for fpath in contents
                   if isfile(realpath(fpath))), dirpath + ': ' + str(contents)

vcf_ignore_lines = [
    '^##bcftools_',
    '^##INFO=',
    '^##FILTER=',
    '^##contig=',
]

# Find and parse all elements containing json data, put data into a list and dumps the result.
# The resulting text is unique per json data, so we can run simple `diff` on them.
html_wrapper = [
    'grep', '-A1', '<div id=".*_json">', '|', 'grep', '-v', '<div id=".*_json">', '|',
    'python', '-c',
        'import sys, json; '
        'sys.stdout.write(json.dumps([json.loads(el) for el in sys.stdin.read().split(\'--\')], '
                                     'indent=2, sort_keys=True))'
]
