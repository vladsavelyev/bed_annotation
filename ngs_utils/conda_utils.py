from os.path import join, abspath, dirname, pardir, exists, isdir
import sys

from ngs_utils.logger import critical


def secondary_conda_env(env_name='pcgr', is_critical=False):
    py_path = sys.executable  # e.g. /miniconda/envs/umccrise/bin/python
    env_path = dirname(dirname(py_path))  # e.g. /miniconda/envs/umccrise
    env_path = env_path + '_' + env_name  # e.g. /miniconda/envs/umccrise_pcgr
    if not isdir(env_path):
        if is_critical:
            critical(f'Can\'t find environment {env_path}')
        else:
            return None
    return env_path

