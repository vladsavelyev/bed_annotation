from contextlib import contextmanager
from traceback import format_exc
import yaml

from ngs_utils.logger import info, err, critical, debug, warn
from ngs_utils.file_utils import verify_file
from ngs_utils.utils import update_dict


def _load_yaml(fname):
    with open(fname) as in_handle:
        config = yaml.load(in_handle, Loader=yaml.FullLoader)
    return config


def load_yaml_config(fpath):
    verify_file(fpath, is_critical=True)
    try:
        dic = _load_yaml(fpath)
    except Exception:
        err(format_exc())
        critical('Could not parse bcbio YAML ' + fpath)
    else:
        return dic


def fill_dict_from_defaults(cur_cnf, defaults_dict):
    for key in defaults_dict:
        if key in cur_cnf:
            if isinstance(cur_cnf[key], dict) and isinstance(defaults_dict[key], dict):
                fill_dict_from_defaults(cur_cnf[key], defaults_dict[key])
        else:
            cur_cnf[key] = defaults_dict[key]
    return cur_cnf


def _join_parent_conf(child_conf, parent_conf):
    bc = parent_conf.copy()
    bc.update(child_conf)
    child_conf.update(bc)
    return child_conf


@contextmanager
def with_cnf(cnf, **kwargs):
    prev_opts = {k: cnf[k] for k in kwargs.keys()}
    try:
        for k, v in kwargs.items():
            cnf[k] = v
        yield cnf
    finally:
        for k, v in prev_opts.items():
            cnf[k] = v


def merge_config_files(fnames):
    """Merge configuration files, preferring definitions in latter files.
    """
    out = _load_yaml(fnames[0])
    for fname in fnames[1:]:
        update_dict(to_update=out, updates=_load_yaml(fname))
    return out

