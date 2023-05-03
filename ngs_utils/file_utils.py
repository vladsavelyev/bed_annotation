"""
Different utilities for dealing with files, directories, and their names.
A lot of functions are grabbed from https://github.com/chapmanb/bcbio-nextgen
"""

import shutil
import os
import string
from os.path import isfile, isdir, getsize, exists, basename, join, abspath, splitext, \
    islink, dirname, realpath, getmtime, getctime
import gzip
import tempfile
import contextlib
import fnmatch
import time

from ngs_utils.logger import info, err, warn, critical, debug
from ngs_utils.utils import is_sequence, is_string


def safe_mkdir(dirpath, descriptive_name=''):
    """ Multiprocessing-safely and recursively creates a directory
    """
    if not dirpath:
        critical(f'Path is empty: {descriptive_name if descriptive_name else ""}')

    if isdir(dirpath):
        return dirpath

    if isfile(dirpath):
        critical(descriptive_name + ' ' + dirpath + ' is a file.')

    num_tries = 0
    max_tries = 10

    while not exists(dirpath):
        # we could get an error here if multiple processes are creating
        # the directory at the same time. Grr, concurrency.
        try:
            os.makedirs(dirpath)
        except OSError as e:
            if num_tries > max_tries:
                raise
            num_tries += 1
            time.sleep(2)
    return dirpath


@contextlib.contextmanager
def chdir(new_dir):
    """Context manager to temporarily change to a new directory.

    http://lucentbeing.com/blog/context-managers-and-the-with-statement-in-python/
    """
    cur_dir = os.getcwd()
    safe_mkdir(new_dir)
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(cur_dir)


def symlink_plus(orig, new):
    """Create relative symlinks and handle associated biological index files.
    """
    for ext in ["", ".idx", ".gbi", ".tbi", ".bai"]:
        if os.path.exists(orig + ext) and not os.path.lexists(new + ext):
            with chdir(os.path.dirname(new)):
                os.symlink(os.path.relpath(orig + ext), os.path.basename(new + ext))
    orig_noext = splitext_plus(orig)[0]
    new_noext = splitext_plus(new)[0]
    for sub_ext in [".bai"]:
        if os.path.exists(orig_noext + sub_ext) and not os.path.lexists(new_noext + sub_ext):
            with chdir(os.path.dirname(new_noext)):
                os.symlink(os.path.relpath(orig_noext + sub_ext), os.path.basename(new_noext + sub_ext))


def open_gzipsafe(f, mode='r'):
    # mode_t = mode.replace('b', '')
    # mode_b = mode if 'b' in mode else mode + 'b'
    if f.endswith('.gz') or f.endswith('.gzip') or f.endswith('.gz.tx') or f.endswith('.gzip.tx'):
        try:
            h = gzip.open(f, mode=mode + 't', encoding='UTF-8')
        except IOError as e:
            err('Error opening gzip ' + f + ': ' + str(e) + ', opening as plain text')
            return open(f, mode=mode)
        else:
            if 'w' in mode:
                return h
            else:
                try:
                    h.read(1)
                except IOError as e:
                    err('Error opening gzip ' + f + ': ' + str(e) + ', opening as plain text')
                    h.close()
                    return open(f, mode=mode)
                else:
                    h.close()
                    h = gzip.open(f, mode=mode + 't')
                    return h
    else:
        return open(f, mode=mode)


def replace_suffix(to_transform, suffix):
    """
    replaces the suffix on a filename or list of filenames
    example: replace_suffix("/path/to/test.sam", ".bam") ->
    "/path/to/test.bam"
    """
    if is_sequence(to_transform):
        transformed = []
        for f in to_transform:
            base, _ = os.path.splitext(f)
            transformed.append(base + suffix)
        return transformed
    elif is_string(to_transform):
        base, _ = os.path.splitext(to_transform)
        return base + suffix
    else:
        raise ValueError("replace_suffix takes a single filename as a string or "
                         "a list of filenames to transform.")


def locate(pattern, root=os.curdir):
    """Locate all files matching supplied filename pattern in and below
    supplied root directory.
    """
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)


def replace_directory(out_files, dest_dir):
    """
    change the output directory to dest_dir
    can take a string (single file) or a list of files
    """
    if is_sequence(out_files):
        filenames = map(os.path.basename, out_files)
        return [os.path.join(dest_dir, x) for x in filenames]
    elif is_string(out_files):
        return os.path.join(dest_dir, os.path.basename(out_files))
    else:
        raise ValueError("in_files must either be a sequence of filenames "
                         "or a string")

def which(program):
    """
    returns the path to an executable or None if it can't be found
    """
    def is_exe(_fpath):
        return os.path.isfile(_fpath) and os.access(_fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def adjust_path(path):
    if path is None: return None

    path = remove_quotes(path)
    if path is None: return None

    path = expanduser(path)
    if path is None: return None

    path = abspath(path)
    if path is None: return None

    return path.replace('//', '/')


code_base_path = abspath(join(dirname(abspath(__file__)), os.path.pardir))


def adjust_system_path(path):
    if path is None: return None

    path = remove_quotes(path)
    if path is None: return None

    path = expanduser(path)
    if path is None: return None

    path = join(code_base_path, path)  # will only join if the tool_path is not absolute:
    if path is None: return None

    path = realpath(path)
    if path is None: return None

    path = abspath(path)
    if path is None: return None

    return path


# Expand paths beginning with '~' or '~user'.
# '~' means $HOME; '~user' means that user's home directory.
# If the path doesn't begin with '~', or if the user or $HOME is unknown,
# the path is returned unchanged (leaving error reporting to whatever
# function is called with the expanded path as argument).
# See also module 'glob' for expansion of *, ? and [...] in pathnames.
# (A function should also be defined to do full *sh-style environment
# variable expansion.)
def expanduser(path):
    """Expand ~ and ~user constructs.

    If user or $HOME is unknown, do nothing."""
    if path[:1] != '~':
        return path
    i, n = 1, len(path)
    while i < n and path[i] not in '/\\':
        i = i + 1

    if 'HOME' in os.environ:
        userhome = os.environ['HOME']
    elif 'USERPROFILE' in os.environ:
        userhome = os.environ['USERPROFILE']
    elif not 'HOMEPATH' in os.environ:
        return path
    else:
        try:
            drive = os.environ['HOMEDRIVE']
        except KeyError:
            drive = ''
        userhome = join(drive, os.environ['HOMEPATH'])

    if i != 1:  # ~user
        userhome = join(dirname(userhome), path[1:i])

    return userhome + path[i:]


def file_exists(fpath):
    """Check if a file exists and is non-empty.
    """
    return fpath and exists(adjust_path(fpath)) and getsize(adjust_path(fpath)) > 0


def _log(msg, silent, is_critical):
    if is_critical:
        critical(msg)
    if not silent:
        warn(msg)


def verify_obj_by_path(path, description='', silent=False, is_critical=False, verify_size=True):
    if path is None:
        msg = (description + ': i' if description else 'I') + 's not specified (None).'
        _log(msg, silent, is_critical)
        return path

    if not path:
        msg = (description + ': f' if description else 'N') + 'ame is empty.'
        _log(msg, silent, is_critical)
        return path

    path = adjust_path(path)
    if not exists(path):
        msg = (description + ': ' if description else '') + path + ' does not exist.'
        _log(msg, silent, is_critical)
        return None

    if isfile(path):
        return verify_file(path, description, silent, verify_size=verify_size)
    elif isdir(path):
        return verify_dir(path, description, silent)
    else:
        msg = (description + ': ' if description else '') + path + ' is not a file or a directory.'
        _log(msg, silent, is_critical)
        return None


def can_reuse(fpath, cmp_f, silent=False):
    """Check if a file `fpath` exists, is non-empty and is more recent than `cmp_f`
    """
    do_reuse = os.environ.get('REUSE', '1')
    if do_reuse == '0':
        return False
    if not fpath or not isfile(fpath):
        return False
    elif verify_file(fpath, cmp_f=cmp_f, silent=True):
        if not silent:
            debug('Reusing ' + fpath)
        return True
    else:
        return False


def verify_file(fpath, description='', silent=False, is_critical=False, verify_size=True, cmp_f=None):
    if fpath is None:
        msg = (description + ': ' if description else ' ') + 'not specified.'
        _log(msg, silent, is_critical)
        return fpath

    if not fpath:
        msg = (description + ': f' if description else 'F') + 'ile name is empty.'
        _log(msg, silent, is_critical)
        return fpath

    fpath = adjust_path(fpath)
    if not exists(fpath):
        msg = (description + ': ' if description else '') + fpath + ' does not exist.'
        _log(msg, silent, is_critical)
        return None

    if not isfile(fpath):
        msg = (description + ': ' if description else '') + fpath + ' is not a file.'
        _log(msg, silent, is_critical)
        return None

    if verify_size and getsize(fpath) <= 0:
        msg = (description + ': ' if description else '') + fpath + ' is empty.'
        _log(msg, silent, is_critical)
        return None

    if cmp_f:
        if isinstance(cmp_f, str):
            cmp_f = [cmp_f]
        try:
            for cmp_f in cmp_f:
                if cmp_f and getmtime(fpath) < getmtime(cmp_f):
                    msg = (description + ': ' if description else '') + fpath + ' is older than ' + cmp_f
                    _log(msg, silent, is_critical)
                    return None
        except OSError as e:
            _log(str(e), silent, is_critical)

    return fpath


def verify_dir(dirpath, description='', silent=False, is_critical=False):
    if dirpath is None:
        msg = (description + ': i' if description else 'I') + 's not specified (None).'
        _log(msg, silent, is_critical)
        return dirpath

    if not dirpath:
        msg = (description + ': d' if description else 'D') + 'ir name is empty.'
        _log(msg, silent, is_critical)
        return None

    dirpath = adjust_path(dirpath)
    if not exists(dirpath):
        msg = (description + ': ' if description else '') + dirpath + ' does not exist.'
        _log(msg, silent, is_critical)
        return None

    if not isdir(dirpath):
        msg = (description + ': ' if description else '') + dirpath + ' is not a directory.'
        _log(msg, silent, is_critical)
        return None

    return dirpath


def verify_module(name):
    try:
        __import__(name)
        return True
    except:
        return False


def num_lines(fpath):
    with open(fpath) as f:
        return sum(1 for _ in f)


def make_tmpdir():
    # base_dir = cnf.tmp_base_dir or cnf.work_dir or os.getcwd()
    # if not verify_dir(base_dir, 'Base directory for temporary files'):
    #     sys.exit(1)
    #
    return tempfile.mkdtemp()


@contextlib.contextmanager
def tmpdir():
    dirpath = make_tmpdir()
    try:
        yield dirpath
    finally:
        try:
            shutil.rmtree(dirpath)
        except OSError:
            warn('Warning: cannot clean up temporary dir ' + dirpath)


@contextlib.contextmanager
def workdir(cnf):
    if cnf.work_dir:
        verify_dir(cnf.work_dir, is_critical=True)
        yield cnf.work_dir
    else:
        cnf.work_dir = make_tmpdir()
        yield cnf.work_dir
        try:
            shutil.rmtree(cnf.work_dir)
        except OSError:
            warn('Warning: cannot clean up temporary dir ' + cnf.work_dir)


def make_tmpfile(work_dir, *args, **kwargs):
    """ Returns tuple (file descriptor and file path)
    """
    return tempfile.mkstemp(work_dir, *args, **kwargs)


@contextlib.contextmanager
def tmpfile(work_dir, *args, **kwargs):
    tmp_file, fpath = make_tmpfile(work_dir, *args, **kwargs)
    try:
        yield fpath
    finally:
        try:
            os.remove(fpath)
        except OSError:
            pass


def splitext_plus(fname):
    """Split on file extensions, allowing for zipped extensions.
    """
    base, ext = splitext(fname)
    if ext in [".gz", ".bz2", ".zip"]:
        base, ext2 = splitext(base)
        ext = ext2 + ext
    return base, ext


def add_suffix(fname, suf):
    base, ext = splitext_plus(fname)
    return base + (('.' + suf) if suf else '') + ext


def intermediate_fname(work_dir, fpath, suf):
    suf_fpath = add_suffix(fpath, suf)
    if work_dir:
        return join(work_dir, basename(suf_fpath))
    else:
        return suf_fpath


def remove_quotes(s):
    if s and s[0] in ['"', "'"]:
        s = s[1:]
    if s and s[-1] in ['"', "'"]:
        s = s[:-1]
    return s


def convert_file(work_dir, input_fpath, convert_file_fn, suffix=None, output_fpath=None,
                 check_result=True, overwrite=False, reuse=True, ctx=None):
    assert output_fpath or suffix, str(output_fpath) + ' ' + str(suffix)
    output_fpath = output_fpath or intermediate_fname(work_dir, input_fpath, suf=suffix)
    if output_fpath.endswith('.gz'):
        debug('output_fpath is .gz, but writing to uncompressed.')
        output_fpath = splitext(output_fpath)[0]
    
    if not overwrite:
        if can_reuse(output_fpath, cmp_f=input_fpath):
            debug('Reusing ' + output_fpath)
            return output_fpath
        if can_reuse(output_fpath + '.gz', cmp_f=input_fpath):
            debug('Reusing ' + output_fpath + '.gz')
            return output_fpath
    
    if islink(output_fpath):
        os.unlink(output_fpath)

    debug('Writing to ' + output_fpath)
    with file_transaction(work_dir, output_fpath) as tx_fpath:
        with open_gzipsafe(input_fpath) as inp_f, open(tx_fpath, 'w') as out_f:
            if ctx:
                convert_file_fn(inp_f, out_f, ctx)
            else:
                convert_file_fn(inp_f, out_f)

    if suffix or output_fpath:
        debug('Saved to ' + output_fpath)

    verify_file(output_fpath, is_critical=check_result)
    return output_fpath


def iterate_file(work_dir, input_fpath, proc_line_fun,
                 suffix=None, output_fpath=None, check_result=True, reuse=False, **kwargs):
    def _proc_file(inp_f, out_f, ctx=None):
        max_bunch_size = 1000 * 1000
        written_lines = 0
        bunch = []

        for i, line in enumerate(inp_f):
            clean_line = line.replace('\n', '')
            if clean_line:
                if ctx:
                    new_l = proc_line_fun(clean_line, i, ctx)
                else:
                    new_l = proc_line_fun(clean_line, i)
                if new_l is not None:
                    bunch.append(new_l + '\n')
                    written_lines += 1
            else:
                bunch.append(line)
                written_lines += 1

            if len(bunch) >= max_bunch_size:
                out_f.writelines(bunch)
                debug('Written lines: ' + str(written_lines))
                bunch = []

        out_f.writelines(bunch)
        debug('Written lines: ' + str(written_lines))

    return convert_file(work_dir, input_fpath, _proc_file,
        suffix=suffix, output_fpath=output_fpath, check_result=check_result, reuse=reuse, **kwargs)


def dots_to_empty_cells(config, tsv_fpath):
    """Put dots instead of empty cells in order to view TSV with column -t
    """
    def proc_line(l, i):
        while '\t\t' in l:
            l = l.replace('\t\t', '\t.\t')
        return l
    return iterate_file(config, tsv_fpath, proc_line, suffix='dots')


#################################################
######## Transaction ############################
@contextlib.contextmanager
def file_transaction(work_dir, *rollback_files):
    """Wrap file generation in a transaction, moving to output if finishes.
    """
    def __remove_files(fnames):
        if isinstance(fnames, str):
            fnames = [fnames]

        for x in fnames:
            if x and os.path.exists(x):
                if os.path.isfile(x):
                    os.remove(x)
                elif os.path.isdir(x):
                    shutil.rmtree(x, ignore_errors=True)

    def _flatten_plus_safe(tmp_dir, rollback_files):
        """Flatten names of files and create temporary file names.
        """
        tx_fpaths, orig_files = [], []
        for fnames in rollback_files:
            if isinstance(fnames, str):
                fnames = [fnames]
            for fname in fnames:
                tx_file = fname + '.tx'
                tx_fpath = join(tmp_dir, tx_file) if tmp_dir else tx_file
                tx_fpaths.append(tx_fpath)
                orig_files.append(fname)
        return tx_fpaths, orig_files

    exts = {".vcf": ".idx", ".bam": ".bai", "vcf.gz": ".tbi"}
    safe_fpaths, orig_names = _flatten_plus_safe(work_dir, rollback_files)
    __remove_files(safe_fpaths)  # remove any half-finished transactions
    try:
        if len(safe_fpaths) == 1:
            yield safe_fpaths[0]
        else:
            yield tuple(safe_fpaths)
    except:  # failure -- delete any temporary files
        __remove_files(safe_fpaths)
        raise
    else:  # worked -- move the temporary files to permanent location
        for safe, orig in zip(safe_fpaths, orig_names):
            if exists(safe):
                shutil.move(safe, orig)
                for check_ext, check_idx in exts.items():
                    if safe.endswith(check_ext):
                        safe_idx = safe + check_idx
                        if exists(safe_idx):
                            shutil.move(safe_idx, orig + check_idx)


@contextlib.contextmanager
def tx_tmpdir(base_dir, rollback_dirpath):
    """Context manager to create and remove a transactional temporary directory.
    """
    # tmp_dir_base = join(base_dir, 'tx', str(uuid.uuid4()))
    # unique_attempts = 0
    # while os.path.exists(tmp_dir_base):
    #     if unique_attempts > 5:
    #         break
    #     tmp_dir_base = join(base_dir, 'tx', str(uuid.uuid4()))
    #     time.sleep(1)
    #     unique_attempts += 1

    # if base_dir is not None:
    #     tmp_dir_base = os.path.join(base_dir, "tx")
    # else:
    #     tmp_dir_base = os.path.join(os.getcwd(), "tx")
    if exists(rollback_dirpath):
        critical(rollback_dirpath + ' already exists')

    tmp_dir = tempfile.mkdtemp(dir=base_dir)
    safe_mkdir(tmp_dir)
    try:
        yield tmp_dir
    finally:
        if tmp_dir and exists(tmp_dir):
            os.rename(tmp_dir, rollback_dirpath)


def safe_symlink_to(fpath, dst_dirpath, rel=False):
    if rel:
        fpath = os.path.relpath(fpath, dst_dirpath)

    dst = join(dst_dirpath, basename(fpath))
    if not exists(dst):
        try:
            if os.lstat(dst):  # broken symlink
                os.remove(dst)
        except OSError:
            pass
        debug('Symlink ' + fpath + ' -> ' + dst)
        os.symlink(fpath, dst)
    return dst


def safe_symlink(src_path, dst_path, rel=False):
    if rel:
        src_path = os.path.relpath(src_path, dirname(dst_path))

    if not exists(dst_path):
        try:
            if os.lstat(dst_path):  # broken symlink
                os.remove(dst_path)
        except:
            pass
        os.symlink(src_path, dst_path)
    return dst_path


def is_gz(fpath, mode='rb'):
    try:
        h = gzip.open(fpath)
    except IOError as e:
        return False
    else:
        try:
            h.read(1)
        except IOError as e:
            h.close()
            return False
        else:
            h.close()
            return True


def str_to_filename(s):
    valid_chars = "-_.%s%s" % (string.ascii_letters, string.digits)
    s = ''.join([c if c in valid_chars else '_' for c in s])
    return s


def get_ungz_gz(fpath):
    if fpath.endswith('.gz'):
        ungz = splitext(fpath)[0]
        gz = fpath
    else:
        ungz = fpath
        gz = fpath + '.gz'
    return ungz, gz
