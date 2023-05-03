import collections
import hashlib
import os
import random
from functools import reduce
from os import environ
import socket
import re
from collections import OrderedDict
import itertools


class OrderedDefaultDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        if not args:
            self.default_factory = None
        else:
            if not (args[0] is None or callable(args[0])):
                raise TypeError('first argument must be callable or None')
            self.default_factory = args[0]
            args = args[1:]
        super(OrderedDefaultDict, self).__init__(*args, **kwargs)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = default = self.default_factory()
        return default

    def __reduce__(self):  # optional, for pickle support
        args = (self.default_factory,) if self.default_factory else ()
        return self.__class__, args, None, None, self.items()


def _tryint(s):
    try:
        return int(s)
    except ValueError:
        return s


def _alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [_tryint(c) for c in re.split('([0-9]+)', s)]


def human_sorted(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=_alphanum_key)
    return l


def format_integer(name, value, unit=''):
    value = int(value)
    if value is not None:
        return '{name}: {value:,}{unit}'.format(**locals())
    else:
        return '{name}: -'.format(**locals())


def format_decimal(name, value, unit=''):
    if value is not None:
        return '{name}: {value:.2f}{unit}'.format(**locals())
    else:
        return '{name}: -'.format(**locals())


def mean(values):
    return float(sum(values)) / len(values) if len(values) > 0 else float('nan')


def median(values):
    values = sorted(values)

    if len(values) % 2 == 1:  # odd number of values
        return values[(len(values) - 1) // 2]
    else:  # even number of values - take the avg of central
        return (values[len(values) // 2] + values[len(values) // 2 - 1]) // 2


def get_numeric_value(string_value):
    """ parses string_value and returns only number-like part
    """
    num_chars = ['.', '+', '-']
    number = ''
    for c in string_value:
        if c.isdigit() or c in num_chars:
            number += c
    return number


def get_hostname():
    return os.environ.get('HOSTNAME') or os.environ.get('HOST') or socket.gethostname()

hostname = get_hostname()

def is_scp():
    return 'scp' in hostname

def is_us():
    return any(name in hostname for name in ['rask', 'chara', 'blue', 'green', 'espo',
                                             'orr', 'usbod', 'bn0', 'pedro', 'papi'])
def is_uk():
    return 'ukap' in hostname

def is_sweden():
    return 'seml' in hostname

def is_china():
    return 'cniclhpc' in hostname

def is_az():
    return is_scp() or is_us() or is_uk() or is_china() or is_sweden()

def is_cloud():
    return 'starcluster' in hostname

def is_cluster():
    return is_az() or is_cloud()

def is_travis():
    return 'TRAVIS' in environ

def is_chihua():
    return hostname == 'chihua'

def is_local():
    return 'Vlads' in hostname or 'Vladislavs' in hostname or 'local' in hostname or 'Home' in hostname or '5180L-135800-M.local' in hostname


def md5(fpath):
    hash = hashlib.md5()
    with open(fpath, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash.update(chunk)
    return hash.hexdigest()


def gray(text):
    return '<span class="gray">' + text + '</span>'


#########################
### Dict utils

def update_dict(to_update, updates):
    """ Recursively updates a nested dict `to_update` from a nested dict `updates`
    """
    for key, val in updates.items():
        if isinstance(val, collections.Mapping):
            to_update[key] = update_dict(to_update.get(key, {}), val)
        else:
            to_update[key] = updates[key]
    return to_update


def get_in(d, t, default=None):
    """
    look up if you can get a tuple of values from a nested dictionary,
    each item in the tuple a deeper layer

    example: get_in({1: {2: 3}}, (1, 2)) -> 3
    example: get_in({1: {2: 3}}, (2, 3)) -> {}
    """
    result = reduce(lambda d, t: d.get(t, {}), t, d)
    if not result:
        return default
    else:
        return result


##########################
### Functional programming

def partition_all(n, iterable):
    """Partition a list into equally sized pieces, including last smaller parts
    http://stackoverflow.com/questions/5129102/python-equivalent-to-clojures-partition-all
    """
    it = iter(iterable)
    while True:
        chunk = list(itertools.islice(it, n))
        if not chunk:
            break
        yield chunk


def partition(pred, iterable):
    'Use a predicate to partition entries into false entries and true entries'
    # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
    t1, t2 = itertools.tee(iterable)
    return itertools.filterfalse(pred, t1), filter(pred, t2)


def flatten(l):
    """
    Flatten an irregular list of lists
    example: flatten([[[1, 2, 3], [4, 5]], 6]) -> [1, 2, 3, 4, 5, 6]
    lifted from: http://stackoverflow.com/questions/2158395/
    """
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, str):
            for sub in flatten(el):
                yield sub
        else:
            yield el


def compose(f, g):
    return lambda x: f(g(x))


##########################
### Misc

def is_sequence(arg):
    """
    check if 'arg' is a sequence

    example: arg([]) -> True
    example: arg("lol") -> False

    """
    return (not hasattr(arg, "strip") and
            hasattr(arg, "__getitem__") or
            hasattr(arg, "__iter__"))


def is_pair(arg):
    """check if 'arg' is a two-item sequence
    """
    return is_sequence(arg) and len(arg) == 2


def is_string(arg):
    return isinstance(arg, str)


def reservoir_sample(stream, num_items, item_parser=lambda x: x):
    """
    samples num_items from the stream keeping each with equal probability
    """
    kept = []
    for index, item in enumerate(stream):
        if index < num_items:
            kept.append(item_parser(item))
        else:
            r = random.randint(0, index)
            if r < num_items:
                kept[r] = item_parser(item)
    return kept


def dictapply(d, fn):
    """
    apply a function to all non-dict values in a dictionary
    """
    for k, v in d.items():
        if isinstance(v, dict):
            v = dictapply(v, fn)
        else:
            d[k] = fn(v)
    return d


def set_locale():
    import locale
    try:
        if 'UTF-8' not in locale.getlocale(locale.LC_ALL):
            try:
                locale.setlocale(locale.LC_ALL, 'en_AU.UTF-8')
            except locale.Error:
                locale.setlocale(locale.LC_ALL, 'C.UTF-8')
            var = '.'.join(locale.getlocale(locale.LC_ALL))
            environ['LC_ALL'] = environ['LANG'] = var
    except TypeError:
        pass


