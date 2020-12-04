"""Common functions and classes used in multiple places in the MDTF code.
Specifically, util.py implements general functionality that's not MDTF-specific.
"""
from __future__ import absolute_import, division, print_function, unicode_literals
import os
import io
from src import six
import collections
import copy
import dataclasses
from distutils.spawn import find_executable
import enum
import errno
import functools
import glob
import json
import re
import shlex
import shutil
import signal
if os.name == 'posix' and six.PY2:
    try:
        import subprocess32 as subprocess
    except ImportError:
        import subprocess
else:
    import subprocess
import threading
import traceback
import typing
import unittest.mock
from six.moves import getcwd, collections_abc

class _Singleton(type):
    """Private metaclass that creates a :class:`~util.Singleton` base class when
    called. This version is copied from `<https://stackoverflow.com/a/6798042>`__ and
    should be compatible with both Python 2 and 3.
    """
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(_Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

class Singleton(_Singleton(six.ensure_str('SingletonMeta'), (object,), {})): 
    """Parent class defining the 
    `Singleton <https://en.wikipedia.org/wiki/Singleton_pattern>`_ pattern. We
    use this as safer way to pass around global state.
    """
    @classmethod
    def _reset(cls):
        """Private method of all :class:`~util.Singleton`-derived classes added
        for use in unit testing only. Calling this method on test teardown 
        deletes the instance, so that tests coming afterward will initialize the 
        :class:`~util.Singleton` correctly, instead of getting the state set 
        during previous tests.
        """
        # pylint: disable=maybe-no-member
        if cls in cls._instances:
            del cls._instances[cls]


class ExceptionPropagatingThread(threading.Thread):
    """Class to propagate exceptions raised in a child thread back to the caller
    thread when the child is join()ed. 
    Adapted from `<https://stackoverflow.com/a/31614591>`__.
    """
    def run(self):
        self.ret = None
        self.exc = None
        try:
            if hasattr(self, '_Thread__target'):
                # Thread uses name mangling prior to Python 3.
                self.ret = self._Thread__target(*self._Thread__args, **self._Thread__kwargs)
            else:
                self.ret = self._target(*self._args, **self._kwargs)
        except BaseException as e:
            self.exc = e

    def join(self, timeout=None):
        super(ExceptionPropagatingThread, self).join(timeout)
        if self.exc:
            raise self.exc
        return self.ret


class MultiMap(collections.defaultdict):
    """Extension of the :obj:`dict` class that allows doing dictionary lookups 
    from either keys or values. 
    
    Syntax for lookup from keys is unchanged, ``bd['key'] = 'val'``, while lookup
    from values is done on the `inverse` attribute and returns a set of matching
    keys if more than one match is present: ``bd.inverse['val'] = ['key1', 'key2']``.    
    See `<https://stackoverflow.com/a/21894086>`__.
    """
    def __init__(self, *args, **kwargs):
        """Initialize :class:`~util.MultiMap` by passing an ordinary :py:obj:`dict`.
        """
        super(MultiMap, self).__init__(set, *args, **kwargs)
        for key in iter(self.keys()):
            super(MultiMap, self).__setitem__(key, coerce_to_iter(self[key], set))

    def __setitem__(self, key, value):
        super(MultiMap, self).__setitem__(key, coerce_to_iter(value, set))

    def get_(self, key):
        if key not in list(self.keys()):
            raise KeyError(key)
        return coerce_from_iter(self[key])
    
    def to_dict(self):
        d = {}
        for key in iter(self.keys()):
            d[key] = self.get_(key)
        return d

    def inverse(self):
        d = collections.defaultdict(set)
        for key, val_set in iter(self.items()):
            for v in val_set:
                d[v].add(key)
        return dict(d)

    def inverse_get_(self, val):
        # don't raise keyerror if empty; could be appropriate result
        inv_lookup = self.inverse()
        return coerce_from_iter(inv_lookup[val])


class NameSpace(dict):
    """ A dictionary that provides attribute-style access.

    For example, `d['key'] = value` becomes `d.key = value`. All methods of 
    :py:obj:`dict` are supported.

    Note: recursive access (`d.key.subkey`, as in C-style languages) is not
        supported.

    Implementation is based on `<https://github.com/Infinidat/munch>`__.
    """

    # only called if k not found in normal places
    def __getattr__(self, k):
        """ Gets key if it exists, otherwise throws AttributeError.
            nb. __getattr__ is only called if key is not found in normal places.
        """
        try:
            # Throws exception if not in prototype chain
            return object.__getattribute__(self, k)
        except AttributeError:
            try:
                return self[k]
            except KeyError:
                raise AttributeError(k)

    def __setattr__(self, k, v):
        """ Sets attribute k if it exists, otherwise sets key k. A KeyError
            raised by set-item (only likely if you subclass NameSpace) will
            propagate as an AttributeError instead.
        """
        try:
            # Throws exception if not in prototype chain
            object.__getattribute__(self, k)
        except AttributeError:
            try:
                self[k] = v
            except Exception:
                raise AttributeError(k)
        else:
            object.__setattr__(self, k, v)

    def __delattr__(self, k):
        """ Deletes attribute k if it exists, otherwise deletes key k. A KeyError
            raised by deleting the key--such as when the key is missing--will
            propagate as an AttributeError instead.
        """
        try:
            # Throws exception if not in prototype chain
            object.__getattribute__(self, k)
        except AttributeError:
            try:
                del self[k]
            except KeyError:
                raise AttributeError(k)
        else:
            object.__delattr__(self, k)

    def __dir__(self):
        return list(self.keys())
    __members__ = __dir__  # for python2.x compatibility

    def __repr__(self):
        """ Invertible* string-form of a Munch.
            (*) Invertible so long as collection contents are each repr-invertible.
        """
        return '{0}({1})'.format(self.__class__.__name__, dict.__repr__(self))

    def __getstate__(self):
        """ Implement a serializable interface used for pickling.
        See `<https://docs.python.org/3.6/library/pickle.html>`__.
        """
        return {k: v for k, v in iter(self.items())}

    def __setstate__(self, state):
        """ Implement a serializable interface used for pickling.
        See `<https://docs.python.org/3.6/library/pickle.html>`__.
        """
        self.clear()
        self.update(state)

    def toDict(self):
        """ Recursively converts a NameSpace back into a dictionary.
        """
        return type(self)._toDict(self)

    @classmethod
    def _toDict(cls, x):
        """ Recursively converts a NameSpace back into a dictionary.
            nb. As dicts are not hashable, they cannot be nested in sets/frozensets.
        """
        if isinstance(x, dict):
            return dict((k, cls._toDict(v)) for k, v in iter(x.items()))
        elif isinstance(x, (list, tuple)):
            return type(x)(cls._toDict(v) for v in x)
        else:
            return x

    @property
    def __dict__(self):
        return self.toDict()

    @classmethod
    def fromDict(cls, x):
        """ Recursively transforms a dictionary into a NameSpace via copy.
            nb. As dicts are not hashable, they cannot be nested in sets/frozensets.
        """
        if isinstance(x, dict):
            return cls((k, cls.fromDict(v)) for k, v in iter(x.items()))
        elif isinstance(x, (list, tuple)):
            return type(x)(cls.fromDict(v) for v in x)
        else:
            return x

    def copy(self):
        return type(self).fromDict(self)
    __copy__ = copy

    def _freeze(self):
        """Return immutable representation of (current) attributes.

        We do this to enable comparison of two Namespaces, which otherwise would 
        be done by the default method of testing if the two objects refer to the
        same location in memory.
        See `<https://stackoverflow.com/a/45170549>`__.
        """
        d = self.toDict()
        d2 = {k: repr(d[k]) for k in d}
        FrozenNameSpace = collections.namedtuple(
            'FrozenNameSpace', sorted(list(d.keys()))
        )
        return FrozenNameSpace(**d2)

    def __eq__(self, other):
        if type(other) is type(self):
            return (self._freeze() == other._freeze())
        else:
            return False

    def __ne__(self, other):
        return (not self.__eq__(other)) # more foolproof

    def __hash__(self):
        return hash(self._freeze())

class MDTFEnum(enum.Enum):
    """Customize :py:class:`~enum.Enum`. 1) Assign (integer) values automatically
    to the members of the enumeration. 2) Provide a ``from_struct`` method to 
    simplify instantiating an instance from a string. To avoid potential 
    confusion with reserved keywords, we use the Python convention that members
    of the enumeration are all uppercase.
    """
    def __new__(cls, *args, **kwargs):
        """AutoNumber recipe from python stdlib docs."""
        value = len(cls.__members__) + 1
        obj = object.__new__(cls)
        obj._value_ = value
        return obj

    def __str__(self):
        return str(self.name).lower()

    def __repr__(self):
        return '<%s.%s>' % (self.__class__.__name__, self.name)

    @classmethod
    def from_struct(cls, str_):
        """Instantiate from string."""
        return cls.__members__.get(str_.upper())

NOTSET = unittest.mock.sentinel.NotSet
NOTSET.__doc__ = """
Sentinel object to detect uninitialized values, in cases where ``None`` is a 
valid value. For implentation, see `python docs 
<https://docs.python.org/3/library/unittest.mock.html#unittest.mock.sentinel>`__.
"""

MANDATORY = unittest.mock.sentinel.Mandatory
MANDATORY.__doc__ = """
Sentinel object to mark :func:`mdtf_dataclass` fields that do not take a default 
value. This is a workaround to avoid errors with non-default fields coming after
default fields in the dataclass-generated ``__init__`` method under 
`inheritance <https://docs.python.org/3/library/dataclasses.html#inheritance>`__:
we use the second solution described in `https://stackoverflow.com/a/53085935`__.

For implentation, see `python docs 
<https://docs.python.org/3/library/unittest.mock.html#unittest.mock.sentinel>`__.
"""

# declaration to allow calling with and without args: python cookbook 9.6
# https://github.com/dabeaz/python-cookbook/blob/master/src/9/defining_a_decorator_that_takes_an_optional_argument/example.py
def mdtf_dataclass(cls=None, **deco_kwargs):
    """Wrap :py:func:`~dataclasses.dataclass` class decorator to customize
    dataclasses to provide (very) rudimentary type checking and conversion. This
    is hacky, since dataclasses don't enforce type annontations for their fields.
    A better solution would be to use a deserialization library like pydantic.

    After the auto-generated ``__init__`` and the class' ``__post_init__``, the
    following tasks are performed:

    1. Verify that mandatory fields have values specified. We have to work around
       the usual :py:func:`~dataclasses.dataclass` way of doing this, because it 
       leads to errors in the signature of the dataclass-generated ``__init__`` 
       method under inheritance (mandatory fields can't come after optional 
       fields.) Mandatory fields must be designated by setting their default to
       ``MANDATORY``, and a ValueError is raised here if mandatory fields are
       uninitialized.

    2. Check each field's value to see if it's consistent with known type info. 
       If not, attempt to coerce it to that type, using a ``from_struct`` method if
       it exists. Raise ValueError if this fails.

    .. warning::
       Unlike :py:func:`~dataclasses.dataclass`, all fields **must** have a 
       *default* or *default_factory* defined. Fields which are mandatory must 
       have their default value set to the sentinel object ``MANDATORY``.

    .. warning::
       Type checking logic used is specific to the ``typing`` module in python 
       3.7. It may or may not work on newer pythons, and definitely will not 
       work with 3.5 or 3.6. See `https://stackoverflow.com/a/52664522`__.
    """
    dc_kwargs = {'init': True, 'repr': True, 'eq': True, 'order': False, 
        'unsafe_hash': False, 'frozen': False}
    dc_kwargs.update(deco_kwargs)
    if cls is None:
        # called without arguments
        return functools.partial(mdtf_dataclass, **dc_kwargs)

    cls = dataclasses.dataclass(cls, **dc_kwargs)
    _old_init = cls.__init__

    @functools.wraps(_old_init)
    def _new_init(self, *args, **kwargs):
        # Execute dataclass' auto-generated __init__ and __post_init__:
        _old_init(self, *args, **kwargs)
        
        for f in dataclasses.fields(self):
            if not f.init:
                # ignore fields that aren't handled at init
                continue
            value = getattr(self, f.name)
            # ignore unset field values, regardless of type
            if value is None or value is NOTSET:
                continue
            if value is MANDATORY:
                raise ValueError((f"{self.__class__.__name__}: No value supplied "
                    f"for mandatory field {f.name}."))
            # guess what types are valid
            new_type = None
            if f.type is typing.Any or isinstance(f.type, typing.TypeVar):
                continue
            elif isinstance(f.type, typing._GenericAlias) \
                or isinstance(f.type, typing._SpecialForm):
                # type is a generic from typing module, eg "typing.List"
                if f.type.__origin__ is typing.Union:
                    new_type = None # can't do coercion, but can test type
                    valid_types = list(f.type.__args__)
                elif issubclass(f.type.__origin__, typing.Generic):
                    continue # can't do anything in this case
                else:
                    new_type = f.type.__origin__
                    valid_types = [new_type]
            else:
                new_type = f.type
                valid_types = [new_type]
            # Get types of field's default value, if present. Dataclass doesn't 
            # require defaults to be same type as what's given for field.
            if not isinstance(f.default, dataclasses._MISSING_TYPE):
                valid_types.append(type(f.default))
            if not isinstance(f.default_factory, dataclasses._MISSING_TYPE):
                valid_types.append(type(f.default_factory()))
            
            try:
                if isinstance(value, tuple(valid_types)):
                    continue
                if new_type is None or hasattr(new_type, '__abstract_methods__'):
                    continue
                    # # can't do type coercion, so print a warning
                    # print((f"\tWarning: {self.__class__.__name__}: type of "
                    #     f" {f.name} is ({f.type}), recieved {repr(value)} of "
                    #     "conflicting type."))
                else:
                    # https://stackoverflow.com/a/54119384 for implementation
                    if hasattr(new_type, 'from_struct'):
                        object.__setattr__(self, f.name, new_type.from_struct(value))
                    else:
                        object.__setattr__(self, f.name, new_type(value))
            except (TypeError, ValueError, dataclasses.FrozenInstanceError) as exc: 
                print(exc)
                raise TypeError((f"{self.__class__.__name__}: Expected {f.name} "
                    f"to be {f.type}, got {type(value)} ({repr(value)}).")) from exc

    cls.__init__ = _new_init
    return cls

def mdtf_dataclass_factory(class_name, *parents, **kwargs):
    """Function that returns an mdtf_dataclass whose fields are the union of 
    the fields specified in its parent classes.

    Args:
        class_name: name of the new class.
        parents: collection of other mdtf_dataclasses to inherit from. Order in
            the collection determines the MRO.
        kwargs: arguments to pass to the mdtf_dataclass decorator for the
            returned class.
    """ 
    def _to_dataclass(self, cls_, **kwargs_):
        f"""Method to create an instance of one of the parent classes of
        {class_name} by copying over the relevant subset of fields.
        """
        new_kwargs = filter_dataclass(self, cls_)
        new_kwargs.update(kwargs_)
        return cls_(**new_kwargs)

    def _from_dataclasses(cls_, *other_dcs, **kwargs_):
        f"""Classmethod to create a new instance of {class_name} from instances
        of its parents, along with any other field values passed in kwargs.
        """
        new_kwargs = dict()
        for dc in other_dcs:
            new_kwargs.update(filter_dataclass(dc, cls_))
        new_kwargs.update(kwargs_)
        return cls_(**new_kwargs)

    methods = {
        'to_dataclass': _to_dataclass,
        'from_dataclasses': classmethod(_from_dataclasses),
    }
    for dc in parents:
        method_nm = 'to_' + dc.__name__.lower()
        methods[method_nm] = functools.partialmethod(_to_dataclass, cls_=dc)
    new_cls = type(class_name, tuple(parents), methods)
    return mdtf_dataclass(new_cls, **kwargs)

class ExceptionQueue(object):
    """Class to retain information about exceptions that were raised, for later
    output.
    """
    def __init__(self):
        self._queue = []

    @property
    def is_empty(self):
        return (len(self._queue) == 0)

    def log(self, exc, exc_to_chain=None):
        wrapped_exc = traceback.TracebackException.from_exception(exc)
        self._queue.append(wrapped_exc)

    def format(self):
        strs_ = [''.join(exc.format()) for exc in self._queue]
        strs_ = [f"***** Caught exception #{i+1}:\n{exc}\n" \
            for i, exc in enumerate(strs_)]
        return "".join(strs_)

# ------------------------------------

def strip_comments(str_, delimiter=None):
    # would be better to use shlex, but that doesn't support multi-character
    # comment delimiters like '//'
    if not delimiter:
        return str_
    s = str_.splitlines()
    for i in list(range(len(s))):
        if s[i].startswith(delimiter):
            s[i] = ''
            continue
        # If delimiter appears quoted in a string, don't want to treat it as
        # a comment. So for each occurrence of delimiter, count number of 
        # "s to its left and only truncate when that's an even number.
        # TODO: handle ' as well as ", for non-JSON applications
        s_parts = s[i].split(delimiter)
        s_counts = [ss.count('"') for ss in s_parts]
        j = 1
        while sum(s_counts[:j]) % 2 != 0:
            j += 1
        s[i] = delimiter.join(s_parts[:j])
    # join lines, stripping blank lines
    return '\n'.join([ss for ss in s if (ss and not ss.isspace())])

def read_json(file_path):
    assert os.path.exists(file_path), \
        "Couldn't find JSON file {}.".format(file_path)
    try:    
        with io.open(file_path, 'r', encoding='utf-8') as file_:
            str_ = file_.read()
    except IOError:
        print('Fatal IOError when trying to read {}. Exiting.'.format(file_path))
        exit()
    return parse_json(str_)

def parse_json(str_):
    str_ = strip_comments(str_, delimiter= '//') # JSONC quasi-standard
    try:
        parsed_json = json.loads(str_, object_pairs_hook=collections.OrderedDict)
    except UnicodeDecodeError:
        print('{} contains non-ascii characters. Exiting.'.format(str_))
        exit()
    return parsed_json

def write_json(struct, file_path, verbose=0, sort_keys=False):
    """Wrapping file I/O simplifies unit testing.

    Args:
        struct (:py:obj:`dict`)
        file_path (:py:obj:`str`): path of the JSON file to write.
        verbose (:py:obj:`int`, optional): Logging verbosity level. Default 0.
    """
    try:
        str_ = json.dumps(struct, 
            sort_keys=sort_keys, indent=2, separators=(',', ': '))
        with io.open(file_path, 'w', encoding='utf-8') as file_:
            file_.write(six.ensure_text(str_, encoding='utf-8', errors='strict'))
    except IOError:
        print('Fatal IOError when trying to write {}. Exiting.'.format(file_path))
        exit()

def pretty_print_json(struct, sort_keys=False):
    """Pseudo-YAML output for human-readable debugging output only - 
    not valid JSON"""
    str_ = json.dumps(struct, sort_keys=sort_keys, indent=2)
    for char in ['"', ',', '}', '[', ']']:
        str_ = str_.replace(char, '')
    str_ = re.sub(r"{\s+", "- ", str_)
    # remove lines containing only whitespace
    return os.linesep.join([s for s in str_.splitlines() if s.strip()]) 

def find_files(src_dirs, filename_globs):
    """Return list of files in `src_dirs` matching any of `filename_globs`. 

    Wraps glob.glob for the use cases encountered in cleaning up POD output.

    Args:
        src_dirs: Directory, or a list of directories, to search for files in.
            The function will also search all subdirectories.
        filename_globs: Glob, or a list of globs, for filenames to match. This 
            is a shell globbing pattern, not a full regex.

    Returns: :py:obj:`list` of paths to files matching any of the criteria.
        If no files are found, the list is empty.
    """
    src_dirs = coerce_to_iter(src_dirs)
    filename_globs = coerce_to_iter(filename_globs)
    files = set([])
    for d in src_dirs:
        for g in filename_globs:
            files.update(glob.glob(os.path.join(d, g)))
            files.update(glob.glob(os.path.join(d, '**', g), recursive=True))
    return list(files)

def recursive_copy(src_files, src_root, dest_root, copy_function=None, 
    overwrite=False):
    """Copy src_files to dest_root, preserving relative subdirectory structure.

    Copies a subset of files in a directory subtree rooted at src_root to an
    identical subtree structure rooted at dest_root, creating any subdirectories
    as needed. For example, `recursive_copy('/A/B/C.txt', '/A', '/D')` will 
    first create the destination subdirectory `/D/B` and copy '/A/B/C.txt` to 
    `/D/B/C.txt`.

    Args:
        src_files: Absolute path, or list of absolute paths, to files to copy.
        src_root: Root subtree of all files in src_files. Raises a ValueError
            if all files in src_files are not contained in the src_root directory.
        dest_root: Destination directory in which to create the copied subtree.
        copy_function: Function to use to copy individual files. Must take two 
            arguments, the source and destination paths, respectively. Defaults 
            to :py:meth:`shutil.copy2`.
        overwrite: Boolean, deafult False. If False, raise an OSError if
            any destination files already exist, otherwise silently overwrite.
    """
    if copy_function is None:
        copy_function = shutil.copy2
    src_files = coerce_to_iter(src_files)
    for f in src_files:
        if not f.startswith(src_root):
            raise ValueError('{} not a sub-path of {}'.format(f, src_root))
    dest_files = [
        os.path.join(dest_root, os.path.relpath(f, start=src_root)) \
        for f in src_files
    ]
    for f in dest_files:
        if not overwrite and os.path.exists(f):
            raise OSError('{} exists.'.format(f))
        os.makedirs(os.path.normpath(os.path.dirname(f)), exist_ok=True)
    for src, dest in zip(src_files, dest_files):
        copy_function(src, dest)

def resolve_path(path, root_path="", env=None):
    """Abbreviation to resolve relative paths.

    Args:
        path (:obj:`str`): path to resolve.
        root_path (:obj:`str`, optional): root path to resolve `path` with. If
            not given, resolves relative to `cwd`.

    Returns: Absolute version of `path`, relative to `root_path` if given, 
        otherwise relative to `os.getcwd`.
    """
    def _expandvars(path, env_dict):
        """Expand quoted variables of the form $key and ${key} in path,
        where key is a key in env_dict, similar to os.path.expandvars.

        See `<https://stackoverflow.com/a/30777398>`__; specialize to not skipping
        escaped characters and not changing unrecognized variables.
        """
        return re.sub(
            r'\$(\w+|\{([^}]*)\})', 
            lambda m: env_dict.get(m.group(2) or m.group(1), m.group(0)), 
            path
        )

    if path == '':
        return path # default value set elsewhere
    path = os.path.expanduser(path) # resolve '~' to home dir
    path = os.path.expandvars(path) # expand $VAR or ${VAR} for shell env_vars
    if isinstance(env, dict):
        path = _expandvars(path, env)
    if '$' in path:
        print("Warning: couldn't resolve all env vars in '{}'".format(path))
        return path
    if os.path.isabs(path):
        return path
    if root_path == "":
        root_path = getcwd()
    assert os.path.isabs(root_path)
    return os.path.normpath(os.path.join(root_path, path))

def check_executable(exec_name):
    """Tests if <exec_name> is found on the current $PATH.

    Args:
        exec_name (:py:obj:`str`): Name of the executable to search for.

    Returns: :py:obj:`bool` True/false if executable was found on $PATH.
    """
    return (find_executable(exec_name) is not None)

def poll_command(command, shell=False, env=None):
    """Runs a shell command and prints stdout in real-time.
    
    Optional ability to pass a different environment to the subprocess. See
    documentation for the Python2 `subprocess 
    <https://docs.python.org/2/library/subprocess.html>`_ module.

    Args:
        command: list of command + arguments, or the same as a single string. 
            See `subprocess` syntax. Note this interacts with the `shell` setting.
        shell (:py:obj:`bool`, optional): shell flag, passed to Popen, 
            default `False`.
        env (:py:obj:`dict`, optional): environment variables to set, passed to 
            Popen, default `None`.
    """
    process = subprocess.Popen(
        command, shell=shell, env=env, stdout=subprocess.PIPE)
    while True:
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()
    return rc

class TimeoutAlarm(Exception):
    # dummy exception for signal handling in run_command
    pass

def run_command(command, env=None, cwd=None, timeout=0, dry_run=False):
    """Subprocess wrapper to facilitate running single command without starting
    a shell.

    Note:
        We hope to save some process overhead by not running the command in a
        shell, but this means the command can't use piping, quoting, environment 
        variables, or filename globbing etc.

    See documentation for the Python2 `subprocess 
    <https://docs.python.org/2/library/subprocess.html>`_ module.

    Args:
        command (list of :py:obj:`str`): List of commands to execute
        env (:py:obj:`dict`, optional): environment variables to set, passed to 
            `Popen`, default `None`.
        cwd (:py:obj:`str`, optional): child processes' working directory, passed
            to `Popen`. Default is `None`, which uses parent processes' directory.
        timeout (:py:obj:`int`, optional): Optionally, kill the command's subprocess
            and raise a CalledProcessError if the command doesn't finish in 
            `timeout` seconds.

    Returns:
        :py:obj:`list` of :py:obj:`str` containing output that was written to stdout  
        by each command. Note: this is split on newlines after the fact.

    Raises:
        CalledProcessError: If any commands return with nonzero exit code.
            Stderr for that command is stored in `output` attribute.
    """
    def _timeout_handler(signum, frame):
        raise TimeoutAlarm

    if isinstance(command, six.string_types):
        command = shlex.split(command)
    cmd_str = ' '.join(command)
    if dry_run:
        print('DRY_RUN: call {}'.format(cmd_str))
        return
    proc = None
    pid = None
    retcode = 1
    stderr = ''
    try:
        proc = subprocess.Popen(
            command, shell=False, env=env, cwd=cwd,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, bufsize=1
        )
        pid = proc.pid
        # py3 has timeout built into subprocess; this is a workaround
        signal.signal(signal.SIGALRM, _timeout_handler)
        signal.alarm(int(timeout))
        (stdout, stderr) = proc.communicate()
        signal.alarm(0)  # cancel the alarm
        retcode = proc.returncode
    except TimeoutAlarm:
        if proc:
            proc.kill()
        retcode = errno.ETIME
        stderr += f"\nKilled by timeout ( > {timeout} sec)."
    except Exception as exc:
        if proc:
            proc.kill()
        stderr += f"\nCaught exception {repr(exc)}."
    if retcode != 0:
        print('run_command on {} (pid {}) exit status={}:{}\n'.format(
            cmd_str, pid, retcode, stderr
        ))
        raise subprocess.CalledProcessError(
            returncode=retcode, cmd=cmd_str, output=stderr)
    if '\0' in stdout:
        return stdout.split('\0')
    else:
        return stdout.splitlines()

def run_shell_command(command, env=None, cwd=None, dry_run=False):
    """Subprocess wrapper to facilitate running shell commands.

    See documentation for the Python2 `subprocess 
    <https://docs.python.org/2/library/subprocess.html>`_ module.

    Args:
        commands (list of :py:obj:`str`): List of commands to execute
        env (:py:obj:`dict`, optional): environment variables to set, passed to 
            `Popen`, default `None`.
        cwd (:py:obj:`str`, optional): child processes' working directory, passed
            to `Popen`. Default is `None`, which uses parent processes' directory.

    Returns:
        :py:obj:`list` of :py:obj:`str` containing output that was written to stdout  
        by each command. Note: this is split on newlines after the fact, so if 
        commands give != 1 lines of output this will not map to the list of commands
        given.

    Raises:
        CalledProcessError: If any commands return with nonzero exit code.
            Stderr for that command is stored in `output` attribute.
    """
    # shouldn't lookup on each invocation, but need abs path to bash in order
    # to pass as executable argument. Pass executable argument because we want
    # bash specifically (not default /bin/sh, and we save a bit of overhead by
    # starting bash directly instead of from sh.)
    bash_exec = find_executable('bash')

    if not isinstance(command, six.string_types):
        command = ' '.join(command)
    if dry_run:
        print('DRY_RUN: call {}'.format(command))
        return
    proc = None
    pid = None
    retcode = 1
    stderr = ''
    try:
        proc = subprocess.Popen(
            command,
            shell=True, executable=bash_exec,
            env=env, cwd=cwd,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            universal_newlines=True, bufsize=1
        )
        pid = proc.pid
        (stdout, stderr) = proc.communicate()
        retcode = proc.returncode
    except Exception as exc:
        if proc:
            proc.kill()
        stderr += f"\nCaught exception {repr(exc)}."
    if retcode != 0:
        print('run_shell_command on {} (pid {}) exit status={}:{}\n'.format(
            command, pid, retcode, stderr
        ))
        raise subprocess.CalledProcessError(
            returncode=retcode, cmd=command, output=stderr)
    if '\0' in stdout:
        return stdout.split('\0')
    else:
        return stdout.splitlines()

def is_iterable(obj):
    return isinstance(obj, collections_abc.Iterable) \
        and not isinstance(obj, six.string_types) # py3 strings have __iter__

def coerce_to_iter(obj, coll_type=list):
    assert coll_type in [list, set, tuple] # only supported types for now
    if obj is None:
        return coll_type([])
    elif isinstance(obj, coll_type):
        return obj
    elif is_iterable(obj):
        return coll_type(obj)
    else:
        return coll_type([obj])

def coerce_from_iter(obj):
    if is_iterable(obj):
        if len(obj) == 1:
            return list(obj)[0]
        else:
            return list(obj)
    else:
        return obj

def remove_prefix(s1, s2):
    if s1.startswith(s2):
        s1 = s1[len(s2):]
    return s1

def remove_suffix(s1, s2):
    if s1.endswith(s2):
        s1 = s1[:-len(s2)]
    return s1

def abbreviate_path(path, old_base, new_base=None):
    """Express path as a path relative to old_base, optionally prepending 
    new_base.
    """
    ps = tuple(os.path.abspath(p) for p in (path, old_base))
    str_ = os.path.relpath(ps[0], start=os.path.commonpath(ps))
    if new_base is not None:
        str_ = os.path.join(new_base, str_)
    return str_

def filter_kwargs(kwarg_dict, function):
    """Given a dict of kwargs, return only those kwargs accepted by function.
    """
    named_args = set(six.get_function_code(function).co_varnames)
    # if 'kwargs' in named_args:
    #    return kwarg_dict # presumably can handle anything
    return dict((k, kwarg_dict[k]) for k in named_args \
        if k in kwarg_dict and k not in ['self', 'args', 'kwargs'])

def filter_dataclass(d, dc):
    """Given a dataclass dc (may be the class or an instance of it), and a dict,
    dataclass or dataclass instance d, return a dict of the subset of fields or 
    entries in d that correspond to the fields in dc.
    """
    assert dataclasses.is_dataclass(dc)
    if dataclasses.is_dataclass(d):
        if isinstance(d, type):
            d = d() # d is a class; instantiate with default field values
        d = dataclasses.asdict(d)
    return {f.name: d[f.name] for f in dataclasses.fields(dc) if f.name in d}

def signal_logger(caller_name, signum=None, frame=None):
    """Lookup signal name from number; `<https://stackoverflow.com/a/2549950>`__.
    """
    if signum:
        sig_lookup = {
            k:v for v, k in reversed(sorted(list(signal.__dict__.items()))) \
                if v.startswith('SIG') and not v.startswith('SIG_')
        }
        sig_name = sig_lookup.get(signum, 'UNKNOWN')
        print(f"DEBUG: {caller_name} caught signal {sig_name} ({signum})")
        # if frame:
        #     traceback.print_stack(f=frame)
