# Utility functions
from __future__ import division, print_function
from future.utils import iteritems
import sys, numpy as np, math, os, signal, fnmatch, pandas as pd
# allows for a nicer ordered dict syntax than found in collections.OrderedDict
# e.g. d = odict['key1': 1, 'key2': 2]
from odictliteral import odict
from collections import OrderedDict, Counter
import traceback, warnings
from decimal import Decimal
from skimage.util.shape import view_as_windows as viewW

class defer_signals(object):
    """
    Context manager to defer signal handling until context exits.
    Takes optional list of signals to defer (default: SIGHUP, SIGINT, SIGTERM, SIGTSTP).
    Signals can be identified by number or by name.
    Allows you to wrap instruction sequences that ought to be atomic and ensure
    that they don't get interrupted mid-way.

    Use:
    with defer-signals():
        do_no_interrupt_func() 
    """

    def __init__(self, signal_list=None):
        # Default list of signals to defer
        if signal_list is None:
            signal_list = [signal.SIGHUP, signal.SIGINT, signal.SIGTERM, signal.SIGTSTP]
        # Accept either signal numbers or string identifiers
        self.signal_list = [
            getattr(signal, sig_id) if isinstance(sig_id, basestring) else sig_id
            for sig_id in signal_list
        ]
        self.deferred = []
        self.previous_handlers = {}

    def defer_signal(self, sig_num, stack_frame):
        self.deferred.append(sig_num)

    def __enter__(self):
        # Replace existing handlers with deferred handler
        for sig_num in self.signal_list:
            # signal.signal returns None when no handler has been set in Python,
            # which is the same as the default handler (SIG_DFL) being set
            self.previous_handlers[sig_num] = (
                signal.signal(sig_num, self.defer_signal) or signal.SIG_DFL)
        return self

    def __exit__(self, *args):
        # Restore handlers
        for sig_num, handler in self.previous_handlers.items():
            signal.signal(sig_num, handler)
        # Send deferred signals
        while self.deferred:
            sig_num = self.deferred.pop(0)
            os.kill(os.getpid(), sig_num)

    def __call__(self):
        """
        If there are any deferred signals pending, trigger them now
        This means that instead of this code:
            for item in collection:
                with defer_signals():
                    item.process()
        You can write this:
            with defer_signals() as handle_signals:
                for item in collection:
                    item.process()
                    handle_signals()
        Which has the same effect but avoids having to embed the context
        manager in the loop
        """
        if self.deferred:
            # Reattach the signal handlers and fire signals
            self.__exit__()
            # Put our deferred signal handlers back in place
            self.__enter__()

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def clrd_print(msg, msg_type):
    """
    Prints a message using a message related color.
    
    Parameters
    ----------
    msg : str
        Message
    msg_type : str
        Message type for color highlighting. Choose from:
            'ok', 'warn', 'error'.
    """
    if msg_type == 'ok':
        c = bcolors.OKGREEN
    elif msg_type == 'warn':
        c = bcolors.WARNING
    elif msg_type == 'error':
        c = bcolors.FAIL
    else:
        raise ValueError("Message type '{}' not implemented.\n".format(msg_type))
        
    print(c+msg+bcolors.ENDC)
    
def warn_with_traceback(message, category, filename, lineno, file = None, line = None):
    """
    Prints warnings with stack traceback. Assign to warnings.showwarning to override
    default behavior.
    """
    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

# ------- functions to handle nested dictionaries and dict helper functions -------
#==================================================================================

def expand_keypath(d, keypath, path_mrkr = '/'):
    """
    Expands given keypath for nested dict d by matching and replacing wildcards.

    Parameters
    ----------
    d : dict
        Nested dict.

    keypath : str
        Keypaths within supplied dict "d" containing Unix-style wildcards.

    path_mrkr : str
        Single-character path marker.

    Returns
    -------
    list of list of str
        Expanded keypath of supplied dict. Keys at each nesting level are sorted for a reproducible output.

    Example
    -------
    d = {
        "a": {
            "b1": 1,
            "b2": 2
        },
        "c": {
            "d1": 1,
            "d2": 2
        }
    }

    print(expand_keypath(d, "*/*"))

    output
    ------

    [['a', 'b1'], ['a', 'b2'], ['c', 'd2'], ['c', 'd1']]

    """
    if not keypath:
        return []
    if not isinstance(keypath, list):
        keypath = parse_dict_paths(keypath, path_mrkr)

    matched_keys = sorted(fnmatch.filter(d.keys(), keypath[0]))
    out = []
    for mk in matched_keys:
        child_kp = expand_keypath(d[mk], keypath[1:], path_mrkr = path_mrkr)
        if child_kp:
            for ckp in child_kp:
                out.append([mk])
                out[-1].extend(ckp)
        else:
            out.append([mk])

    return out

def get_dpath(source, keypaths, path_mrkr = '/', ignore_broken_paths = False):
    """
    Obtains values of items in a nested source dict based on a list of key paths, which are placed in a new ordered dict.

    Parameters
    ----------
    source : dict
        Nested source dictionary.
    keypaths : str, iterable of str or dict
        Key paths containing a path marker, e.g. '/' used to look up the value of an item in the nested source dict.
        If dict, dict keys are keypaths and values are str used to rename (e.g longer) keys in the original dict.
    path_mrkr : str
        Path marker for key names in the target dict.
    ignore_missing_paths : bool
        If True, broken or missing paths are not retrieved, otherwise, if False (default) an exception is generated.

    Returns
    -------
    odictliteral.odict
        Ordered dict with source dict values and supplied key paths.

    If ignore_broken_paths == False, and key path is broken, raises KeyError.
    """
    out = odict()
    if isinstance(keypaths, str):
        keypaths = [keypaths]
    elif isinstance(keypaths, dict):
        keypaths = list(keypaths)
    for keypath in keypaths:
        assert keypath # empty string not allowed
        path = parse_dict_paths(keypath, path_mrkr)
        item = source
        item_found = False
        for src_key in path:
            if src_key in item:
                item_found = True
                item = item[src_key]
            elif ignore_broken_paths == False:
                raise KeyError
            else:
                item_found = False
                break
        if item_found:        
            out[keypath] = item

    if isinstance(keypaths, dict):
        renamed_out = odict()
        for k_old, k_new in iteritems(keypaths):
            renamed_out[k_new] = out[k_old]
        out = renamed_out

    return out

def get_dpath_val(source, keypath, path_mrkr = '/'):
    """
    Obtains a single value in a nested source dict given a key path.

    Parameters
    ----------
    source : dict
        Nested source dictionary.
    keypath : str
        Key path containing a path marker, e.g. '/' used to look up the value of the item in the nested source dict.
    path_mrkr : str
        Path marker for key names in the target dict.
    """
    try:
        out = get_dpath(source, [keypath], path_mrkr)[keypath]
    except KeyError:
        raise KeyError("Could not access keypath '{}' in source dict.".format(keypath))
    return out

def set_dpath(source, target, path_mrkr = '/'):
    """
    Sets target nested dictionary keys to the values in a source dict that has keys coded as a path in the nested target dictionary.
    If path keys described in the source dict do not exist in target dict, then target dict keys are created. If path marker is contained
    in the dict key, when specifying the path name, it must be doubled, e.g. 'key1//key2' for setting a key 'key1/key2'.

    Parameters
    ----------
    source : dict
        Single level dictionary with keys describing a path in the target dictionary using a path marker.

    target : dict
        Nested dictionaries.

    path_mrkr : str
        Path marker.

    Returns
    -------
    None
    """
    for src_key, src_val in iteritems(source):
        path = parse_dict_paths(src_key, path_mrkr)
        # visit nested dict items
        item = target
        for trgt_key in path[:-1]:
            # if dict key exists, then assign value
            if trgt_key in item:
                item = item[trgt_key]
            else:
            # if key does not exist, create key and assign value
                item[trgt_key] = {}
                item = item[trgt_key]
        item[path[-1]] = src_val
        
def get_keyval_from_dpath(source, keypath, key):
    """
    Moves from deeper to more superficial keypaths and tries to obtain a given key at each level.

    Parameters
    ----------
    source : dict
        Nested source dictionary
    keypath : str
        Dict keypath as strings separated by '/'. If key contains '/', it must be escaped as '//'.
    key : str
        Key to return value for.

    Returns
    -------
    None if key not found in path, otherwise object of key if found.
    """
    kp = parse_dict_paths(keypath)
    curr_kp = keypath # current keypath
    for level in range(len(kp), 0, -1):
        obj = get_dpath_val(source, curr_kp)
        if isinstance(obj, dict) and key in obj:
            return obj[key]
        curr_kp = '/'.join(kp[:level])

    return None

def parse_dict_paths(path, path_mrkr = '/'):
    """
    Splits a nested dict key path into individual keys. If path_mrkr is contained by a key, it must be doubled, e.g. '//'

    Parameters
    ----------
    path : str
        Nested key path e.g. /key1/key2_1//key2_2 corresponding to a dict of the form {'key1': {'key2_1/key2_2': val}}
    path_mrkr : str
        Path marker, commonly '/'.

    Returns
    -------
    list of str
        Nested dict keys e.g. ['key1', 'key2_1/key2_2']
    """
    out = path.strip(path_mrkr).split(path_mrkr)
    idx = 0
    while idx < len(out):
        if not out[idx]:
            out[idx] = out[idx-1] + path_mrkr + out[idx+1]
            del out[idx-1]
            idx -= 1
            del out[idx+1]
        idx += 1
    return out

def chk_dict_path(d, path, path_mrkr = '/'):
    """
    Check if nested dict key path exists.

    Parameters
    ----------
    d : dict
        Nested dict.
    path : str
        Keys paths separated by path_marker.
    path_marker : str
        Single character path marker of nested keys.

    Returns
    -------
    bool
        True if dict key path exists and False otherwise.
    """
    keys = parse_dict_paths(path, path_mrkr = path_mrkr)
    d_cursor = d
    for k in keys:
        if k in d_cursor:
            d_cursor = d_cursor[k]
        else:
            return False

    return True

def parse_time_interval(trange, delimiter = ','):
    """
    Parses a variable comma delimited time interval range <t start>,<t end> to produce a uniform output.
    If <t start> is omitted, tstart=0, and if <t end> is omitted, tend = None (meaning till end).

    Parameters
    ----------
    trange : str
        Delimited time interval as e.g. "<t start>,<t end>", when using a comma.

    delimiter : str, 1 character
        Interval delimiter.

    Returns
    -------
    tuple
    (tstart,tend)
        tstart : float
            Start time, 0 if start time ommited or provided value.
        tend : float or None
            End time, float or None if end time ommited. 
    """
    if trange:
        split_interval = trange.split(delimiter)
        assert len(split_interval) == 2
        if split_interval[0]:
            tstart = float(split_interval[0])
        else:     
            tstart = 0
        if split_interval[1]:
            tend = float(split_interval[1])
        else:     
            tend = None
    else:
        tstart = 0
        tend = None

    return tstart, tend

def set_default_keys(src, trgt):
    """
    Helper function to set key values in a target dict if keys don't exist using
    values from a source dict.

    Parameters
    ----------
    src, trgt : dict
        Source and target dict.

    Returns
    -------
    None
    """
    for src_key, src_val in iteritems(src):
        if src_key not in trgt:
            trgt[src_key] = src_val

# misc
#========================================================

def round_up_to_odd_int(f):
    """
    Round up to nearest odd integer.

    Parameters
    ----------
    f : float

    Returns
    -------
    int
    """
    return int(np.ceil(f) // 2 * 2 + 1)

def round_to_nearest_odd_int(x):
    """
    Round number to nearest odd int. If x=0, returns 1.

    Parameters
    ----------
    x : float
        Number to round.

    Returns
    -------
    int
    """
    r = x%2
    if r < 1:
       return int(math.floor(x)+1)
    else:
       return int(math.floor(x))

def distribute_subplots(max_nrows, max_ncols, nsubplots, layout = 'flexible'):
    """
    Distributes n subplots into multiple subplot grids of maximum size nrows x ncols.
    Distribution is done first by row, then by column with [0,0] corresponding to
    the upper left corner.

    Parameters
    ----------
    max_nrows, max_ncols : int
        Maximum number of rows and columns in a figure.

    nsubplots : int
        Total number of subplots to distribute.

    layout : str
        Choose 'fixed' to distribute plots on a fixed grid size on each page or
        choose 'flexible' to adjust the grid size depending on the number of plots
        to maximize space filling on each page.
    Returns
    -------
    list of tuple:
        Number of list elements is the number of figures needed to distribute plots. Each tuple is of the form
        ((nrows, ncols), [(plt_idx, row_1, col_1),...(plt_idx, row_N, col_N)])
        where:
        - first tuple element is the current plot index [0, nsubplots).
        - second tuple element is the grid size in # rows and # columns
        - third element is a list with distributed subplot 0-index coordinates
    """
    out = []
    nleft_to_assign = nsubplots
    curr_plt_idx = 0
    while nleft_to_assign:
        if layout == 'flexible':
            ncols = min(max_ncols, nleft_to_assign)
            nrows = min(max_nrows, int(math.ceil(nleft_to_assign/ncols)))
        elif layout == 'fixed':
            ncols = max_ncols
            nrows = max_nrows
        else:
            raise ValueError("Grid layout can be 'flexible' or 'fixed'.")

        # counter for number of subplots assigned for current figure
        plt_ctr = 0
        fig_plots = []
        nplts = min(nleft_to_assign, ncols*nrows)
        while plt_ctr < nplts:
            col_idx = plt_ctr%ncols
            row_idx = int(plt_ctr/ncols)
            fig_plots.append((curr_plt_idx,row_idx,col_idx))
            plt_ctr += 1
            curr_plt_idx += 1

        out.append(((nrows,ncols),fig_plots))
        nleft_to_assign -= nplts

    return out

def get_duplicates(x):
    """
    Returns duplicate items.

    Parameters
    ----------
    x : iterable

    Return
    ------
    list
    """
    return [item for item, count in Counter(x).items() if count > 1]

def del_list_idxs(l, idxs):
    """
    Removes elements from a list in place with given indices, which may repeat.

    Parameters
    ----------
    l : list
        List to remove elements from.
    idxs : list
        List indices to remove.

    Returns
    -------
    None
    """
    for idx in sorted(set(idxs), reverse = True):
        del l[idx]

def resample_rounded_fs(sig, fs, fs_resolution):
    """
    Resamples a signal with a rounded sampling rate using linear interpolation.

    Parameters
    ----------
    sig : 1D numpy.ndarray
        Signal to resample.
    fs : float
        Original sampling rate in [Hz].
    fs_resolution : float
        Sampling rate resolution to use for rounding.

    Returns
    -------
    tuple
        (rounded_fs, resampled_sig)
    """
    rounded_fs = float(Decimal(round(fs/fs_resolution))*Decimal(fs_resolution))
    if sig is None or not len(sig):
        return rounded_fs, np.array([])
    else:
        original_t_vec = np.arange(0, len(sig)/fs, 1./fs)
        adjusted_t_vec = np.arange(0, len(sig)/fs, 1./rounded_fs)
        interp_nsamp = min(len(adjusted_t_vec), len(sig))
        resampled_sig = np.interp(adjusted_t_vec, original_t_vec[:interp_nsamp], sig[:interp_nsamp])
        return rounded_fs, resampled_sig

def group_elements(x, max_n_group):
    """
    Groups an iterable x into groups with a maximum number of elements.
    Parameters
    ----------
    x : iterable
        Input array of elements.
    max_n_group : int
        Maximum number of elements in a group

    Returns
    -------
    list of list
        List of grouped elements
    """
    out = []
    n_full_groups = int(len(x)/max_n_group)
    if n_full_groups:
        out.extend([x[gidx*max_n_group:(gidx+1)*max_n_group] for gidx in range(n_full_groups)])
    # add incomplete group
    incomplete_group = x[n_full_groups*max_n_group:]
    if incomplete_group:
        out.append(incomplete_group)

    return out

flatten_list = lambda l: [item for sublist in l for item in sublist]

def df_to_xls(dfs, fpath, float_format = "%.2f"):
    """
    Prints multiple dataframes to an excel file.

    Parameters
    ----------
    dfs : list of iterable of 2 elements of str and pandas.DataFrame types
        List of sheet name and pandas dataframe pairs.

    fpath : str
        File path to save data.
    """
    writer = pd.ExcelWriter(fpath, engine = 'xlsxwriter')
    for sheetname, df in dfs:
        df.to_excel(writer, sheet_name = sheetname, float_format = float_format)  # send df to writer
        worksheet = writer.sheets[sheetname]  # pull worksheet object
        for idx, col in enumerate(df):  # loop through all columns
            series = df[col]
            max_len = max((
                series.astype(str).map(len).max(),  # len of largest item
                len(str(series.name))  # len of column name/header
                )) + 1  # adding a little extra space
            worksheet.set_column(idx, idx, max_len)  # set column width
    writer.save()

def find_nearest(a, val):
    """
    Finds nearest value and index thereof in an array from a reference value.

    Parameters
    ----------
    a :  iterable
        Array
    val : int, float
        Reference value

    Returns
    -------
    tuple
        (idx, val)
        idx - array index
        val - value
    """
    a = np.asarray(a)
    idx = (np.abs(a-val)).argmin()
    return idx, a[idx]

def query_yes_no(question, default = None):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

def make_subfolders(parent_path, child_path):
    """
    Combines a child and parent folder path and creates the full path if needed.

    Parameters
    ----------
    parent_path : str
        Parent folder path.
    child_path : str
        Child folder path.

    Returns
    -------
    str
        Merged folder path.
    """
    full_path = os.path.join(parent_path, child_path)
    if not os.path.exists(full_path):
        os.makedirs(full_path)
    return full_path

def remove_nans(x):
    """
    Removes nans from an array.

    Parameters
    ----------
    x : numpy.ndarray
        Input array

    Returns
    -------
    numpy.ndarray
        New array with nan elements removed
    """
    return x[~np.isnan(x)]

def shift_rows(m, s):
    """
    Applies varying shifts to each row of a matrix.
    Parameters
    ----------
    m : 2D numpy.ndarray
        Input matrix.
    s : 1D numpy.array
        Array of shifts to apply to each row.
    Returns
    -------
    2D numpy.ndarray
        Matrix with shifted rows and nan padding.
    """
    p = np.full((m.shape[0], m.shape[1]-1), np.nan)
    m_ext = np.concatenate((p,m,p), axis = 1)

    # get sliding windows
    n = m.shape[1]
    return viewW(m_ext,(1,n))[np.arange(len(s)), -s + (n-1),0]