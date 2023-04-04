from __future__ import division
from future.utils import iteritems
import numpy as np, math, numbers
from collections import OrderedDict
from itertools import cycle as iter_cycle
import matplotlib as mpl
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib import pyplot as plt
import matplotlib.offsetbox
from matplotlib.lines import Line2D
import process as proc, util
import os, copy
from scipy import stats
import mne
import statsmodels.tsa.stattools as tsa_stat
import seaborn as sb

# color palettes
# color name, pyplot named color/0-1 rgb values/hex codes
cpal = {
    # seaborn colorblind palette
    "sb-cb": 
        OrderedDict([
            ("darkblue" , (0.00392156862745098, 0.45098039215686275, 0.6980392156862745)),
            ("darkorange", (0.8705882352941177, 0.5607843137254902, 0.0196078431372549)),
            ("seagreen", (0.00784313725490196, 0.6196078431372549, 0.45098039215686275)),
            ("chocolate", (0.8352941176470589, 0.3686274509803922, 0.0)),
            ("darkpurple", (0.8, 0.47058823529411764, 0.7372549019607844)),
            ("wood", (0.792156862745098, 0.5686274509803921, 0.3803921568627451)),
            ("lightpurple", (0.984313725490196, 0.6862745098039216, 0.8941176470588236)),
            ("grey", (0.5803921568627451, 0.5803921568627451, 0.5803921568627451)),
            ("yellow", (0.9254901960784314, 0.8823529411764706, 0.2)),
            ("lightblue", (0.33725490196078434, 0.7058823529411765, 0.9137254901960784))
        ])
}

# default plotting parameters
letter_portrait = (8.5, 11) # width, height, portrait
letter_landscape = (11, 8.5) # width, height, landscape
default_figsize = (20, 5)
default_print_fig_width = 8
default_label_fontsize = 16
default_title_fontsize = 16
default_plt_markersize = 2 # markersize in pt
default_legend_fontsize = 12
default_tick_label_fontsize = 14
default_nyticks = 5
default_linewidth = 1.5
default_tick_label_pad = 5
default_fig_dpi = 300

# Global plotting settings to have uniform figures
# =====================================================================================
# matplotlib RC parameters
mpl_rc_context = {
    "paper": {
        "font.size": 8,
        "axes.titlesize": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "pdf.fonttype": 42
    }
}
# width and height of a single 89 mm column figure with golden aspect ratio in inches
golden_ratio = 1.618
single_column_89mm_width = 89/25.4 
double_column_183mm_width = 183/25.4
single_column_golden_size_89mm_wide = (3.50394, 3.50394/golden_ratio)
single_column_golden_size_89mm_wide_portrait = (3.50394, 3.50394*golden_ratio)
double_column_golden_size_183mm_wide = (7.20472, 7.20472/golden_ratio)
double_column_golden_size_183mm_wide_portrait = (7.20472, 7.20472*golden_ratio)

class ArbitraryColors:
    def __init__(self, palette = None):
        """No guarantee these colors are aesthetically pleasing or
        do not repeat
        """
        if palette is None:
            self.palette = sb.color_palette()
        else:
            self.palette = palette
        self.n_colors = len(self.palette)
    
    def __getitem__(self, key):
        return self.palette[hash(key)%self.n_colors]
    
    def __call__(self, key):
        return self[key]

RANDOM_COLORS = ArbitraryColors(sb.color_palette('bright'))       

class CircularColors:
    """
    Class used to return colors from a circular named list.
    """
    def __init__(self, palette = 'sb-cb', exclude = [], alpha = None):
        """
        Initializes circular colors object.

        Parameters
        ----------
        palette : str
            Name of color palette. Choose from 'sb-cb'.
        exclude : str
            Exclude a certain color name from the palette.
        alpha : None, float
            Alpha trasparency channel. If specified should be [0,1]. Returned colors will be 4element tuples.
        """
        self.alpha  = alpha
        self._cdict = cpal[palette]
        for ex in exclude:
            del self._cdict[ex]

        self._colcycle = iter_cycle(self._cdict.values())

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._colcycle)

    # for python 2.7 compatibility
    def next(self):
        if self.alpha is None:
            return self.__next__()
        else:
            return self.__next__()+(self.alpha,)

    def __getitem__(self, key):
        if self.alpha is None:
            return self._cdict[key]
        else:
            return self._cdict[key]+(self.alpha,)

    def showpal(self):
        """
        Show current color palette and print color names.
        """
        print(self._cdict.keys())
        sb.palplot(self._cdict.values())

class AnchoredHScaleBar(matplotlib.offsetbox.AnchoredOffsetbox):
    """ size: length of bar in data units
        extent : height of bar ends in axes units """
    def __init__(self, size = 1, extent = 0.03, label = "", loc = 2, ax = None,
                 pad = 0.4, borderpad = 0.5, ppad = 0, sep = 2, prop = None, 
                 frameon = True, linekw = {}, label_fontsize = 8, **kwargs):
        if not ax:
            ax = plt.gca()
        trans = ax.get_xaxis_transform()
        size_bar = matplotlib.offsetbox.AuxTransformBox(trans)
        line = Line2D([0,size],[0,0], **linekw)
        vline1 = Line2D([0,0],[-extent/2.,extent/2.], **linekw)
        vline2 = Line2D([size,size],[-extent/2.,extent/2.], **linekw)
        size_bar.add_artist(line)
        size_bar.add_artist(vline1)
        size_bar.add_artist(vline2)
        txt = matplotlib.offsetbox.TextArea(label, minimumdescent = False, textprops = dict(fontsize = label_fontsize))
        self.vpac = matplotlib.offsetbox.VPacker(children = [size_bar,txt],  
                                 align = "center", pad = ppad, sep = sep) 
        matplotlib.offsetbox.AnchoredOffsetbox.__init__(self, loc, pad = pad, 
                 borderpad = borderpad, child = self.vpac, prop = prop, frameon = frameon,
                 **kwargs)

def plot_vspaced_waveforms(sig, pset, fig_title = '', sig_f_name = ''):
    """
    Plots grouped vertically stacked waveforms.

    Parameters
    ----------
    sig : dict
        Aligned signals dict.

    pset : dict
        Plotting parameter settings, with keys:
            "max_groups_per_page" : optional, int, default 1
                Maximum number of signal groups to place on a page
            "min_groups_per_page" : optional, int, default 1
                Minimum number of plot groups to create a plot. Useful when goal is to focus on only more than one group and see relationship between signals.
            "fs": str
                Path to sampling rate info within the aligned signals file.
            "plot_groups" : list
                Description of plot groups.
            "figsize": iterable of 2 elements
                Page width and height in inches.
            "xgrid" : bool
                Adds a faint vertical grid at major and minor x-axis locators.
            "tslice" : optional, list of 2
                Beginning and end time interval to use for cropping signals in [s]. Use None for end time if plotting till end of signal. Default is no crop, i.e. [0, None]

    fig_title : str
        Figure title.

    sig_f_name : str
        Aligned signals file name.

    Returns
    -------
    list of matplotlib.Figure
        Generates multiple figures with requested signals.
    """
    if 'fs' in pset:
        fs = util.get_dpath_val(sig, pset['fs'])
    else:
        fs = None

    figs = []
    # expand signal groups
    expanded_groups = expand_imaging_plot_groups(sig, pset['plot_groups'])
    
    # crop to time interval
    if 'tslice' in pset:
        crop_interval = pset['tslice']
    else:
        crop_interval = [0, None]

    # maximum number of groups per page
    if 'max_groups_per_page' in pset:
        max_groups_per_page = pset['max_groups_per_page']
    else:
        max_groups_per_page = 1
    # minimum number of groups per page to include a page
    if 'min_groups_per_page' in pset:
        min_groups_per_page = pset['min_groups_per_page']
    else:
        min_groups_per_page = 1

    # return no figures if the minimum number of groups condition is not satisfied
    if len(expanded_groups) < min_groups_per_page: 
        return []

    # chunk signal groups up to maximum number of groups per page
    chunked_groups = util.group_elements(expanded_groups, max_groups_per_page)
    for cg_idx, cg in enumerate(chunked_groups):
        # iterate over signals in group
        sig_idx = 0 # index each signal within the group
        plot_data = {}
        groups = []
        for g in cg:
            groups.append(list(range(sig_idx,sig_idx+len(g))))
            for s in g:
                # prepare signals to plot
                plot_data[sig_idx] = copy.deepcopy(s)
                
                # replace keypaths with signal and highlight
                if 'sig' in plot_data[sig_idx] and plot_data[sig_idx]['sig'] is not None:
                    if util.chk_dict_path(sig, plot_data[sig_idx]['sig']):
                        plot_data[sig_idx]['sig'] = util.get_dpath_val(sig, plot_data[sig_idx]['sig'])
                        if "invert_around_unit_baseline" in plot_data[sig_idx] and plot_data[sig_idx]["invert_around_unit_baseline"]:
                            plot_data[sig_idx]['sig'] = 1-plot_data[sig_idx]['sig']
                    else:
                        if sig_f_name:
                            util.clrd_print("Warning: Found broken signal path '{}' in file '{}'.".format(plot_data[sig_idx]['sig'], sig_f_name), 'warn')
                        else:
                            util.clrd_print("Warning: Found broken signal path '{}'.".format(plot_data[sig_idx]['sig']), 'warn')
                        del plot_data[sig_idx]['sig']
                
                # replace keypaths with signal and highlight
                if 'sig2' in plot_data[sig_idx] and plot_data[sig_idx]['sig2'] is not None:
                    if util.chk_dict_path(sig, plot_data[sig_idx]['sig2']):
                        plot_data[sig_idx]['sig2'] = util.get_dpath_val(sig, plot_data[sig_idx]['sig2'])
                        if "invert_around_unit_baseline" in plot_data[sig_idx] and plot_data[sig_idx]["invert_around_unit_baseline"]:
                            plot_data[sig_idx]['sig2'] = 1-plot_data[sig_idx]['sig2']
                    else:
                        if sig_f_name:
                            util.clrd_print("Warning: Found broken signal path '{}' in file '{}'.".format(plot_data[sig_idx]['sig2'], sig_f_name), 'warn')
                        else:
                            util.clrd_print("Warning: Found broken signal path '{}'.".format(plot_data[sig_idx]['sig2']), 'warn')
                        del plot_data[sig_idx]['sig2']

                # used for plotting a horizontal line marking a a level (useful for showing a threshold)
                if 'level_marker' in plot_data[sig_idx]:
                    if plot_data[sig_idx]['level_marker'] is None:
                        del plot_data[sig_idx]['level_marker']
                    else:
                        if util.chk_dict_path(sig, plot_data[sig_idx]['level_marker']):
                            plot_data[sig_idx]['level_marker'] = util.get_dpath_val(sig, plot_data[sig_idx]['level_marker'])
                        else:
                            if sig_f_name:
                                util.clrd_print("Warning: Found broken level path '{}' in file '{}'.".format(plot_data[sig_idx]['level_marker'], sig_f_name), 'warn')
                            else:
                                util.clrd_print("Warning: Found broken level path '{}'.".format(plot_data[sig_idx]['level_marker']), 'warn')
                            del plot_data[sig_idx]['level_marker']

                        
                if 'highlights' in plot_data[sig_idx]:
                    for h in plot_data[sig_idx]['highlights']:
                        if 'epochs' in h:
                            if util.chk_dict_path(sig, h['epochs']) and plot_data[sig_idx]['sig'] is not None:
                                epochs = util.get_dpath_val(sig, h['epochs'])
                                if isinstance(epochs, np.ndarray) and epochs.dtype == "bool":
                                    h['epochs'] = epochs
                                else:
                                    h['epochs'] = proc.intervals_to_bool(epochs, fs, len(plot_data[sig_idx]['sig']))
                            else:
                                if sig_f_name:
                                    util.clrd_print("Warning: Found broken highlight epochs path '{}' in file '{}'.".format(h['epochs'], sig_f_name), 'warn')
                                else:
                                    util.clrd_print("Warning: Found broken highlight epochs path '{}'.".format(h['epochs']), 'warn')
                                del h['epochs']
                        if 'events' in h and h['events'] is not None:
                            if util.chk_dict_path(sig, h['events']):
                                h['events'] = util.get_dpath_val(sig, h['events'])
                            else:
                                if sig_f_name:
                                    util.clrd_print("Warning: Found broken event marker path '{}' in file '{}'.".format(h['events'], sig_f_name), 'warn')
                                else:
                                    util.clrd_print("Warning: Found broken event marker path '{}'.".format(h['events']), 'warn')
                                del h['events']
                sig_idx += 1
        # generate figure
        fig, _ = vspaced_waveforms(plot_data, fs = fs, groups = groups, fig = None, bg_to_wg_margin = 7.5, wg_to_ax_margin = 0.1,
            figsize = pset['figsize'] if 'figsize' in pset else letter_landscape,
            xgrid = pset['xgrid'] if 'xgrid' in pset else False,
            crop = crop_interval)

        fig.suptitle("%s (page %d/%d)"%(fig_title, cg_idx+1,len(chunked_groups)), fontsize = default_title_fontsize)
        figs.append(fig)

    return figs

def expand_imaging_plot_groups(sig, plot_groups, path_mrkr = '/'):
    """
    Expands plot groups given aligned signals dict.

    Parameters
    ----------
    sig : dict
        Aligned signals

    plot_groups : list
        List of plot groups of list of dict:
            Signal plotting specification

        Example
        -------

        [
            // first plot group, if containing imaging channels, will be expanded by iteration over ROI group and ROI target and filtering
            [
                {
                    // * can be used for ROI group and ROI target, allows all 
                    // specific ROI groups can be named singly or separated by ",".
                    "sig": "imaging/ROIs/spine,spine-shaft/*/Ch1_db_lp_filt/sig",
                    "sig2": "imaging/ROIs/spine,spine-shaft/*/Ch2_db_lp_filt/sig", // plotted on second y-axis
                    "color": "b",
                    "y_axis_label": "{ROItarget}",
                    "highlights": [
                        {
                            "epochs": "imaging/ROIs/spine,spine-shaft/*/Ch1_db_lp_filt/analysis/events/depol",
                            "color": "y"
                        }
                    ],
                    "level_marker": "imaging/ROIs/spine,spine-shaft/*/Ch1_db_lp_filt/analysis/level"
                },
                {   "sig": "imaging/ROIs/spine,spine-shaft/*/Ch2_db_lp_filt/sig",
                    "highlights": []
                }
                // more signals can be added to the group
            ]
            // more groups can be added here
        ]

        Notes
        -----
        "y_axis_label":
            If using "{ROItarget}", and signal originates from an imaging channel, it will be replaced after expansion.

    Returns
    -------
    list
        List of plot groups with expanded imaging keypaths.
    """
    def needs_expansion(keypath, path_mrkr):
        """
        Checks if keypath expansion is needed.

        Returns
        -------
        bool
            True if expansion is needed.
        """
        if keypath.startswith('imaging/ROIs/'):
            parsed_keypath = util.parse_dict_paths(keypath, path_mrkr = path_mrkr)
            if len(parsed_keypath) >= 4 and ('*' in [parsed_keypath[2], parsed_keypath[3]] or ',' in [parsed_keypath[2], parsed_keypath[3]]):
                return True
            else:
                return False
        else:
            return False

    def adjust_keypath(keypath, roi_group_name, roi_target_name, path_mrkr):
        """
        Adjusts imaging keypath by filtering given ROI group and ROI target names.
        """
        expanded = False
        if keypath.startswith('imaging/ROIs/'):
            parsed_keypath = util.parse_dict_paths(keypath, path_mrkr = path_mrkr)
            filtered_ROI_group_name = None
            filtered_ROI_target_name = None
            if needs_expansion(keypath, path_mrkr):
                # filter ROI group name
                if parsed_keypath[2] == '*' or roi_group_name in parsed_keypath[2].split(','):
                    # use all available ROI group names
                    filtered_ROI_group_name = roi_group_name
                
                # filter ROI target name
                if parsed_keypath[3] == '*' or roi_target_name in parsed_keypath[3].split(','):
                    # use all available ROI group names
                    filtered_ROI_target_name = roi_target_name

                if filtered_ROI_group_name is not None and filtered_ROI_target_name is not None:
                    out = path_mrkr.join(['imaging', 'ROIs', filtered_ROI_group_name, filtered_ROI_target_name]+parsed_keypath[4:])
                    expanded = True
                else:
                    out = None
            elif len(parsed_keypath) >= 4:
                filtered_ROI_group_name = parsed_keypath[2]
                filtered_ROI_target_name = parsed_keypath[3]
                out = keypath

        else:
            out = keypath

        return out, expanded
     
    # expanded plot groups
    expanded_pg = []
    for pg in plot_groups:
        # check if plot group needs to be expanded
        expansion_needed = False
        # ps_set - plot signal settings
        for ps_set in pg:
            # check if any signal within the group is an imaging channel
            if 'sig' in ps_set and needs_expansion(ps_set['sig'], path_mrkr):
                expansion_needed = True
            if 'sig2' in ps_set and needs_expansion(ps_set['sig2'], path_mrkr):
                expansion_needed = True
            # check if any signal highlight within the group is from an imaging channel
            if 'highlights' in ps_set and (any([needs_expansion(d['epochs'],path_mrkr) for d in ps_set['highlights'] if 'epochs' in d]) or
                any([needs_expansion(d['events'],path_mrkr) for d in ps_set['highlights'] if 'events' in d])):
                expansion_needed = True
            if 'level_marker' in ps_set and needs_expansion(ps_set['level_marker'], path_mrkr):
                expansion_needed = True

        # expand imaging keypaths
        if expansion_needed:
            for roi_group_name, roi_group in sorted(iteritems(sig['imaging']['ROIs'])):
                for roi_target_name, roi_target in sorted(iteritems(roi_group)):
                    # ps_set - plot signal settings
                    adjusted_pset = []
                    expanded = False # this is to ensure that only expanded groups are collected
                    for ps_set in pg:
                        ps_set_copy = copy.deepcopy(ps_set)

                        # expand keypaths
                        adjusted_chan_keypath = None
                        if 'sig' in ps_set:
                            adjusted_chan_keypath, _expanded = adjust_keypath(ps_set['sig'], roi_group_name, roi_target_name, path_mrkr)
                            expanded = expanded or _expanded

                        adjusted_chan_keypath2 = None
                        if 'sig2' in ps_set:
                            adjusted_chan_keypath2, _expanded = adjust_keypath(ps_set['sig2'], roi_group_name, roi_target_name, path_mrkr)
                            expanded = expanded or _expanded                            

                        # expand level marker
                        adjusted_level_marker_keypath = None
                        if 'level_marker' in ps_set:
                            adjusted_level_marker_keypath, _expanded = adjust_keypath(ps_set['level_marker'], roi_group_name, roi_target_name, path_mrkr)
                            expanded = expanded or _expanded

                        adjusted_highlight_keypaths = []
                        adjusted_event_marker_keypaths = []
                        if 'highlights' in ps_set:
                            for psh in ps_set['highlights']:
                                if 'epochs' in psh:
                                    adjusted_highlight_keypath, _expanded = adjust_keypath(psh['epochs'], roi_group_name, roi_target_name, path_mrkr)
                                    expanded = expanded or _expanded
                                    adjusted_highlight_keypaths.append(adjusted_highlight_keypath)
                                else:
                                    adjusted_highlight_keypaths.append(None)
                                if 'events' in psh:
                                    adjusted_event_marker_keypath, _expanded = adjust_keypath(psh['events'], roi_group_name, roi_target_name, path_mrkr)
                                    expanded = expanded or _expanded
                                    adjusted_event_marker_keypaths.append(adjusted_event_marker_keypath)
                                else:
                                    adjusted_event_marker_keypaths.append(None)

                        # keep only if either signal or at least one highlight channel are present
                        if adjusted_chan_keypath is not None or any([x is not None for x in adjusted_highlight_keypaths]):
                            ps_set_copy['sig'] = adjusted_chan_keypath
                            ps_set_copy['level_marker'] = adjusted_level_marker_keypath
                            if 'highlights' in ps_set_copy:
                                for psh_idx, psh in enumerate(ps_set_copy['highlights']):
                                    if adjusted_highlight_keypaths[psh_idx] is not None:
                                        psh['epochs'] = adjusted_highlight_keypaths[psh_idx]
                                for psh_idx, psh in enumerate(ps_set_copy['highlights']):
                                    if adjusted_event_marker_keypaths[psh_idx] is not None:    
                                        psh['events'] = adjusted_event_marker_keypaths[psh_idx]
                            if adjusted_chan_keypath2 is not None:
                                ps_set_copy['sig2'] = adjusted_chan_keypath2
                            adjusted_pset.append(ps_set_copy)

                    if adjusted_pset and expanded:
                        expanded_pg.append(adjusted_pset)
        else:
            expanded_pg.append(pg)

    # replace placeholder names in "y_axis_label"
    for pg in expanded_pg:
        for plot_item in pg:
            if 'y_axis_label' in plot_item and plot_item['y_axis_label'] and plot_item['sig'] is not None:
                if plot_item['sig'].startswith('imaging'+path_mrkr+'ROIs'):
                    split_sig_path = plot_item['sig'].split(path_mrkr)
                    if len(split_sig_path) >= 4:
                        plot_item['y_axis_label'] = plot_item['y_axis_label'].replace(r'{ROIgroup}', split_sig_path[2])
                        plot_item['y_axis_label'] = plot_item['y_axis_label'].replace(r'{ROItarget}', split_sig_path[3])
                    
    return expanded_pg

def vspaced_waveforms(plot_data, fs = None, groups = [], crop = (0, None), fig = None, bg_to_wg_margin = 7.5, wg_to_ax_margin = 0.1,
    fig_margins = (0.1, 0.025, 0.1, 0.1), xlabel = '', figsize = (21,3), y_axis_label_offset = -0.04, xgrid = True):
    """
    Plot multiple grouped waveforms and associated highlights that share the same x-axis.

    Parameters
    ----------
    plot_data : dict of dict
        Signals and associated highlights, dict with keys:
            'sig' : 1D numpy.ndarray
                First y-axis signal.
            'sig2' : optional, 1D numpy.ndarray  
                Second y-axis signal.
            'color': str
                Signal color, use pyplot compatible names.
            'y_axis_label': str
                Y-axis label, can contain placeholder names: "{ROItarget}", "{ROIgroup}" which will be replaced
                by imaging ROI target and group if imaging channel is accessed.
            'highlights': list of dict
                Signal highlights, dict with keys:
                    'epochs': 
                        1) iterable of 2 element iterable
                            Intervals to highlight, e.g. list of 2 element tuples (start, end).
                            If sampling frequency is supplied, itervals are measured in [s], otherwise, array index.
                        2) iterable of bool type, e.g. 1D numpy bool array
                            Highlighting is applied to True array elements.
                    'color': optional, str, defaults to 'b' blue
                        Highlight color, pyplot specifier.
                    'alpha': optional, float, default 0.5
                        Highlight transparency.
                Note: highlights are plotted first and waveforms plotted over them.
            'level_marker': float
                Plots a faint transparent horizontal grey line marking a signal level.

    fs : None or float
        If specified, sampling frequency in [Hz], otherwise array index is used.
        
    crop : tuple of int or float
        If sampling frequency provided, start and end times in [s] to plot as (start_t, end_t), otherwise tuple of array indices.
        If waveform has no end boundary, pass None to the right crop tuple.
        
    groups : list of iterable of plot_data keys
        Name of signals to plot and way in which to group them. If empty, all are plotted evenly spaced.
    
    bg_to_wg_margin : float
        Ratio of margin between groups to within groups. Should be >=1.
        
    wg_to_ax_margin : float
        Within group margin relative to axes height.
        
    fig_margins : tuple 
        Figure margins as (left, right, top, bottom).

    xgrid : bool
        Add a faint vertical grid to the plot at major and minor locators.

    Returns
    -------
    (fig, ax) : tuple
        Figure object and axes as dict using same keys as in groups (or all keys from sig if no groups specified).
        
    """
    assert bg_to_wg_margin >= 1
    def count_n_axes(x):
        """
        Count number of axes from the groups specifier.
        """
        if isinstance(x, basestring):
            return 1
        else:
            try:
                return len(x)
            except:
                return 1
    # construct default groups
    if not groups:
        groups = (tuple(plot_data.keys()),)
        no_groups = True
    else:
        no_groups = False
        
    # default figure
    if fig is None:
        fig = plt.figure(figsize = figsize, dpi = default_fig_dpi)
            
    # total number of axes to use for plotting
    n_ax = sum([count_n_axes(x) for x in groups])
    # height of axes relative to figure height (all axes have same height)
    ax_height = (1-fig_margins[2]-fig_margins[3])/(bg_to_wg_margin*wg_to_ax_margin*(len(groups)-1)+n_ax+wg_to_ax_margin*(n_ax-len(groups)))
    # margin within groups relative to figure height
    mwg = wg_to_ax_margin*ax_height
    # margin between groups relative to figure height
    mbg = bg_to_wg_margin*mwg
    # group bottom in normalized figure coordinates
    g_bottom = [mbg*(len(groups)-i-1)+fig_margins[3]+sum([count_n_axes(groups[j])*(ax_height+mwg)-mwg for j in range(i+1,len(groups))]) for i in range(len(groups))]
    # axes bottom in normalized figure coordinates
    ax_bottom = [[g_bottom[g_idx]+(count_n_axes(groups[g_idx])-ax_idx-1)*(ax_height+mwg) for ax_idx in range(count_n_axes(groups[g_idx]))] for g_idx in range(len(groups))]
        
    # plot signals
    ax = {}
    for g_idx in range(len(groups)):
        for ax_idx in range(len(ax_bottom[g_idx])):
            if ax_idx:
                ax[groups[g_idx][ax_idx]] = fig.add_axes([fig_margins[0], ax_bottom[g_idx][ax_idx], 1-fig_margins[0]-fig_margins[1], ax_height], sharex = ax[groups[0][0]])
            else:
                ax[groups[g_idx][ax_idx]] = fig.add_axes([fig_margins[0], ax_bottom[g_idx][ax_idx], 1-fig_margins[0]-fig_margins[1], ax_height])

            # plot waveforms and their highlights
            if 'sig' in plot_data[groups[g_idx][ax_idx]]:
                waveform(
                    sig = plot_data[groups[g_idx][ax_idx]]['sig'],
                    sig2 = plot_data[groups[g_idx][ax_idx]]['sig2'] if 'sig2' in plot_data[groups[g_idx][ax_idx]] else None,
                    fs = fs,
                    highlights = plot_data[groups[g_idx][ax_idx]]['highlights'] if 'highlights' in plot_data[groups[g_idx][ax_idx]] else [],
                    level_marker = plot_data[groups[g_idx][ax_idx]]['level_marker'] if 'level_marker' in plot_data[groups[g_idx][ax_idx]] else None,
                    crop = crop, title = '',
                    ylabel = plot_data[groups[g_idx][ax_idx]]['y_axis_label'] if 'y_axis_label' in plot_data[groups[g_idx][ax_idx]] else '',
                    xlabel = xlabel,
                    color = plot_data[groups[g_idx][ax_idx]]['color'] if 'color' in plot_data[groups[g_idx][ax_idx]] else 'b',
                    alpha = 1, linewidth = 1, ax = ax[groups[g_idx][ax_idx]], xgrid = xgrid)
            
            # set tight plotting on y-axis
            #ax[groups[g_idx][ax_idx]].autoscale(enable = True, axis = 'y', tight = True)
            # hide top, bottom and right spines
            ax[groups[g_idx][ax_idx]].spines["top"].set_visible(False)
            ax[groups[g_idx][ax_idx]].spines["right"].set_visible(False)
            if not (g_idx == len(groups)-1 and ax_idx == len(ax_bottom[g_idx])-1):
                ax[groups[g_idx][ax_idx]].tick_params(
                    axis = 'x',          # changes apply to the x-axis
                    which = 'both',      # both major and minor ticks are affected
                    bottom = False,      # ticks along the bottom edge are off
                    top = False,         # ticks along the top edge are off
                    labelbottom = False) # labels along the bottom edge are off
                ax[groups[g_idx][ax_idx]].xaxis.label.set_visible(False)
                ax[groups[g_idx][ax_idx]].spines["bottom"].set_visible(False)
            else:
                # adjust bottom axis position to detach it from y axis
                ax[groups[g_idx][ax_idx]].spines['bottom'].set_position(('outward', 12))
            # adjust fontsize as fraction of axes height
            ax[groups[g_idx][ax_idx]].tick_params(axis = 'both', which = 'major', labelsize = 10)
            
            # align y-axis labels between stacked sublots
            ax[groups[g_idx][ax_idx]].yaxis.set_label_coords(y_axis_label_offset, 0.5)

            # adjust plot limits
            if 'axis_ylim' in plot_data[groups[g_idx][ax_idx]]:
                ax[groups[g_idx][ax_idx]].set_ylim(plot_data[groups[g_idx][ax_idx]]['axis_ylim'])

    return (fig, ax)

def highlight_waveform(epochs, fs, ax, waveform, box_level = None, box_height = None, crop = (0, None), color = 'y', alpha = 0.5):
    """
    Highlights portions of a waveform.

    Parameters
    ----------
    epochs : 1D numpy.ndarray of bool or iterable of 2 element iterable
        Box highlighting epochs as either:
        - Binary signal, with False no highlight and True otherwise. Binary signal mut be timed w.r.t. the waveform to be highlighted.
        - Iterable of array indices or time intervals

    fs : float or None
        If provided, waveform sampling frequency in [Hz] to match time axes between highlighting boxes and target waveform.

    ax : matplotlib.axes._subplots.AxesSubplot
        Axes object of an existing plot containing the target waveform to be highlighted.

    waveform : None or 1D numpy.ndarray
        If provided, waveform to highlight. If highlight is automatic, min and max of waveform will determine highlighting rectangle placement.

    box_level : None or float
        Center vertical extent of the box. If None and waveform is provided, box center is midway between the signal max and min, otherwise it is set to 0.5.

    box_height : None or float
        Height of the highligting box matching the units and scale of the target signal. If None and waveform is provided, it is set to include waveform min and max,
        otherwise it is set to 1.

    crop : tuple
        Start and end time in [s] or array indices used to crop the highlighting signal and waveform to be considered.
    
    color : str
        Box color.

    alpha : float
        Filled box alpha transparency.
    """
    

    # adjust highlighting box display based on waveform extent
    if waveform is not None:
        cropped_waveform = proc.tslice(sig = waveform, fs = fs, t_slice = crop)
        if len(cropped_waveform) and np.any(np.isfinite(cropped_waveform)):
            if box_level is None:
                box_level = (np.nanmax(cropped_waveform)+np.nanmin(cropped_waveform))/2.
            if box_height is None:
                box_height = np.nanmax(cropped_waveform)-np.nanmin(cropped_waveform)
        else:
            # default highlight shape
            if box_level is None:
                box_level = 0.5
            if box_height is None:
                box_height = 1
    else:
        # default highlight shape
        if box_level is None:
            box_level = 0.5
        if box_height is None:
            box_height = 1
             
    # epochs must be 1D numpy.ndarray
    if isinstance(epochs, np.ndarray):
        # get indices for False->True transitions
        t = (epochs[:-1] < epochs[1:]).nonzero()[0]
        if epochs[0]:
            t = np.insert(t, 0, 0)
        # draw boxes
        for t_idx in t:
            # get box width
            w = 0
            while t_idx+w+1 < len(epochs) and epochs[t_idx+w+1]:
                w += 1
            # add highlighting box
            if fs is None:
                x = t_idx+1
                width = w
            else:
                x = (t_idx+1)/fs
                width = w/fs
            
            if x >= crop[0] and (crop[1] is None or x+width <= crop[1]):
                box = mpl.patches.Rectangle((x,box_level-box_height/2.), width, height = box_height, fill = True, color = color, alpha = alpha, linewidth = 0)
                ax.add_patch(box)
    else:
        # intervals
        for ep in epochs:
            if ep[0] >= crop[0] and (crop[1] is None or ep[1] <= crop[1]):
                box = mpl.patches.Rectangle((ep[0],box_level-box_height/2.), ep[1]-ep[0], height = box_height, fill = True, color = color, alpha = alpha, linewidth = 0)
                ax.add_patch(box)

def waveform(sig, sig2 = None, fs = None, highlights = [], level_marker = None, crop = (0, None), title = '', ylabel = '', ylabel2 = '',
    xlabel = '', color = 'b', color2 = 'b', marker = '', marker2 = '', linestyle = '-', linestyle2 = '-', alpha = 1, alpha2 = 1,
    linewidth = 1, linewidth2 = 1, figsize = (21,3), ax = None, use_second_yaxis = False, legend_label = '', legend_label2 = '', xgrid = True):
    """
    Plots a waveform.
    
    Parameters
    ----------
    sig : np.array
        Signal to plot using first (left) y-axis.

    sig2 : np.array
        Signal to plot using second (right) y-axis.
        
    fs : float or None
        Sampling frequency in [Hz].
        
    highlights : list of dict
        First (left) y-axis signal highlights, list of dict with keys:
            'epochs':
                1) iterable of 2 element iterable
                    Intervals to highlight, e.g. list of 2 element tuples (start, end).
                    If sampling frequency is supplied, itervals are measured in [s], otherwise, array index.
                2) 1D numpy.ndarray of bool type
                    Highlighting is applied to True array elements.
            'events': iterable of float
                If sampling frequency is provided, adds a marker to time points on the waveform. Otherwise, if sampling
                frequency is not provided, array indices are used.
            'color': optional, str, default yellow
                Use for same epoch highlight and event marker color, pyplot specifier.
            'epoch_color': optional, str, default yellow
                Epoch highlight color only. Overriden by 'color'.
            'event_color': optional, str, default yellow
                Event marker color only. Overriden by 'color'.
            'alpha': optional, float, default 0.5
                Epoch highlight and event marker transparency.
            'epoch_color': optional, float, default 0.5
                Epoch highlight alpha only. Overriden by 'alpha'.
            'event_color': optional, float, default 0.5
                Event marker alpha only. Overriden by 'alpha'.
            'event_marker_level': str or float
                Choose between 'min', 'max', a float or combination of both 'min' or 'max' and a float in the form 'min<+|-><offset>' to specify event marker level.
            'event_marker_size': float, default 3
                marker size in points.

    level_marker : float or None
        Adds a faint horizontal semi-transparent level marker to first (left) y-axis.

    crop : tuple of int or float
        If sampling frequency provided, start and end times in [s] to plot as (start_t, end_t), otherwise tuple of array indices.
        If waveform has no end boundary, pass None to the right crop tuple.
        
    title : str
        Plot title.
        
    ylabel : str
        Y-axis label.
        
    xlabel : str
        X-axis label. If fs is given, it defaults to "Time (s)".
    
    color : str
        Line color, pyplot keyword.

    marker : str
        Marker, pyplot keyword.

    linestyle : str
        Line style, pyplot keyword.
    
    linewidth : float
        Plot line width.
        
    figsize : tuple
        Figure size if no axis is provided.
        
    ax : matplotlib.axes._subplots.AxesSubplot
        If provided, will use axis to plot, otherwise will plot in a new figure.

    use_second_yaxis : bool
        If True, plot 'sig' on second y-axis that shares the same x-axis as the provided axis.

    legend_label : str
        Waveform legend label.
        
    xgrid : bool
        Plots vertical grid lines.
    
    """
    cropped_sig, (crop_left, crop_right) = proc.tslice(sig, fs, crop, return_intervals = True)
    if sig2 is not None:
        cropped_sig2 = proc.tslice(sig2, fs, crop, return_intervals = False)
    
    # cannot plot on second axis if second axis signal has been specified
    if use_second_yaxis and sig2 is not None:
        raise Exception("Cannot use second y-axis if second signal is specified.")

    if ax is None:
        fig, ax = plt.subplots(figsize = figsize, dpi = default_fig_dpi)
    else:
        if use_second_yaxis:
            ax = ax.twinx()
        fig = None
        
    if sig2 is not None:
        ax2 = ax.twinx()

    # add vertical grid if needed
    if xgrid:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.xaxis.grid(True, which = 'major', linestyle = '--', color = (0.5, 0.5, 0.5), linewidth = 1)
        ax.xaxis.grid(True, which = 'minor', linestyle = '--', color = (0.75, 0.75, 0.75), linewidth = 0.5)

    
    # add highlights
    for h in highlights:
        if 'color' in h:
            # if a single 'color' specifier is given, set same for event and for epoch
            epoch_color = h['color']
            event_color = h['color']
        else:
            # set epoch and event colors separately, use yellow by default
            epoch_color = h['epoch_color'] if 'epoch_color' in h else 'y'
            event_color = h['event_color'] if 'event_color' in h else 'y'

        if 'alpha' in h:
            # if a single 'alpha' specifier is given, set same for event and for epoch
            epoch_alpha = h['alpha']
            event_alpha = h['alpha']
        else:
            # set highlight and marker alpha separately, use 0.5 by default
            epoch_alpha = h['epoch_alpha'] if 'epoch_alpha' in h else 0.5
            event_alpha = h['event_alpha'] if 'event_alpha' in h else 0.5        

        if 'epochs' in h:
            highlight_waveform(h['epochs'], fs = fs, ax = ax, waveform = sig, crop = crop, color = epoch_color,
                alpha = epoch_alpha)
        if 'events' in h:
            if isinstance(h['event_marker_level'], numbers.Number):
                marker_level = h['event_marker_level']
            elif h['event_marker_level'] == 'min':
                marker_level = np.nanmin(cropped_sig)
            elif h['event_marker_level'] == 'max':
                marker_level = np.nanmax(cropped_sig)
            elif h['event_marker_level'].startswith('min'):
                marker_level = np.nanmin(cropped_sig)+float(h['event_marker_level'][3:])
            elif h['event_marker_level'].startswith('max'):
                marker_level = np.nanmax(cropped_sig)+float(h['event_marker_level'][3:])

            events_tp = proc.crop_timepoints(tp = h['events'], crop = crop)
            if events_tp:
                ax.plot(events_tp, [marker_level]*len(events_tp), 'o', c = event_color, markersize = h['event_marker_size'] if 'event_marker_size' in h else 3,
                    alpha = event_alpha)

    # plot first y-axis signal
    # ------------------------------------------------------------------------------
    # plot waveform
    t = np.linspace(crop_left, crop_right, len(cropped_sig))
    # add level marker
    if level_marker is not None:
        ax.axhline(level_marker, c = 'k', alpha = 0.2)
    ax.plot(t, cropped_sig, color = color, linestyle = linestyle,
        marker = marker, linewidth = linewidth, label = legend_label, alpha = alpha)
    # plot whole x-axis range regardless of starting or ending nans
    ax.set_xlim(t[0], t[-1])
    # plot second y-axis signal
    # ------------------------------------------------------------------------------
    if sig2 is not None:
        t2 = np.linspace(crop_left, crop_right, len(cropped_sig2))
        ax2.plot(t2, cropped_sig2, color = color2, linestyle = linestyle2,
        marker = marker2, linewidth = linewidth2, label = legend_label2, alpha = alpha2)
        # plot whole x-axis range regardless of starting or ending nans
        ax2.set_xlim(t2[0], t2[-1])

    if fs is not None:
        ax.set_xlabel('Time (s)',)
    else:
        if xlabel:
            ax.set_xlabel(xlabel,)
        
    # first y-axis label
    if ylabel != '':
        ax.set_ylabel(ylabel,)

    # second y-axis label
    if ylabel2 != '':
        ax.set_ylabel(ylabel2,)


    if title != '':
        ax.set_title(title,)
        
    ax.tick_params(axis = 'both', pad = default_tick_label_pad)
    if sig2 is not None:
        ax2.tick_params(axis = 'both', pad = default_tick_label_pad)
    
    return ax
    
def auto_corr(sig, fs, nlags, alpha = 0.05, linecolor = (0,0,1,1), bandcolor = (0,0,1,0.1), ylim = (None, None), ax = None, title = '',
    yticks = 5):
    """
    Plots the autocorrelation of a signal with confidence intervals.
    """
    out = tsa_stat.acf(sig, nlags = nlags, alpha = alpha , qstat = True, fft = True)
    # plot
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize = (10,3), facecolor = 'w')
    ax.plot(np.arange(0, len(out[0]))*1e3/fs, out[0], c = linecolor)
    ax.fill_between(np.arange(0, len(out[0]))*1e3/fs, out[1][:,0], out[1][:,1], edgecolor = bandcolor, facecolor = (bandcolor))
    ax.axhline(linestyle = '--', c = 'k')
    ax.set_xlabel('Time (ms)', fontsize = default_label_fontsize)
    ax.set_ylabel('Corr', fontsize = default_label_fontsize)
    ax.set_title(title, fontsize = default_title_fontsize)
    ax.tick_params(axis = 'both', labelsize = default_tick_label_fontsize, pad = default_tick_label_pad)
    ax.set_ylim(ylim)
    ax.locator_params(axis = 'y', tight = True, nbins = yticks)

    return (fig, ax)

def wavelet(sig, fs, wt_freq = (1, 500, 10), wt_ncycles = 5, t_slice = (0, None), ax = None, xlabel = True):
    """
    Plots wavelet transform.
    
    Parameters
    ----------
    sig : np.array
        Signal to analyze.
        
    fs : float
        Sampling frequency in [Hz].
        
    wt_freq : tuple
        Wavelet transform lower, upper and frequency bin as (f_lower, f_upper, f_bin).
        
    wt_ncycles : int
        Number of cycles in a wavelet at all frequencies. This is the tradeoff between temporal and frequency estimation.
        
    t_slice : tuple
        Start and end times in [s] from signal to plot as (start_t, end_t).
        
    ax : matplotlib.axes._subplots.AxesSubplot
        If provided, will use axis to plot, otherwise will plot in a new figure.
        
    xlabel : bool
        If True, "Time (s)" label is added to the x axis.
    """
    sig, (t_start, t_end) =  proc.tslice(sig, fs, t_slice, return_intervals = True)    
    
    # wavelet frequency bins
    wt_freqs = np.arange(*wt_freq)
    # generate wavelets
    w = mne.time_frequency.tfr.morlet(sfreq = fs, freqs = wt_freqs, n_cycles = wt_ncycles, sigma = None, zero_mean = False)
    # perform continuous wavelet transform with set of wavelets
    wt = np.reshape(mne.time_frequency.tfr.cwt(np.reshape(sig,(1, -1)), w, use_fft = True, mode = 'same', decim = 1), (len(wt_freqs), -1))
    # calculate power spectral density
    wt_psd = np.abs(wt)**2
    
    # plot log10 power spectrum
    if ax is None:
        fig, ax = plt.subplots(figsize = (15,5), facecolor = 'w', dpi = default_fig_dpi)
    else:
        fig = None
        
    ax.imshow(np.log10(wt_psd), extent = [t_start, t_end, wt_freq[0], wt_freq[1]], interpolation = 'nearest', origin = 'lower',
              aspect = 'auto')
    if xlabel:
        ax.set_xlabel('Time (s)', fontsize = default_label_fontsize)
    ax.set_ylabel('Frequency (Hz)', fontsize = default_label_fontsize)
    ax.tick_params(axis = 'both', labelsize = default_tick_label_fontsize, pad = default_tick_label_pad)
  
def cross_corr(sig1, sig2, nlag, fs, alpha = 0.05, stat_type = 'spearmanr', title = '', shuffle = False, bp_filter = None):
    """
    Plots the cross-correlation of two signals with statistical power at different lags by shifting (rolling) sig1
    w.r.t. sig2.
    
    Parameters
    ----------
    sig1, sig2 : numpy.array
        Simultaneously sampled signals.
    
    nlag : int
        Number of lags to perform shift; shift range is 2*nlag between [-nlag, nlag]. 
    
    fs : float
        Sampling frequency in [Hz].
        
    alpha : float
        Significance level. Lagged testing is done by adjusting the significance level according to Bonferroni's criterion.
            
    stat_type : str
        Type of cross-correlation statistic to use. Choose between 'pearsonr' and 'spearmanr'.
    
    title : str
        Plot title.
        
    shuffle : bool
        If True, also plots cross-correlation of sig1 with shuffled version of sig2.
        
    bp_filter : tuple
        Band-pass filter to apply to original and shuffled sig2 before cross-correlating with sig1 as (f_low, f_high) in [Hz].
    """
    if bp_filter is None:
        sig2_shuffled = np.random.permutation(sig2.copy())
    else:
        # shuffle original sig2 and then apply filter
        sig2_shuffled = proc.filt_FIR_bp(np.random.permutation(sig2.copy()), fs, bp_filter[0], bp_filter[1])
        # apply filter to original sig2
        sig2 = proc.filt_FIR_bp(sig2, fs, bp_filter[0], bp_filter[1])
              
    l = np.empty((2*nlag, 2))
    l_shuffled = np.empty((2*nlag, 2))
    for idx, s in enumerate(range(-nlag, nlag)):
        if -nlag+s:
            r_idx = -nlag+s
        else:
            r_idx = None 
        if stat_type == 'spearmanr':
            stat_func = stats.spearmanr
        elif stat_type == 'pearsonr':
            stat_func = stats.pearsonr
        else:
            raise Exception('Statistic not implemented/supported.')
        stat = stat_func(sig1[nlag:-nlag], sig2[nlag+s:r_idx])
        l[idx, 0] = stat[0]
        if stat[1]<=0:
            print(stat[1])     
        l[idx, 1] = stat[1] # p-value
        stat_shuffled = stat_func(sig1[nlag:-nlag], sig2_shuffled[nlag+s:r_idx])    
        l_shuffled[idx, 0] = stat_shuffled[0]
        l_shuffled[idx, 1] = stat_shuffled[1]
            
    # plot result
    if shuffle:
        norm = mpl.colors.LogNorm(vmin = min(alpha/(2*nlag)*1e-1, l[:,1].min(), l_shuffled[:,1].min()), vmax = alpha/(2*nlag))
    else:
        norm = mpl.colors.LogNorm(vmin = min(alpha/(2*nlag)*1e-1, l[:,1].min()), vmax = alpha/(2*nlag))
    mapper = mpl.cm.ScalarMappable(norm = norm, cmap = mpl.cm.jet_r)
    fig, ax = plt.subplots(figsize = (15,5), facecolor = 'w')
    ax.margins(x=0)
    ax.plot(np.arange(-nlag, nlag)*1e3/fs, l[:,0], '-', c = 'k')
    if shuffle:
        ax.plot(np.arange(-nlag, nlag)*1e3/fs, l_shuffled[:,0], '-.', c = (0.5, 0.5, 0.5))
    ax.set_xlabel('Time (ms)', fontsize = default_label_fontsize)
    ax.set_ylabel('Corr', fontsize = default_label_fontsize)
    ax.set_title(title, fontsize = default_title_fontsize)
    ax.tick_params(axis = 'both', labelsize = default_tick_label_fontsize, pad = default_tick_label_pad)
    ax.axvline(0, linestyle = '--', color = (0.5, 0.5, 0.5))
    cmap = mpl.cm.get_cmap('jet_r')
    
    for idx, s in enumerate(np.arange(-nlag, nlag)*1e3/fs):
        if l[idx,1] == 0:
            continue
        if l[idx,1] <= alpha/(2*nlag):
            sc = plt.scatter(s, l[idx,0], s = 50, c = mapper.to_rgba(l[idx,1]), marker = 'o', lw = 0)
        else:
            sc = plt.scatter(s, l[idx,0], s = 50, c = 'k', marker = '*', lw = 0)
        if shuffle:
            if l_shuffled[idx,1] <= alpha/(2*nlag):
                sc = plt.scatter(s, l_shuffled[idx,0], s = 50, c = mapper.to_rgba(l_shuffled[idx,1]), marker = 'o', lw = 0)
            else:
                sc = plt.scatter(s, l_shuffled[idx,0], s = 50, c = 'k', marker = '*', lw = 0, alpha = 0.5)    
    # add colorbar
    cax, _ = mpl.colorbar.make_axes(ax)
    cb = mpl.colorbar.ColorbarBase(cax, cmap = cmap, norm = norm)
    tick_locator = mpl.ticker.LogLocator()
    cb.locator = tick_locator
    cb.update_ticks()
    cb.ax.set_title('p-val', size = default_tick_label_fontsize)
    cb.ax.tick_params(labelsize = default_tick_label_fontsize)
    
def theta_LFP(lfp, fs, theta_band = (6, 10)):
    """
    Plots local-field potential (LFP) in the theta pass-band.
    
    Parameters
    ----------
    lfp : array_like
        Local field potential.
        
    fs : float
        Sampling frequency in [Hz].
        
    theta_band : tuple
        Theta band frequency limits as (low, high) tuple in [Hz].
    """
    filt_lfp = proc.filt_theta_band(lfp, fs, theta_band)
    t = np.arange(0, len(filt_lfp))/fs
    plt.figure(figsize = default_figsize, dpi = default_fig_dpi, facecolor = 'w', edgecolor = 'k')
    plt.plot(t, filt_lfp*1e3, '-b')
    plt.xlabel('Time (s)')
    plt.ylabel('LFP (mV)')
    plt.title('Theta-band LFP', fontsize = default_title_fontsize)
    plt.show()
    
def two_chan_fluorescence(data_path, sig, roi_name, ratiometric_fl = True, n_mov_avg = 1, max_fl_change = 0.5):
    """
    Plots fluorescence data from up to two channels for given ROI. 
    
    Parameters
    ----------
    data_path : str
        Path to folder containing data where plots will be saved in subfolder /analysis/fluorescence.pdf
    sig : dict
        Signals to plot on a common time axis. Dict must have at least the following keys:
           'fs' : float
              Common sampling frequency determined by the acquisition frame rate in [Hz].
           'ROIs' : dict
              Extracted fluorescence regions of interest (ROIs) as dict with keys denoting ROI name and values as dict with channel name keys and array_like values.         
    roi_name : str
        Name of ROI used to plot signals.
        
    ratiometric_fl : bool
        If True, plot also the ratiometric signal.
        
    n_mov_avg : int
        If > 1, apply n-point moving average.
    
    Returns
    -------
    fig : matplotlib.Figure  
    """
    # try to take number of samples from Ch1, otherwise from Ch2.
    if 'Ch1' in sig['ROIs'][roi_name]:
        nsamp = len(sig['ROIs'][roi_name]['Ch1'])
    elif 'Ch2' in sig['ROIs'][roi_name]:
        nsamp = len(sig['ROIs'][roi_name]['Ch2'])
    else:
        raise Exception("Fluorescence ROI '{}' does not contain 'Ch1' or 'Ch2' data".format(roi_name))
    # get number of channels
    nchan = sum([ch in sig['ROIs'][roi_name] for ch in ['Ch1', 'Ch2']])
    
    # time vector     
    t = np.arange(0, nsamp)/sig['fs']
    
    fig, axs = plt.subplots(1+nchan, 1, sharex = True, \
                            figsize = (default_figsize[0], default_figsize[1]*(1+nchan)), \
                             dpi = default_fig_dpi, facecolor = 'w', edgecolor = 'k')
    
    plt_idx = 0
    
    processed_sig = {}
    if ratiometric_fl:
        alpha = 0.5
    else:
        alpha = 1
    
    for ch in [('Ch1', 'g'), ('Ch2', 'r')]:
        if ch[0] in sig['ROIs'][roi_name]:
            # apply only moving average if needed, but leave bleaching unchanged
            processed_sig[ch[0]] = proc.mov_avg(sig['ROIs'][roi_name][ch[0]], n_mov_avg)
            # apply bleaching correction
            processed_sig[ch[0]+' debleached'] = proc.debleach(processed_sig[ch[0]], maxfev = 50000)      
            # plot only moving averaged signal    
            axs[plt_idx].plot(t[:len(processed_sig[ch[0]])], processed_sig[ch[0]], ch[1], label = ch[0], linewidth = default_linewidth)
            axs[plt_idx].set_ylabel(ch[0], fontsize = default_label_fontsize)
            axs[plt_idx].grid(which='major', axis='x', linestyle='--')
            axs[plt_idx].tick_params(axis='y', labelsize = default_tick_label_fontsize)
            axs[plt_idx].locator_params(axis = 'y', nbins = default_nyticks)
            plt_idx += 1    
    
    # plot bleaching corrected signals
    if processed_sig:
        for ch in [('Ch1', 'g'), ('Ch2', 'r')]:
            axs[plt_idx].plot(t[:len(processed_sig[ch[0]+' debleached'])], processed_sig[ch[0]+' debleached'], ch[1], label = ch[0], alpha = alpha, \
                              linewidth = default_linewidth)
    
    if 'Ch1' in processed_sig and 'Ch2' in processed_sig and ratiometric_fl:
        ratio_sig = processed_sig['Ch1 debleached']/processed_sig['Ch2 debleached']
        axs[plt_idx].plot(t[:len(ratio_sig)], ratio_sig, '-k', label = 'Ch1/Ch2', linewidth = default_linewidth)
        
        #axs[plt_idx].set_ylim([1-max_fl_change, 1+max_fl_change])
        axs[plt_idx].legend(loc = "upper right", fontsize = default_legend_fontsize)
        axs[plt_idx].set_ylabel('Norm. fluorescence', fontsize = default_label_fontsize)
        axs[plt_idx].grid(which='major', axis='x', linestyle='--')
        axs[plt_idx].locator_params(axis = 'y', nbins = default_nyticks)
     
    # adjust time axis label and ticks
    axs[plt_idx].set_xlabel('Time (s)', fontsize = default_label_fontsize)
    axs[plt_idx].tick_params(axis='both', labelsize = default_tick_label_fontsize)  
    # adjust and show figure
    fig.suptitle("ROI {}".format(roi_name)+' fluorescence', fontsize = default_title_fontsize)
    fig.tight_layout()
    fig.subplots_adjust(top=0.95)
    fig.show()
    # save figure
    fig_filename = os.path.join(data_path, "analysis/ROI {} fluorescence.pdf".format(roi_name))
    if not os.path.exists(os.path.dirname(fig_filename)):
        try:
            os.makedirs(os.path.dirname(fig_filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise

    fig.savefig(fig_filename, orientation = 'landscape', papertype = 'letter')   

    return fig
    
def aligned_signals(data_path, sig, ratiometric_fl = True, n_mov_avg = 1, max_fl_change = None):
    """
    Plots time-aligned signals.
    
    Parameters
    ----------
    data_path : str
        Path to folder containing data where plots will be saved in subfolder /analysis/aligned_ROI_LFP_pos.pdf
    sig : dict
        Signals to plot on a common time axis. Dict must have the following keys:
          'fs' : float
              Common sampling frequency determined by the acquisition frame rate in [Hz].
          'pos' : array_like
              Wrapped belt position in [mm].
          'LFP' : array_like
              Local field potential in [V].
          'ROIs' : dict
              Extracted fluorescence regions of interest (ROIs) as dict with keys denoting ROI name and values as dict with channel name keys and array_like values.    
      
    ratiometric_fl : bool
        If True, plot ratiometric Ch1/Ch2 signal.
        
    n_mov_avg : int
        If > 1, a moving n-sample average is applied to all signals. Note that the moving average filter will both delay and reduce the amplitude of signals.
        
    max_fl_change : float
        Maximum relative fluorescence change used for adjusting y-axis plotting range.
        
    Returns
    -------
    fig : matplotlib.Figure
    """ 
    # determine if behavior plot is added
    use_belt_position = 'pos' in sig
    
    # use number of LFP elements to construct time vector
    t = np.arange(0, len(sig['LFP']))/sig['fs']
    fig, axs = plt.subplots(1+1+use_belt_position+len(sig['ROIs']), 1, sharex = True, \
                            figsize = (default_figsize[0], default_figsize[1]*(1+1+use_belt_position+len(sig['ROIs']))), \
                             dpi = default_fig_dpi, facecolor = 'w', edgecolor = 'k')
    
    plt_idx = 0
    #==== plot wrapped belt position ====
    if use_belt_position:
      belt_pos_sig = proc.mov_avg(sig['pos'], n_mov_avg)
      axs[plt_idx].plot(t[:len(belt_pos_sig)], belt_pos_sig, '-b')
      axs[plt_idx].set_title('Treadmill', fontsize = default_title_fontsize)
      axs[plt_idx].set_ylabel('Belt position (mm)', fontsize = default_label_fontsize)
      axs[plt_idx].grid(which='major', axis='x', linestyle='--')
      plt_idx += 1
      
    #==== plot LFP ====
    lfp_sig = proc.mov_avg(sig['LFP']*1e3, n_mov_avg)  
    axs[plt_idx].plot(t[:len(lfp_sig)], lfp_sig, '-b')
    axs[plt_idx].set_title('LFP', fontsize = default_title_fontsize)
    axs[plt_idx].set_ylabel('Amplitude (mV)', fontsize = default_label_fontsize)
    axs[plt_idx].tick_params(axis='y', labelsize = 12)
    axs[plt_idx].grid(which='major', axis='x', linestyle='--')
    plt_idx += 1
    
    #==== plot LFP wavelet ====
    plot_wavelet(sig['LFP'], sig['fs'], ax = axs[plt_idx], xlabel = False)
    axs[plt_idx].grid(which = 'major', axis = 'x', linestyle = '--')
    plt_idx += 1
    
    #==== plot ROIs ====
    for roi_group in sig['ROIs']:
        for roi_name in sig['ROIs'][roi_group]:
            fl_sig_found = False
            ch1_sig = None
            ch2_sig = None
            if ratiometric_fl:
                alpha = 0.5
            else:
                alpha = 1
                
            if 'Ch1' in sig['ROIs'][roi_group][roi_name]:
                ch1_sig = proc.mov_avg(proc.debleach(sig['ROIs'][roi_group][roi_name]['Ch1'], maxfev = 50000), n_mov_avg)    
                axs[plt_idx].plot(t[:len(ch1_sig)], ch1_sig, '-g', label = 'Ch1', alpha = alpha)
                fl_sig_found = True
            if 'Ch2' in sig['ROIs'][roi_group][roi_name]:
                ch2_sig = proc.mov_avg(proc.debleach(sig['ROIs'][roi_group][roi_name]['Ch2'], maxfev = 50000), n_mov_avg) 
                axs[plt_idx].plot(t[:len(ch2_sig)], ch2_sig, '-r', label = 'Ch2', alpha = alpha)
                fl_sig_found = True
                
            if ch1_sig is not None and ch2_sig is not None and ratiometric_fl:
                ratio_sig = ch1_sig/ch2_sig
                axs[plt_idx].plot(t[:len(ratio_sig)], ratio_sig, '-k', label = 'Ch1/Ch2')
                    
            if not fl_sig_found:
                raise Exception("No fluorescence channels 'Ch1' or 'Ch2' were found in the ROIs signals.")
            
            if max_fl_change is not None:
                axs[plt_idx].set_ylim([1-max_fl_change, 1+max_fl_change])
                
            axs[plt_idx].tick_params(axis='y', labelsize = 12)
            axs[plt_idx].legend(loc = "upper right", fontsize = default_legend_fontsize)
            axs[plt_idx].set_title("ROI {}/{}".format(roi_group, roi_name)+' fluorescence', fontsize = default_title_fontsize)
            axs[plt_idx].set_ylabel('Fluorescence', fontsize = default_label_fontsize)
            axs[plt_idx].grid(which='major', axis='x', linestyle='--')
            
            plt_idx += 1
             
    axs[plt_idx-1].set_xlabel('Time (s)', fontsize = default_label_fontsize)
    axs[plt_idx-1].tick_params(axis='x', labelsize = 14)  
    # change fig size for printing and save
    #fig.set_size_inches(default_print_fig_width, default_figsize[1]*(1+use_belt_position+len(sig['ROIs'])))
    #fig.tight_layout()
    fig.savefig(os.path.join(data_path, "analysis/aligned_ROI_LFP_pos.pdf"), orientation = 'landscape', papertype = 'letter',
                dpi = fig.dpi)
    # change back figure size for display
    #fig.set_size_inches(default_figsize[0], default_figsize[1]*(1+use_belt_position+len(sig['ROIs'])))
    #fig.tight_layout()
    #fig.show()
    
    return fig
       
def spectral_coherence(sig1, sig2, fs, sig1_name = 'signal 1', sig2_name = 'signal 2', nFFT = 256):
    """
    Plot spectral coherence of two signals.
    
    Parameters
    ----------
    sig1, sig2 : array_like
        Signals. Signals need not be of equal length, analysis is carried out only for the common portion.
        
    fs : float
        Sampling frequency in [Hz].
        
    sig1_name, sig2_name : str
        Name of signals 1 & 2.
    """
    nsamp = min(len(sig1), len(sig2))
    t = np.arange(0, nsamp)/fs
    
    fig, axs = plt.subplots(2, 1, figsize = (default_figsize[0], default_figsize[1]*2), \
                            dpi = default_fig_dpi, facecolor = 'w', edgecolor = 'k')
    axs[0].plot(t, sig1, '-b', label = sig1_name)
    axs[0].plot(t, sig2, '-g', label = sig2_name) 
    axs[0].set_xlabel('Time (s)', fontsize = default_label_fontsize)
    axs[0].set_ylabel('Signals', fontsize = default_label_fontsize)
    axs[0].grid(which='major',axis ='x', linestyle = '--')
    axs[0].legend(loc = "upper right", fontsize = default_legend_fontsize)
    
    cxy, f = axs[1].cohere(sig1[:nsamp], sig2[:nsamp], nFFT, fs)
    axs[1].set_ylabel('Coherence')
    
    fig.show()
        
def theta_phase_bin_average(sig, lfp, fs, nbins = 15, sig_bp = (2, 20)):
    """
    Plots a phase-bin average of provided signals referenced to the theta band-passed local field potential. 
    
    Parameters
    ----------
    sig : list of 1D numpy.ndarray
        Signals used for averaging.
            
    lfp : list of 1D numpy.ndarray
        Local field potential used to determine theta-band oscillation throughs for averaging.
        
    
        
    """
    # theta-band oscillation through times
    throughs = [an.proc.get_theta_throughs(s, fs) for s in sig]
    
    phase_bins = an.proc.phase_bin_avg2([an.proc.filt_theta_band(sig['ROIs']['dendrite']['d1']['Ch1/Ch2'], sig['fs'], theta_band = (2, 20))],
                                           sig['fs'], throughs, nbins = nbins)
    plt.figure(figsize= (5,5))
    phase_bins = np.array(phase_bins*2)
    plt.plot(np.arange(phase_bins.shape[0]), phase_bins[:,0])
    plt.fill_between(np.arange(phase_bins.shape[0]), phase_bins[:,1], phase_bins[:,2], color = '#539caf', alpha = 0.4, label = '95% CI')

def _regression_fit_ci_manual(t, s_err, n, x, x2, y2, fill_color = "#b9cfe7", ax = None):
    """
    Return an axes of confidence bands using a simple approach.

    Notes
    -----
    .. math:: \left| \: \hat{\mu}_{y|x0} - \mu_{y|x0} \: \right| \; \leq \; T_{n-2}^{.975} \; \hat{\sigma} \; \sqrt{\frac{1}{n}+\frac{(x_0-\bar{x})^2}{\sum_{i=1}^n{(x_i-\bar{x})^2}}}
    .. math:: \hat{\sigma} = \sqrt{\sum_{i=1}^n{\frac{(y_i-\hat{y})^2}{n-2}}}

    References
    ----------
    .. [1] M. Duarte.  "Curve fitting," Jupyter Notebook.
       http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/CurveFitting.ipynb

    """
    if ax is None:
        ax = plt.gca()

    ci = t * s_err * np.sqrt(1/n + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))
    ax.fill_between(x2, y2 + ci, y2 - ci, color = fill_color, edgecolor = "")

    return ax

def _regression_fit_ci_bootstrap(xs, ys, order, resid, nboot = 500, ax = None):
    """
    Return an axes of confidence bands using a bootstrap approach.

    Notes
    -----
    The bootstrap approach iteratively resampling residuals.
    It plots `nboot` number of straight lines and outlines the shape of a band.
    The density of overlapping lines indicates improved confidence.

    Returns
    -------
    ax : axes
        - Cluster of lines
        - Upper and Lower bounds (high and low) (optional)  Note: sensitive to outliers

    References
    ----------
    .. [1] J. Stults. "Visualizing Confidence Intervals", Various Consequences.
       http://www.variousconsequences.com/2010/02/visualizing-confidence-intervals.html

    """ 
    if ax is None:
        ax = plt.gca()

    bootindex = sp.random.randint

    for _ in range(nboot):
        resamp_resid = resid[bootindex(0, len(resid) - 1, len(resid))]
        # Make coeffs of for polys
        pc = sp.polyfit(xs, ys + resamp_resid, deg = order)                   
        # Plot bootstrap cluster
        ax.plot(xs, sp.polyval(pc, xs), "b-", linewidth = 2, alpha = 3.0 / float(nboot))

    return ax

def _data_picked_event(event, labels = None):
    """
    Print info about the closest datapoint to the mouse clicked position on the plot.

    Parameters
    ----------
    event :
        Matplotlib event picker event.
    labels : dict of list of str
        Datapoint labels as dict with keys being assigned Artist ID and values list of str with
        same order as plotted data points.
    """
    # pick closest data point to mouse click location
    ind = event.ind
    datax,datay = event.artist.get_data()
    if len(ind) > 1:
        _datax,_datay = [datax[i] for i in ind],[datay[i] for i in ind]
        msx, msy = event.mouseevent.xdata, event.mouseevent.ydata
        dist = np.sqrt((np.array(_datax)-msx)**2+(np.array(_datay)-msy)**2)
        ind = ind[np.argmin(dist)]
    else:
        ind = ind[0]

    gid = event.artist.get_gid()
    if gid is None or not labels:
        print('\ndatapoint (x,y) = {},{} picked: no info is available.'.format(datax[ind], datay[ind]))
        return
    print(ind)
    print('\ndatapoint (x,y) = {},{} picked:\n{}'.format(datax[ind], datay[ind], labels[gid][ind]))
    
def errorbars_with_regression_fit(x, y, gid = None, order = 1, ax = None, xerr = None, yerr = None,
    ebar_color = (0.5,0.5,0.5,0.5), marker_facecolor = 'k', marker_edge_color = 'k', scatter_alpha = None,
    marker_size = default_plt_markersize, label = None, fit_ci = 'manual'):
    """
    Plots a polynomial regression fit.

    Parameters
    ----------
    x : 1D numpy.ndarray
        Independent variable.

    y : 1D numpy.ndarray
        Dependent variable.

    gid: int
        User assigned Artist group ID. Use this to look up datapoint labels on user click.

    order : int
        Polynomial fit order.

    ax : matplotlib.axes._subplots.AxesSubplot
        Use provided axes for plotting, otherwise plot in a new figure.

    xerr, yerr : 2D numpy.ndarray
        Independent and dependent data error bars. First row, lower limit, second row upper limit.

    ebar_color : tuple or str
        Color of error bars as (R,G,B,A) components ranging [0,1] or pyplot color name.

    marker_facecolor : tuple or str
        Color of marker face as (R,G,B,A) components ranging [0,1] or pyplot color name.

    marker_edge_color : tuple or str
        Color of marker edge as (R,G,B,A) components ranging [0,1] or pyplot color name.

    scatter_alpha : float
        Alpha value of scatter markers and their error bars.

    marker_size : float
        Scatter plot marker size in points.

    label : str or None
        Dataset label.

    fit_ci : str
        Fit confidence interval type. Choose from 'manual' and 'bootstrap'.
    """
    # holds regression statistics
    stat_out = {}

    if ax is None:
        fig, ax = plt.subplots(figsize = (5,5), facecolor = 'w', dpi = default_fig_dpi)
    else:
        fig = ax.get_figure()

    def poly_func(a, b):
        """Return a 1D polynomial."""
        return np.polyval(a, b)

    # scatter plot with error bars
    # --------------------------------------------------------------------------------------------------------------------------------------------
    ax.errorbar(x, y, yerr = yerr, xerr = xerr, fmt = 'ok', ecolor = ebar_color, markersize = marker_size,
        markerfacecolor = marker_facecolor, markeredgecolor = marker_edge_color, alpha = scatter_alpha,
        label = label, picker = 5, gid = gid)

    # add regression fit if there is enough data
    # --------------------------------------------------------------------------------------------------------------------------------------------
    if len(x)<3:
        util.clrd_print("Regression fit skipped.\n", 'warn')
        return fig, ax, stat_out

    # fit polynomial to data
    p, cov = np.polyfit(x, y, deg = order, cov = True)         # parameters and covariance from the fit of 1-D polynom.
    y_model = poly_func(p, x)                                  # model using the fit parameters; NOTE: parameters here are coefficients

    # compute statistics
    stat_out['n'] = y.size                                     # number of observations
    stat_out['m'] = p.size                                     # number of parameters
    stat_out['dof'] = stat_out['n'] - stat_out['m']            # degrees of freedom
    stat_out['t'] = stats.t.ppf(0.975, stat_out['dof'])        # used for CI and PI bands, correction for two sided 95% confidence t-statistic

    # estimate error in data and model
    resid = y - y_model                           
    stat_out['chi2'] = np.sum((resid / y_model)**2)            # chi-squared; estimates error in data
    stat_out['chi2_red'] = stat_out['chi2'] / stat_out['dof']  # reduced chi-squared; measures goodness of fit
    s_err = np.sqrt(np.sum(resid**2) / stat_out['dof'])        # standard deviation of the error

    
    # plot fit
    ax.plot(x, y_model, "-", color = (0, 0, 0, 0.5), linewidth = 1, label = label+' fit' if label is not None else "fit")
    # plot fit confidence interval
    x2 = np.linspace(np.min(x), np.max(x), 100)
    y2 = poly_func(p, x2)
    if fit_ci == 'manual':
        _regression_fit_ci_manual(stat_out['t'], s_err, stat_out['n'], x, x2, y2, fill_color = (0.75, 0.75, 0.75, 0.5), ax = ax)
    elif fit_ci == 'bootstrap':
        _regression_fit_plot_ci_bootstrap(x, y, order, resid, ax = ax)
    else:
        raise ValueError("Invalid fit confidence interval method '{}'.\n".format(fit_ci))

    # plot prediction interval
    pi = stat_out['t'] * s_err * np.sqrt(1 + 1/stat_out['n'] + (x2 - np.mean(x))**2 / np.sum((x - np.mean(x))**2))   
    ax.fill_between(x2, y2 + pi, y2 - pi, color = "None", linestyle = "--")
    ax.plot(x2, y2 - pi, "--", color = "0.75", label = "95% Prediction Limits")
    ax.plot(x2, y2 + pi, "--", color = "0.75")
    
    return fig, ax, stat_out

def histplot(data, x = None, bins = None, kde = False, shared_bins = True, hue = None, vertical = False,
    palette = RANDOM_COLORS, ax = None, legend = False, **seaborn_kws):
    """Extends sb.distplot to work with dataframes, hues 
    
    Parameters
    ----------
    data : dataframe, wide-form
        Wide-form dataframe where `x` is a subset of the columns
    x : str, optional
        Column of dataframe to plot as histograms
    bins : int or bins, optional
        If int, all data will be plotted on same bins unless 
        `shared_bins` is False
    shared_bins : bool, optional
        Whether to plot different labels on same auto-generated bins. 
        Defaults to True.
    ax : Axes instance, optional
        Axes to plot on
    hue : str, column of data
        Column to color data by
    vertical : bool
        If True, observed values are on y-axis
    palette : dict or palette, optional
        Colors to use by hue
    legend : bool
        Whether to include a legend
    **seaborn_kws 
    
    Returns
    -------
    ax : Axes instance
        Axes plotted on
    """

    if x is None:
        hist_data = data
    else:
        hist_data = data[x]
        
    if bins is None:
        bins = 'auto'

    if shared_bins:
        bins = np.histogram_bin_edges(hist_data, bins = bins)
    
    if ax is None:
        fig, ax = plt.subplots()

    # if kde is True:
    #     seaborn_kws['norm_hist'] = True
    #     if 'kde_kws' in seaborn_kws:
    #         kde_kws = seaborn_kws['kde_kws']
    #     else:
    #         kde_kws = {}
    #     # if 'cumulative' in seaborn_kws:
    #     #     kde_kws['cumulative'] = seaborn_kws['cumulative']
    #     densityplot(data=data, x=x, hue=hue, palette=palette, 
    #                 ax=ax, **kde_kws)

    
    if hue is None:
        sb.distplot(hist_data, bins = bins, ax = ax, 
                    vertical = vertical,
                    **seaborn_kws)
    else:
        for hue_key, hue_data in hist_data.groupby(hue):
            sb.distplot(hue_data, bins = bins, ax = ax, 
                         kde = kde,
                         color = palette[hue_key], 
                         label = hue_key,
                         vertical = vertical,
                         **seaborn_kws)

    if legend:
        ax.legend()
    
    return ax

def barplot(df, stacked = False, max_groups = 15, columns = [], xlabel = '', ylabel = '', title = '',
    figsize = (9,3), legend_labels = {}, legend_loc = 'best', legend_ncol = 1, legend_title = '',
    xtick_label_fontsize = None):
    """
    Plots columns of a pandas DataFrame as barplot on multiple figures.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with columns to plot. 
    
    stacked : bool
        If True, stack barplots from selected columns.

    max_groups : int
        Maximum number of barplot groups per figure. If len(df) > max_groups, additional barplots
        are placed on multiple figures.

    columns : list of str
        Column names to include.

    xlabel, ylabel : str
        X and Y axis labels.

    figsize : tuple
        Figure size as tuple (width, height) in inches.

    legend_labels : dict
        Use custom legend names for selected dataframe columns. Pass as {'<old col name>': '<new col name>'}.

    legend_loc : str
        Legend location, see matplotlib specs.

    legend_ncol : int
        Number of legend columns.

    legend_title : str
        Legend title.

    xtick_label_fontsize : float
        X-axis tick label fontsize.

    Returns
    -------
    list of matplotlib Figure
    """
    out = []
    nfigs = math.ceil(len(df)/max_groups)
    # select columns
    if columns:
        df = df[columns]
    # rename columns e.g. to have nice short names for legend
    if legend_labels:
        df.rename(columns = legend_labels, inplace = True) 
    for fig_idx in range(nfigs):
        fig, ax = plt.subplots(1, 1, figsize = tuple(figsize))
        out.append(fig)
        # fit max rows per figure
        df_sel = df.iloc[fig_idx*max_groups:min((fig_idx+1)*max_groups, len(df))]
        # stacked barplot
        df_sel.plot(ax = ax, kind = 'bar', stacked = stacked)
        # rotate index x-axis labels
        ax.set_xticklabels(ax.get_xticklabels(), rotation = 30, horizontalalignment = "center")
        # change x-tick label font size
        if xtick_label_fontsize is not None:
            ax.tick_params(axis = "x", labelsize = xtick_label_fontsize)
            
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title+' (plot %d/%d)'%(fig_idx+1, nfigs))
        # adjust legend
        ax.legend(loc = legend_loc, ncol = legend_ncol, title = legend_title, framealpha = 0.5)

    return out

def apply_category_plot_grid1(fig, ax, par):
    """
    Applies a minimalist category grid plot layout by keeping spines and ticks only for left and bottom most plot edges.
    The subplot grid has:
        - same column-wise and row-wise subplot labels.
        - common x-axis label.

    Parameters
    ----------
    fig: matplotlib.figure.Figure
        Figure object.
    ax : 2D iterale of matplotlib.axes._subplots.AxesSubplot
        Subplot axes grid with:
            axis 0 = subplot row
            axis 1 = subplot column
    par : dict
        Plotting parameters, dict with keys:
            'suptitle' : optional, str
                Adds a title to the whole figure.
            'col_labels' : optional, list of str
                Column-wise labels. Default, no labels.
            'row_labels' : optional, list of str
                Row-wise labels. Default, no labels.
            'shared_x_label' : optional, str
                Shared x-axis label for all columns. Default, no label.
            'shared_y_label' : optional, str
                Shared y-axis label for all rows. Default, no label.
            'shared_x_label_pad', 'shared_y_label_pad' : optional, float
                Extra padding added to shared axis labels so they don't overlap with row-wise labels or with tick labels.
                default shared_x_label_pad = 10
                default shared_y_label_pad = 10 if there are no row labels, otherwise doubles this ammount.
            'legend' : optional, list of dict:
                Add legend to one or more subplots. List of dict with keys:
                    'subplot_id' : tuple of 2 elements
                        Row and column 0-index of subplot to which to attach legend.
                    'loc' : str
                        Legend location, e.g. 'upper right'. Default 'best'.
                    'ncol' : optional, int
                        Number of columns to use for legend, default 1.
    """
    # number of rows
    nrows = len(ax)
    ncols = [len(x) for x in ax]
    assert all([x == ncols[0] for x in ncols]) # to be sure this is a rectangular grid
    # number of columns
    ncols = ncols[0]
    # default parameters
    util.set_default_keys( {
        'shared_x_label_pad': 10,
        'shared_y_label_pad': 10,
        }, par)

    # set subplot column-wise labels
    if 'col_labels' in par:
        for col_idx, col_label in enumerate(par['col_labels']):
            ax[0][col_idx].set_title(col_label)
    
    # set subplot row-wise labels
    if 'row_labels' in par:
        for row_idx, row_label in enumerate(par['row_labels']):
            ax[row_idx][0].set_ylabel(row_label)
        par['shared_y_label_pad'] *= 2
    
    # add shared x- and y-axis labels
    if 'shared_x_label' in par or 'shared_y_label' in par:
        frame_ax = fig.add_subplot(111, frame_on = False)
        frame_ax.tick_params(axis = 'both', which = 'both', labelcolor = "none", bottom = False, left = False, right = False, top = False)

    if 'shared_x_label' in par:
        frame_ax.set_xlabel(par['shared_x_label'], labelpad = par['shared_x_label_pad'])
    if 'shared_y_label' in par:
        frame_ax.set_ylabel(par['shared_y_label'], labelpad = par['shared_y_label_pad'])
    
    # add legend
    if 'legend' in par:
        for l in par['legend']:
            ax[l['subplot_id'][0]][l['subplot_id'][1]].legend(loc = l['loc'] if 'loc' in l else 'best',
                ncol = l['ncol'] if 'ncol' in l else 1)
    
    # add suptitle
    if 'suptitle' in par:
        fig.suptitle(par['suptitle'], y = 1)

    # show only left spines with ticks for first column and bottom spines with ticks for last row
    spine_locations = ['top', 'right', 'bottom', 'left']
    for row_idx in range(nrows):
        for col_idx in range(ncols):
            for spine_loc in spine_locations:
                ax[row_idx][col_idx].spines[spine_loc].set_visible(False)
                ax[row_idx][col_idx].tick_params(top = 'off', right = 'off', bottom = 'off', left = 'off')
    for row_idx in range(nrows):
        ax[row_idx][0].spines['left'].set_visible(True)
        ax[row_idx][0].tick_params(left = 'on')
    for col_idx in range(ncols):
        ax[1][col_idx].spines['bottom'].set_visible(True)
        ax[1][col_idx].tick_params(bottom = 'on')

    # adjust bottom axis position to detach it from y axis
    for col_idx in range(ncols):
        ax[1][col_idx].spines['bottom'].set_position(('axes', -0.05))

def add_combined_legend(fig, **kwargs):
    """
    Combines legends from multiple axes in the same figure.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Figure object.

    **kwargs :
        Variable named arguments passed to plt.legend()

    Returns
    -------
    None
    """
    handles, labels = [],[]
    for ax in fig.axes:
        for h,l in zip(*ax.get_legend_handles_labels()):
            handles.append(h)
            labels.append(l)
    plt.legend(handles,labels, **kwargs)

def sort_legend_labels(ax):
    """
    Sorts order of legend entries by label name.

    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
        Plot axis.

    Returns
    -------
    tuple
        (handles, labels)
        Pass on to ax.legend()
    """
    handles, labels = ax.get_legend_handles_labels()
    # sort both labels and handles by labels
    if handles and labels:
        labels, handles = zip(*sorted(zip(labels, handles), key = lambda t: t[0]))
    return handles, labels