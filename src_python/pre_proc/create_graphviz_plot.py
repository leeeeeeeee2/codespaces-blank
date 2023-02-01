#!/usr/bin/env python
# encoding: utf-8

"""
File Name   : create_graphviz_plot
Project Name: MPR
Description : creates a graphviz visualization of data array interdependencies as read from namelist mpr.nml
Author      : Robert Schweppe (robert.schweppe@ufz.de)
Created     : 2019-03-08 13:30
"""

# IMPORTS
import argparse
from itertools import zip_longest
import f90nml
from graphviz import Digraph, ENGINES, FORMATS
import warnings
import pathlib as pl
from textwrap import dedent
from src_python.pre_proc.mpr_interface import OPTIONS
import re

DOT2TEX_IS_AVAILABLE = True
try:
    import dot2tex as d2t
except ImportError:
    d2t = None
    DOT2TEX_IS_AVAILABLE = False

DOT2TEX_IS_AVAILABLE = True
try:
    import dot2tex as d2t
except ImportError:
    d2t = None
    DOT2TEX_IS_AVAILABLE = False

# GLOBAL VARIABLES
OUTPUT_FILE_BASENAME = 'MPR_graphviz'
DOT2TEX_GUIDE_HTML = 'https://dot2tex.readthedocs.io/en/latest/installation_guide.html'
DEFAULT_PATH_TO_ROOT = '../../'
DEFAULT_CONFIG_FILE = 'predefined_nmls/mpr_mhm_visual.nml'
# DEFAULT_PARAM_FILE = 'mpr_global_parameter.nml' # NOT USED
DEFAULT_OUTPUT_FORMAT = 'pdf'
DEFAULT_ENGINE = 'dot'
# TODO: this is broken?
DEFAULT_TRANSFER_ONLY = False
DEFAULT_USE_NODE_INDEX = False
DEFAULT_USE_DOT2TEX = False
# TODO: this is broken?
DEFAULT_GROUP_EDGES = True
DEFAULT_BLANK = False
FONT = 'Arial'
COORD_STYLE_OFFSET = 5
NODE_STYLE_NAMES = ['da_node_from_to_file', 'da_node_to_file', 'da_node_from_file', 'da_node']
NODE_STYLE_COLORS = ['red!80', 'red!80', 'blue!80', 'black!80']
EDGE_STYLE_NAMES = ['tf', 'uo', 'tf', 'uo_tf']
EDGE_STYLE_PENS = ['solid,black!20', 'solid,black!20', 'solid,black!20', 'solid,black!20']
NODE_COLORS_DICT = {k: v for k, v in zip(NODE_STYLE_NAMES, NODE_STYLE_COLORS)}
LEGEND_MARKER = 'LEGEND_MARKER'

# FUNCTIONS


# CLASSES
class DataArray(object):
    """basic container for all relevant information on data array"""
    # how we draw the edges based on DataArray properties
    # keys: (transfer_func is None, target_coord_names is None, use_dot2tex)
    DEFAULT_EDGE_ATTRS = {
        # this first combination should technically not happen, all DataArrays shall either use UO or TF
        (True, True, True): {'style': EDGE_STYLE_NAMES[0]},
        (True, False, True): {'style': EDGE_STYLE_NAMES[1]},
        (False, True, True): {'style': EDGE_STYLE_NAMES[2]},
        (False, False, True): {'style': EDGE_STYLE_NAMES[3]},
        # this first combination should technically not happen, all DataArrays shall either use UO or TF
        (True, True, False): {'style': EDGE_STYLE_PENS[0]},
        (True, False, False): {'style': EDGE_STYLE_PENS[1]},
        (False, True, False): {'style': EDGE_STYLE_PENS[2]},
        (False, False, False): {'style': EDGE_STYLE_PENS[3]},
    }
    # how we draw the nodes based on DataArray properties
    # keys: (to_file, from_file, use_dot2tex)
    DEFAULT_NODE_ATTRS = {
        (True, True, True): {'style': NODE_STYLE_NAMES[0]},
        (True, False, True): {'style': NODE_STYLE_NAMES[1]},
        (False, True, True): {'style': NODE_STYLE_NAMES[2]},
        (False, False, True): {'style': NODE_STYLE_NAMES[3]},
        (True, True, False): {'shape': 'ellipse', 'color': NODE_STYLE_COLORS[0].split('!')[0]},
        (True, False, False): {'shape': 'ellipse', 'color': NODE_STYLE_COLORS[1].split('!')[0]},
        (False, True, False): {'shape': 'ellipse', 'color': NODE_STYLE_COLORS[2].split('!')[0]},
        (False, False, False): {'shape': 'ellipse', 'color': NODE_STYLE_COLORS[3].split('!')[0]},
    }
    # how we draw the tf nodes based on transfer func properties
    # keys: use_dot2tex
    DEFAULT_TF_NODE_ATTRS = {
        True: {'style': 'tf_node'},
        False: {'shape': 'rectangle', 'color': 'grey'},
    }

    def __init__(self, name, from_file=None, to_file=None, target_coord_names=None, from_data_arrays=None,
                 transfer_func_label=None, transfer_func=None):
        self.name = name
        self.from_file = from_file
        # if None is set, this defaults to True in MPR
        if to_file is None:
            self.to_file = True
        else:
            self.to_file = to_file
        if target_coord_names is None:
            self.target_coord_names = []
        else:
            self.target_coord_names = [_ for _ in target_coord_names if _ is not None]
        if from_data_arrays is None:
            # we set an empty list for cases without from_data_arrays
            self.from_data_arrays = []
        else:
            # we only want the valid entries (None is inserted due to namelist structure)
            self.from_data_arrays = [_ for _ in from_data_arrays if _ is not None]
        self.transfer_func_label = transfer_func_label
        self.transfer_func = transfer_func
        # we do not want the appear the *.*_bound variables on the diagram as they are inferred from dimensions
        if self.transfer_func is not None and self.transfer_func.endswith(('.lower_bound', '.upper_bound')):
            self.transfer_func = None
            self.from_data_arrays = []

    def update_node_attrs(self, node_attrs, use_dot2tex=False):
        """return dict for node_attrs based on from_file and to_file"""
        node_attrs = dict(node_attrs)
        node_attrs.update(self.DEFAULT_NODE_ATTRS[(
            self.to_file,
            self.from_file is not None,
            use_dot2tex,
        )])

        return node_attrs

    def update_tf_node_attrs(self, node_attrs, update_label=False, use_dot2tex=False):
        """return dict for node_attrs based on from_file and to_file"""
        node_attrs = dict(node_attrs)
        node_attrs.update(self.DEFAULT_TF_NODE_ATTRS[use_dot2tex])
        if update_label and self.transfer_func_label:
            # the name for TF gets a nice little offset by inserting a whitespace
            # tried some other techniques, but none worked
            # if use_dot2tex: node_attrs.update(dict(texlbl=self.transfer_func_label))
            if use_dot2tex:
                label = 'TF: {}'.format(self.transfer_func_label.replace(r'\n', r'\\'))
            else:
                # only use leading whitespace for label if using dot directly
                label = ' TF: {}'.format(self.transfer_func_label)
        else:
            # label = self.transfer_func_label
            label = ''
        node_attrs.update(dict(label=label))

        return node_attrs

    def update_edge_attrs(self, edge_attrs, use_dot2tex=False):
        """return dict for edge_attrs based on transfer_func and target_coord_names"""
        edge_attrs = dict(edge_attrs)
        edge_attrs.update(self.DEFAULT_EDGE_ATTRS[(
            self.transfer_func is None,
            not self.target_coord_names,
            use_dot2tex,
        )])

        return edge_attrs

    def get_coord_labels(self, coords):
        # coords is a list of lists
        labels = []
        if self.target_coord_names:
            for coord_name in self.target_coord_names:
                for coord_group in coords:
                    if coord_name in coord_group:
                        labels.append(coord_group[0])
        return labels

class Redirecter(object):
    """The Redirecter class handles the redirection of edges,
    if only data_arrays with transfer functions shall be plotted.
    It also stores a dict with integer indices for each data array"""

    def __init__(self, use_node_index):
        # flag whether to use integer index or name (default)
        self.use_node_index = use_node_index
        # contains all redirections as {array.name: [redirected from_data_array entries]}
        self.redirect_dict = {}
        # contains all integer indices as {array.name: index}
        self.index_dict = {}
        # counter storing the current index for index_dict
        self.current_index = 1

    def add(self, key, value=None):
        """add a new key to the dicts"""
        if value is not None:
            # this means that a reference needs to be set
            values = []
            [values.extend(self.get(_)) for _ in value]
            # now all original data_arrays along from_data_array paths are in one flat list
            self.redirect_dict[key] = values
        # set index and advance counter
        self.index_dict[key] = str(self.current_index)
        self.current_index += 1

    def get(self, key):
        """get all the root from_data_array names for a given array name"""
        if key in self.redirect_dict:
            keys = []
            for _key in self.redirect_dict[key]:
                # get the root from_data_array names of the key, follow each root recursively
                keys.extend(self.get(_key))
            return keys
        else:
            # that key is not in the list, return self as a list
            return [key]

    def get_name(self, name):
        """return a name or its index based on the global flag"""
        if self.use_node_index:
            return self.index_dict[name]
        else:
            return name


class GraphvizPlotter(object):
    """main interface of graphviz plot creation for MPR data arrays"""

    # options for dot2tex
    DOT2TEX_OPTIONS = {
        'format': 'tikz',
        'texmode': 'verbatim',
        #'tikzedgelabels': True,
        'debug': True,
        #'valignmode': 'center',
        'encoding': 'latin1',
        'styleonly': True,
        #'autosize': True,
        #'prog': 'circo',
        'template': 'graphviz_dot2tex_template.tex',
    }
    # options for all nodes and edges
    DOT2TEX_NODE_EDGE_OPTIONS = {
        'node': 'shape="ellipse"',
        'edge': 'lblstyle="auto"',
    }
    # template options
    TEMPLATE_OPTIONS = {
        'graphstyle': '-latex',
        'd2tdocpreamble': dedent(r'''\
        % get x and y coordinate from tikz coordinate
        % copied from https://tex.stackexchange.com/a/33765/207160
        \makeatletter
        \newcommand{\gettikzxy}[3]{%
            \tikz@scan@one@point\pgfutil@firstofone#1\relax
            \edef#2{\the\pgf@x}%
            \edef#3{\the\pgf@y}%
        }
        \makeatother
        % produce an elliptic arc based on its center coordinate
        % based on https://tex.stackexchange.com/a/66220/207160
        \newcommand{\centerarc}[6]{
            % Syntax: \centerarc{draw options}{reference node}{initial angle}{final angle}{swell size}{label}
            \gettikzxy{(#2.center)}{\cx}{\cy}
            \gettikzxy{(#2.east)}{\ex}{\ey}
            \gettikzxy{(#2.north)}{\nx}{\ny}
            \ifthenelse { \equal {#6} {{}} } { \def\postactionarg {} }   % if #6 == blank
                                             { \def\postactionarg {[postaction={decorate,
                                    decoration={markings,
                                                raise=-10,
                                                mark=at position .5 with {\node {\contour{white}{#6}};}
                                                }
                                    }
                       ]} }   % else (not blank)
            \draw[#1] \postactionarg
                ($({\cx+(\ex-\cx+#5)*cos(#3)},{\cy+(\ny-\cy+#5)*sin(#3)})$) 
                arc[start angle=#3,end angle=#4,x radius=(\ex-\cx+#5),y radius=(\ny-\cy+#5)]
        }
        % ......................................................................
        '''),
        # optional tikzstyle option: double distance=1pt
        'd2tfigpreamble':
            '\\tikzstyle{tf_node}=[rounded corners,draw=black!20,thick,align=center]\n' +
            ''.join(['\\tikzstyle{{{}}}=[ellipse split,thick,align=center,draw={}]\n'.format(name, pen) for name, pen in zip(NODE_STYLE_NAMES, NODE_STYLE_COLORS)]) +
            ''.join(['\\tikzstyle{{{}}}=[{}]\n'.format(name, pen) for name, pen in zip(EDGE_STYLE_NAMES, EDGE_STYLE_PENS)]),
        'd2tfigpostamble': dedent(r'''\
        % now the legend ''' + LEGEND_MARKER + r'''
        \matrix[row sep=10bp,column sep=10bp] {
            % DataArrays
            \node (m11) {DataArrays}; &
            \node (m12) [da_node_from_file] {\nodepart{lower} }; &
            \node (m13) [da_node] {\footnotesize name \nodepart{lower} }; &
            \node (m14) [da_node_to_file] {\nodepart{lower} }; & \\\
            % Coordinates
            \node (m21) [align=center] {Coordinates (only shown when scaled)}; & &
            \node (da_node4)[da_node] { \nodepart{lower} $4D$}; \\\
            % coordinate labels are entered dynamically later on
            % Coordinate types
            \node (m31) {Coordinate types}; &
            \draw [dotted](-10bp,0) -- node [above=4bp,anchor=base]{\footnotesize high res.} (10bp,0); &
            \draw [dashed](-10bp,0) -- node [above=4bp,anchor=base]{\footnotesize low res.} (10bp,0); &
            \draw [solid](-10bp,0) -- node [above=4bp,anchor=base]{\footnotesize one cell} (10bp,0); & \\\
            % upscale
            \node (m41) {scaling only}; &
            \node (tf4_source) [da_node] { \nodepart{lower} }; & &
            % \centerarc{black!80,dotted}{tf4_source}{95}{175}{2}{}; & &
            \node (tf4_target) [da_node] { \nodepart{lower} };
            \centerarc{black!80,solid}{tf4_target}{95}{175}{2}{};\\\
            % broadcasting
            \node (m51) {broadcasting}; &
            \node (tf3_source) [da_node] { \nodepart{lower} $2D$}; & &
            % \centerarc{black!80,dotted}{tf3_source}{95}{175}{2}{}; & &
            \node (tf3_target) [da_node] { \nodepart{lower} $3D$};
            \centerarc{black!80,dashed}{tf3_target}{185}{265}{2}{};\\\
            % transfer
            \node (m61) {transfer only}; &
            \node (tf2_source) [da_node] { \nodepart{lower} }; &
            % \centerarc{black!80,dotted}{tf2_source}{95}{175}{2}{}; &
            \node (tf2) [tf_node] {TF}; &
            \node (tf2_target) [da_node] { \nodepart{lower} };
            \\\
            % upscale and transfer
            \node (m71) {scaling \& transfer}; &
            \node (tf1_source) [da_node] { \nodepart{lower} };
            % \centerarc{black!80,dotted}{tf1_source}{95}{175}{2}{}; 
            &
            \node (tf1) [tf_node] {TF}; &
            \node (tf1_target) [da_node] { \nodepart{lower} };
            \centerarc{black!80,solid}{tf1_target}{95}{175}{2}{};\\\
        };
        \draw [->,uo,shorten <=5, shorten >=5] (tf4_source) -- (tf4_target);
        \draw [->,uo,shorten <=5, shorten >=5] (tf3_source) -- (tf3_target);
        \draw [->,tf,shorten >=5] (tf2) -- (tf2_target);
        \draw [tf,shorten <=5] (tf2_source) -- (tf2);
        \draw [->,uo_tf,shorten >=5] (tf1) -- (tf1_target);
        \draw [uo_tf,shorten <=5,] (tf1_source) -- (tf1);
        \node (m12l) [above of=m12,yshift=-10] {\footnotesize from file};
        \node (m14l) [above of=m14,yshift=-10] {\footnotesize to file};% now the background
        \begin{scope}[on background layer]
            \node[
                draw=black,
                %fill=black!20,
                rounded corners,
                fit=(m11)(m12)(m13)(m14)(m21)(m31)(m41)(m51)(m61)(m71)(m12l)(m14l)(tf1_target)] (legend_box) {};
            %\node [above=of legend_box, yshift=-5pt] (legend_box_label) [draw=black!75, fill=white, rounded corners] {};
        \end{scope}
        '''),
    }
    # keys: use_dot2tex
    DATAARRAY_NODE_NAME_DICT = {
        False: '{}',
        True: '{} \\nodepart{{lower}} ${}$',
    }

    CENTERARC_TEMPLATE = \
        '{spacing}\\centerarc{{{formats}}}{{{ref_node}}}{{{start_angle}}}{{{end_angle}}}{{{offset}}}{{{label}}};\n'

    def __init__(self):
        # all the attributes
        self.commandLineArgs = None
        self.nml_dict = {}
        self.arrays = []
        self.graph = None
        self.redirecter = None
        self.node_attrs = {}
        self.coord_aliases = []
        self.coords = {}

    def parse_args(self):
        """parse the arguments passed from the console and display help with '-h' switch"""
        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                         description='''Visualizing aid for mpr.nml files.

            author: Robert Schweppe
            created: Mar 2019''')
        parser.add_argument('-c', '--config_file', action='store',
                            default=DEFAULT_PATH_TO_ROOT + DEFAULT_CONFIG_FILE, dest='config_file',
                            metavar='config_file',
                            help="path to config file for MPR (Default: {})".format(
                                DEFAULT_PATH_TO_ROOT + DEFAULT_CONFIG_FILE))
        parser.add_argument('-o', '--output_file', action='store',
                            default=DEFAULT_PATH_TO_ROOT + OUTPUT_FILE_BASENAME, dest='output_file',
                            help="output path (base and filename stem) for various output files (Default: {})".format(
                                DEFAULT_PATH_TO_ROOT + OUTPUT_FILE_BASENAME))
        parser.add_argument('-f', '--format', action='store',
                            default=DEFAULT_OUTPUT_FORMAT, dest='output_format', choices=FORMATS,
                            help="output format for graphviz plot (Default: {})".format(DEFAULT_OUTPUT_FORMAT))
        parser.add_argument('-e', '--engine', action='store',
                            default=DEFAULT_ENGINE, dest='engine', choices=ENGINES,
                            help="engine for rendering plot (Default: {})".format(DEFAULT_ENGINE))
        parser.add_argument('-t', '--transfer_only', action='store_true',
                            default=DEFAULT_TRANSFER_ONLY, dest='transfer_only',
                            help="plot only nodes with transfer function? (Default: {})".format(DEFAULT_TRANSFER_ONLY))
        parser.add_argument('-u', '--use_node_index', action='store_true',
                            default=DEFAULT_USE_NODE_INDEX, dest='use_node_index',
                            help="use numbers instead of names for node labelling? (Default: {})".format(
                                DEFAULT_USE_NODE_INDEX))
        parser.add_argument('-d', '--use_dot2tex', action='store_true',
                            default=DEFAULT_USE_DOT2TEX, dest='use_dot2tex',
                            help="postprocess dot output to produce tex? (Default: {})".format(
                                DEFAULT_USE_DOT2TEX))
        parser.add_argument('-g', '--group_edges', action='store_true',
                            default=DEFAULT_GROUP_EDGES, dest='group_edges',
                            help="group edges and label transfer function? (Default: {})".format(DEFAULT_GROUP_EDGES))
        parser.add_argument('-b', '--blank', action='store_true',
                            default=DEFAULT_BLANK, dest='blank',
                            help="omit any labels, make plot blank? (Default: {})".format(DEFAULT_BLANK))

        self.commandLineArgs = parser.parse_args()

        # sanity check for dot2tex option
        if self.commandLineArgs.use_dot2tex and not DOT2TEX_IS_AVAILABLE:
            warnings.warn('Module dot2tex is currently not installed, visit {} for installation guide. '.format(
                DOT2TEX_GUIDE_HTML) +
                          'Continuing without tex export.')
            self.commandLineArgs.use_dot2tex = False

    def read_and_parse_namelist(self):
        """IO for mpr.nml"""
        self._read_namelist()

        self._parse_data_arrays()

        self._parse_coordinate_aliases()

        self._parse_coordinates()

    def _parse_coordinate_aliases(self):
        # get the correct name for the key (coordinate_aliases)
        key = list(OPTIONS.keys())[0]
        self.coord_aliases = [_ for _ in self.nml_dict[key] if not all([item is None for item in _])]

    def _parse_coordinates(self):
        # get the correct name for the key (coordinate_aliases)
        names = self.nml_dict[list(OPTIONS.keys())[5]]
        counts = self._get_option(list(OPTIONS.keys())[12], default=None, target_length=len(names))
        steps = self._get_option(list(OPTIONS.keys())[11], default=None, target_length=len(names))
        values = self._get_option(list(OPTIONS.keys())[8], default=[None], target_length=len(names))
        for iname, name in enumerate(names):
            if name is None:
                continue
            if counts[iname] == 1:
                # "integral" coordinate with one cell
                self.coords[name] = 'solid'
            elif counts[iname] or steps[iname] or not all([_ is None for _ in values[iname]]):
                # presumably low-res
                self.coords[name] = 'dashed'
            else:
                # not specified, presumably read from file dynamically, presumably high-res
                self.coords[name] = 'dotted'

    def _parse_data_arrays(self):
        # gather all the relevant properties occurring in namelist
        all_props_keys = [
            list(OPTIONS.keys())[20],
            list(OPTIONS.keys())[21],
            list(OPTIONS.keys())[22],
            list(OPTIONS.keys())[23],
            list(OPTIONS.keys())[24],
            list(OPTIONS.keys())[27],
            list(OPTIONS.keys())[35],
        ]
        # not all of the properties need to appear in namelist file, so we check which exist and broadcast
        section = self.nml_dict[all_props_keys[0][0]]
        da_keys, da_props = zip(*[(key[1], self.nml_dict[key]) for key in all_props_keys if key[1] in section])
        # loop over all entries in namelist and create DataArray objects and set to self
        # not all properties are set, so we want default values (None) to be filled in unset properties
        for da_prop in zip_longest(*da_props):
            # names do not need to be set consecutively in namelists, we thus need to ignore gaps
            if da_prop[0] is not None:
                # init DataArray by zipping information and key and append to global list
                self.arrays.append(DataArray(**{k: v for k, v in zip(da_keys, da_prop)}))

    def _read_namelist(self):
        # read namelists into dict, use start_index for implicitly sized arrays
        nml_parser = f90nml.Parser()
        nml_parser.global_start_index = 1
        self.nml_dict = nml_parser.read(self.commandLineArgs.config_file)

    def _format_dot2tex_arg(self, key, value):
        if value is True:
            return '--{}'.format(key)
        elif value is False:
            raise Exception
        else:
            return '--{} {}'.format(key, value)

    def init_graph(self):
        """init the graphviz.Digraph instance and set some default properties"""
        # initiate Graph with settings
        self.graph = Digraph(
            name=self.nml_dict[list(OPTIONS.keys())[1]],
            filename=self.commandLineArgs.output_file,
            format=self.commandLineArgs.output_format,
            engine=self.commandLineArgs.engine,
        )

        # setting the font works only that way
        if self.commandLineArgs.use_dot2tex:
            for key, value in self.TEMPLATE_OPTIONS.items():
                self.graph.body.append('\t{} = "{}"'.format(key, value))
            new_string = '\td2toptions = "{}"'.format(' '.join([self._format_dot2tex_arg(key, value) for key, value in self.DOT2TEX_OPTIONS.items() if value is not False]))
            self.graph.body.append(new_string)
            for key, value in self.DOT2TEX_NODE_EDGE_OPTIONS.items():
                self.graph.body.append('\t{} [{}]'.format(key, value))
        else:
            self.graph.attr('graph', fontname=FONT)
            self.graph.attr('node', fontname=FONT, penwidth='3')
            self.graph.attr('edge', fontname=FONT)

        # this is something standard, escpecially for neato
        if self.commandLineArgs.engine == 'neato':
            self.graph.attr(ordering='out')
        self.graph.attr(overlap='false')

        # this is prepared for node initialization
        self.redirecter = Redirecter(self.commandLineArgs.use_node_index)
        # node label is intentionally left empty
        if self.commandLineArgs.blank:
            self.node_attrs['label'] = ''
        # else: label is taken from node name by default (assigned dynamically)

    def build_graph(self):
        """iteratively loop over all DataArrays and set the nodes, prepare edges and add edges"""
        # loop over all arrays
        for i_array, array in enumerate(self.arrays):
            # in case of plotting only DataArrays resulting from transfer functions (e.g. only upscaling is performed):
            # do not plot the node, but register the redirection of from_data_arrays of all child nodes
            if self.commandLineArgs.transfer_only and \
                    array.transfer_func_label is None and \
                    array.from_data_arrays:
                # add it to the redirecter instance
                self.redirecter.add(array.name, array.from_data_arrays)
                continue
            else:
                # register only the index
                self.redirecter.add(array.name)
            # add the node and set attributes accordingly
            self.node_attrs = array.update_node_attrs(
                self.node_attrs,
                use_dot2tex=self.commandLineArgs.use_dot2tex
            )
            self.graph.node(self.redirecter.get_name(array.name),
                            **self.node_attrs)

            # prepare the edges for the node (name and properties)
            edge_name, edge_attrs = self._prepare_edges(array)

            # now draw the edges
            self._add_edges(array, edge_name, edge_attrs)

    def _prepare_edges(self, array):
        """prepare the edges from the from_data_array arrays to array"""
        edge_attrs = array.update_edge_attrs(
            dict(),
            use_dot2tex=self.commandLineArgs.use_dot2tex
        )

        # do we want to group all incoming edges into one edge?
        #if len(array.from_data_arrays) > 0 and self.commandLineArgs.group_edges:
        if array.transfer_func is not None:
            # create invisible temporary node with unique name (prefix underscore)
            temp_node_name = '_{}'.format(self.redirecter.get_name(array.name))
            # add the dummy node and set attributes accordingly
            node_attrs = array.update_tf_node_attrs(
                self.node_attrs,
                update_label=not self.commandLineArgs.blank,
                use_dot2tex=self.commandLineArgs.use_dot2tex
            )
            #merged_attrs = node_attrs
            #merged_attrs['style'] = edge_attrs['style']

            self.graph.node(temp_node_name, **node_attrs)

            # update general edge attrs
            # apply blank label also here for labeling edges
            temp_edge_attrs = array.update_edge_attrs(
                dict(),
                use_dot2tex=self.commandLineArgs.use_dot2tex
            )

            # draw edge from TF node to target node now!
            # later on we only loop over all predictors, so it is necessary to do it now
            self.graph.edge(
                temp_node_name,
                self.redirecter.get_name(array.name),
                **temp_edge_attrs
            )
            # the edges leading to the invisible node do not get a direction arrow
            edge_attrs.update(dict(arrowhead='none'))
        else:
            # default: draw a direction arrow directly to the target with some minimum length
            temp_node_name = self.redirecter.get_name(array.name)
            edge_attrs.update(dict(
                minlen='1.5'
            ))
        return temp_node_name, edge_attrs

    def _add_edges(self, array, target_edge, edge_attrs):
        """add the edges to the plot"""
        prev_src = ''
        # we need to iterate over all from_data_arrays (edge from_data_array -> array)
        for src in array.from_data_arrays:
            # there might be some redirection going on for each from_data_array, we need to cover each
            for red_src in self.redirecter.get(src):
                redirected_src = self.redirecter.get_name(red_src)
                # connect each node only once (duplicate connection may occur due to redirections)
                if redirected_src != prev_src:
                    # add each edge to its predictor
                    self.graph.edge(redirected_src, target_edge, **edge_attrs)
                    prev_src = redirected_src

    def render_graph(self):
        """Save to source, target format and optional export formats (tex)"""

        output_file_dotx = pl.Path('{}.{}'.format(self.commandLineArgs.output_file, 'gv'))
        output_file_final = output_file_dotx.with_suffix('.{}'.format(self.commandLineArgs.output_format))
        print('Saving to {}'.format(output_file_final))
        self.graph.render(output_file_dotx)

        # directly call dot2tex on dotx file
        if self.commandLineArgs.use_dot2tex:
            # open the output dotx file and read stream
            with open(output_file_dotx) as dotx_file:
                dotx = dotx_file.read()

            # apply conversion
            texcode = d2t.dot2tex(dotx,
                                  **self.DOT2TEX_OPTIONS,
                                  #**self.DOT2TEX_D2TNODEOPTIONS_DICT,
                                  #**self.DOT2TEX_D2TEDGEOPTIONS_DICT,
                                  #**self.DOT2TEX_D2TGRAPHSTYLE_DICT,
                                  )

            # save to .tex file
            texcode = self.insert_coordinate_notations(texcode)
            output_file_tex = pl.Path('{}.{}'.format(self.commandLineArgs.output_file, 'tex'))
            print('Exporting to {}'.format(output_file_tex))
            with open(output_file_tex, 'w') as tex_file:
                tex_file.write(texcode)

            if self.DOT2TEX_OPTIONS.get('debug', False):
                logstream = d2t.get_logstream()
                if logstream is not None:
                    print(logstream.getvalue())

    def insert_coordinate_notations(self, texcode):
        splitcode = texcode.splitlines(keepends=True)
        addons = ''
        marker_found = False
        for iline, line in enumerate(splitcode):
            if LEGEND_MARKER in line:
                marker_found = True
            if '[da_node' in line and not marker_found:
                line, addon = self._modify_nodes(line=line)
                addons += addon
            if '\\draw' in line and not marker_found:
                line = self._modify_paths(line)
            if '\end{tikzpicture}' in line:
                # enter the coordinate label for the legend
                # da_node 4, default coordinate
                _, addon = self._modify_nodes(coord_names=[aliases[0] for aliases in self.coord_aliases],
                                              used_labels=[aliases[0] for aliases in self.coord_aliases],
                                              ref_node='da_node4')
                addons += addon
                line = addons + line
            splitcode[iline] = line
        return ''.join(splitcode)

    def _modify_paths(self, line):
        # extract all node names
        nodes = re.findall(r'\((\w+)\)', line)
        array_names = [arr.name for arr in self.arrays]
        style_addons = []
        # check for each node on the path, if it refers to a node modified by our script
        for inode, node in enumerate(nodes):
            if node in array_names:
                if inode == 0:
                    style_addons.append('shorten <={}'.format(COORD_STYLE_OFFSET))
                else:
                    style_addons.append('shorten >={}'.format(COORD_STYLE_OFFSET))

        if style_addons:
            insert_index = re.search(r'\[([\S]+)\]', line).end() - 1
            line = ','.join([line[:insert_index], *style_addons]) + line[insert_index:]
        return line

    def _modify_nodes(self, line='', used_labels=None, coord_names=None, ref_node=None):
        addons = ''
        if line:
            # extract node name
            ref_node = re.findall(r'\((\w+)\)', line)[0]
            # get the DataArray with the same name
            arr = self.get_array_by_name(ref_node)
            addon_format = NODE_COLORS_DICT[arr.update_node_attrs({}, True)['style']]
            # generate coordinate label
            used_labels, all_labels = self.get_coord_labels(arr)
            # get current label
            current_label = re.findall(r'\{([\w\\]+)\}', line)[0]
            # insert new combined label
            lower_node_part = '?'
            if all_labels:
                lower_node_part = '{}D'.format(len(all_labels))
            line = line.replace(current_label, self.DATAARRAY_NODE_NAME_DICT[True].format(current_label, lower_node_part))
            coord_names = arr.target_coord_names
        else:
            if used_labels is None:
                used_labels = []
            if coord_names is None:
                coord_names = []
            addon_format = NODE_COLORS_DICT[NODE_STYLE_NAMES[3]]
            all_labels = used_labels

        if not coord_names:
            addons = self.CENTERARC_TEMPLATE.format(
                spacing=' ' * 2,
                formats=addon_format + ',dotted',
                ref_node=ref_node,
                start_angle=95,
                end_angle=445,
                offset=COORD_STYLE_OFFSET,
                label='',
            )
        else:
            for name, label in zip(coord_names, all_labels):
                if label not in used_labels:
                    continue
                # start and end angle of arc
                start_angle, end_angle = self.get_angles_coord_visuals(name)
                addons += self.CENTERARC_TEMPLATE.format(
                    spacing=' ' * 2,
                    formats=addon_format + ',{}'.format(self.coords.get(name, 'dotted')),
                    ref_node=ref_node,
                    start_angle=start_angle,
                    end_angle=end_angle,
                    offset=COORD_STYLE_OFFSET,
                    label='${}$'.format(label)
                )
        return line, addons

    def get_array_by_name(self, name):
        for arr in self.arrays:
            if arr.name == name:
                return arr
        return None

    def get_angles_coord_visuals(self, coord_name):
        for i_alias, alias in enumerate(self.coord_aliases):
            if coord_name in alias:
                break
        return (int(i_alias * 360 / len(self.coord_aliases) + 5) + 90), \
               (int((i_alias + 1) * 360 / len(self.coord_aliases) - 5) + 90)

    def _get_option(self, keys, default, target_length):
        d = self.nml_dict
        for key in keys:
            d = d.get(key, [default])
            if d is None:
                break
        return d + [default] * (target_length - len(d))

    def get_coord_labels(self, arr):
        all_labels = arr.get_coord_labels(self.coord_aliases)
        used_labels = all_labels
        pred = arr
        # infer the labels from its parent dataarrays recursively
        while len(all_labels) == 0:
            if not pred.from_data_arrays:
                break
            # all from_data_arrays must have same shape in a valid mpr.nml, so it is okay to select the first only
            pred = self.get_array_by_name(pred.from_data_arrays[0])
            all_labels = pred.get_coord_labels(self.coord_aliases)
            if len(all_labels) > 0:
                arr.target_coord_names = pred.target_coord_names

        if arr.from_data_arrays:
            pred = self.get_array_by_name(arr.from_data_arrays[0])
            used_labels = [label for i, label in enumerate(all_labels) if arr.target_coord_names[i] not in pred.target_coord_names]

        return used_labels, all_labels


# SCRIPT
if __name__ == '__main__':
    plot = GraphvizPlotter()
    # parse the arguments
    plot.parse_args()
    # read the namelist
    plot.read_and_parse_namelist()
    # init the graphviz Digraph
    plot.init_graph()
    # build all nodes and edges
    plot.build_graph()
    # render the graph
    plot.render_graph()
