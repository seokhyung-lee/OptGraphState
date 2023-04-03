import time
import os
import sys

import networkx as nx

from .graph_tools import *
from .plotting import *

try:
    import parmap
except ModuleNotFoundError:
    pass


def _max_seed():
    return 2**32


class OptGraphState:
    graph: ig.Graph
    graph_info: set
    unraveled_graph: ig.Graph
    fusion_network: ig.Graph
    unraveled_bcss: set
    unraveled_cliques: set
    data: dict

    def __init__(self,
                 graph=None,
                 edges=None,
                 shape=None,
                 prms=None,
                 cliffords=None,
                 unraveled_graph=None,
                 fusion_network=None):
        """
        Class for calculating and optimizing the resource overhead of the fusion-based generation of a graph
        state.

        The graph of the concerned graph state can be given by the following three ways:
        1. Given explicitly as igraph.Graph or networkx.Graph.
        2. Given by a list of edges.
        3. Chosen among predefined graphs.

        An unraveled graph and a fusion network of the graph can be explicitly given by the parameters
        'unraveled_graph' and 'fusion_network', respectively, which should be igraph.Graph.

        Parameters
        ----------
        graph : None or igraph.Graph or networkx.Graph
            Graph of the concerned graph. If it is given, 'edges', 'shape', and 'prms' are ignored.
            If it is networkx.Graph, it is internally converted to igraph.Graph.

        edges : None or list of 2-tuples of integers
            List of edges that form the concerned graph. Each integer in the tuples indicates a vertex label.
            If it is given and 'graph' is None, 'shape' and 'prms' are ignored.

        shape : None or str.
            Shape of the concerned graph chosen among predefined graphs.
            One of [None, 'random', 'complete', 'star', 'linear', 'cycle', 'lattice', 'tree', 'rhg',
            'repeater', 'parity_encoding', 'ptqc'].

            > shape='random' : Random graph for fixed numbers of vertices and edges, sampled by the Erdos-Renyi
            model.
                - prms[0] (int) : Number of vertices.
                - prms[1] (int) : Number of edges.
            > shape='complete','star','linear','cycle' : Complete, star, linear, and cycle graph, respectively.
                - prms[0] (int) : Number of vertices.
            > shape='lattice' : 2D lattice graph.
                - prms (2-tuple ints) : Numbers of repeated vertices along the two axes.
            > shape='tree' : Tree graph where all branches in each generation has an equal number of children.
                - prms[0] (int) : Degree of the root vertex.
                - prms[i] (int, i >= 1) : Number of the children of each ith-generation branch.
            > shape='rhg' : Raussendorf-Harrington-Goyal lattice where all the three bondaries are primal.
                - prms (3-tuple of ints) : Size of the lattice along the three axes in the unit of a cell.
            > shape='repeater' : Repeater graph with 4m vertices.
                - prms[0] (int) : Parameter 'm'.
            > shape='parity_encoding' : (n, m) parity-encoded graph.
                - prms[0] (ig.Graph) : Logical-level graph.
                - prms[1], prms[2] (int) : Parameters 'n' and 'm' of the parity encoding.
            > shape='ptqc' : Microcluster for parity-encoding-based topological quantum computing protocol.
                - prms[0], prms[1] (int) : Parameter 'n' and 'm' of the parity encoding.
                - prms[2] (bool) : H-configuration.
                - prms[3] (bool) : Whether the microcluster is central (True) or side (False) one.

        prms : None or tuple of (ints or ig.Graph or bool)
            Parameters for a predefined graph. See the description for 'shape'.

        cliffords : None or list of str
            Local clifford operations applied on the qubits of the graph state.

        unraveled_graph : None or igraph.Graph or networkx.Graph
            Pregiven unraveled graph. The code does not check the validity of the given unraveled graph.
            If it is networkx.Graph, it is internally converted to igraph.Graph.

        fusion_network : None or igraph.Graph or networkx.Graph
            Pregiven fusion network. The code does not check the validity of the given fusion network.
            If it is networkx.Graph, it is internally converted to igraph.Graph.
        """

        def convert_type(g, varname):
            if g is None:
                return None
            elif isinstance(g, ig.Graph):
                return g
            elif isinstance(g, nx.Graph):
                return ig.Graph.from_networkx(g)
            else:
                raise TypeError(f'Parameter {varname} should be igraph.Graph or networkx.Graph.')

        if graph is not None:
            self.graph = convert_type(graph, 'graph')
            self.graph_info = None

        elif edges is not None:
            self.graph = ig.Graph(edges=edges)
            self.graph_info = None

        elif shape is not None:
            self.graph, self.graph_info = get_sample_graph(shape, *prms)

        else:
            raise ValueError("At least one of graph, edges, and shape should be given.")

        if self.graph_info is None:
            self.graph_info = {'num_vertices': self.graph.vcount(), 'num_edges': self.graph.ecount()}

        self.graph.vs['name'] = [str(vid) for vid in range(self.graph.vcount())]

        if cliffords is not None and edges is None:
            self.graph.vs['clifford'] = cliffords

        self.unraveled_graph = convert_type(unraveled_graph, 'unraveled_graph')
        self.unraveled_bcss = set()
        self.unraveled_cliques = set()

        self.fusion_network = convert_type(fusion_network, 'fusion_network')

        self.data = {}

    def initialize(self):
        """
        Initialize the created unraveled graph and fusion network.
        """
        self.unraveled_graph = None
        self.fusion_network = None
        self.unraveled_bcss = set()
        self.unraveled_cliques = set()
        self.data = {}

    def get_vertex_data(self, vertex):
        """
        Get the data of a vertex of the unraveled graph or a node of the fusion network.

        For a given vertex name, it returns the attribute data of the vertex or node with the same name in the
        unraveled graph and fusion network.

        Parameters
        ----------
        vertex : str or int
            Name of the vertex.

        Return
        ------
        data : dict
            Data of the vertex or node in the unraveled graph and fusion network.
        """
        name = vertex if isinstance(vertex, str) else str(vertex)
        data = {}

        try:
            qubit_data = self.unraveled_graph.vs.find(name=name).attributes()
            del qubit_data['name']
            data['unraveled_graph'] = qubit_data
        except (ValueError, AttributeError):
            pass

        try:
            fusion_data = self.fusion_network.vs.find(name=name).attributes()
            del fusion_data['name']
            data['fusion_network'] = fusion_data
        except (ValueError, AttributeError):
            pass

        return data

    def unravel_graph(self, unravel_bcs_first='random', plot=False, verbose=False):
        """
        Unravel bipartitely-complete subgraphs (BCSs) and cliques of the graph.

        The unraveled graph is saved in self.unraveled_graph as igraph.Graph.

        Parameters
        ----------
        unravel_bcs_first : one of [True, False, 'random'] (default: 'random')
            If it is True, BCSs are unraveled first.
            If it is False, cliques are unraveled first.
            If it is 'random', the order is randomly chosen.
        plot : bool (default: False)
            Whether to plot the unraveled graph after unraveling.
        verbose : bool (default: False)

        """
        if unravel_bcs_first == 'random':
            unravel_bcs_first = np.random.choice([True, False])

        if unravel_bcs_first:
            bcss = self.unravel_bipartitely_complete_subgraphs(verbose=verbose)
            cliques = self.unravel_cliques(verbose=verbose)
        else:
            cliques = self.unravel_cliques(verbose=verbose)
            bcss = self.unravel_bipartitely_complete_subgraphs(verbose=verbose)

        if plot or verbose:
            if verbose:
                print('[Final]')
            self.plot_graph(unraveled=True)
            plt.show()

        self.data['unravel'] = True
        self.data['unravel_bcs_first'] = unravel_bcs_first

        return bcss, cliques

    # @logging_time
    def unravel_bipartitely_complete_subgraphs(self, verbose=False):
        if self.unraveled_graph is None:
            self.unraveled_graph = self.graph.copy()

        graph = self.unraveled_graph

        vs_attrs = graph.vs.attributes()
        if 'ext_fusion' not in vs_attrs:
            graph.vs['ext_fusion'] = None
        if 'clifford' not in vs_attrs:
            graph.vs['clifford'] = None

        unraveled_bcss = {}

        new_vertex_name = graph.vcount()

        stage = 1
        while True:
            # Repeat until there are no bipartitely-complete subgraphs
            bcs_exist = False
            while True:
                bcss = find_all_nonoverlapping_bcs(graph, get_name=True)
                if bcss:
                    bcs_exist = True
                else:
                    break

                unraveled_bcss[stage] = bcss
                eids_to_remove = []
                for part1, part2 in bcss:
                    if verbose:
                        graph.delete_edges(eids_to_remove)
                        eids_to_remove.clear()
                        print('bcs to unravel =', part1, '&', part2)
                        vertex_color = []
                        for v in graph.vs:
                            if v['name'] in part1:
                                vertex_color.append('orange')
                            elif v['name'] in part2:
                                vertex_color.append('blue')
                            else:
                                vertex_color.append('white')
                        self.plot_graph(unraveled=True,
                                        vertex_color=vertex_color)
                        plt.show()

                    eids_to_remove.extend([graph.get_eid(vname1, vname2) for vname1, vname2
                                           in itertools.product(part1, part2)])

                    vname1 = str(new_vertex_name)
                    vname2 = str(new_vertex_name + 1)
                    new_v1 = graph.add_vertex(name=vname1, ext_fusion=vname2, clifford=None)
                    new_v2 = graph.add_vertex(name=vname2, ext_fusion=vname1, clifford=None)
                    new_vertex_name += 2
                    # ext_fusions.append((new_v1['name'], new_v2['name'], bcs))

                    graph.add_edges(itertools.product([new_v1], part1))
                    graph.add_edges(itertools.product([new_v2], part2))

                graph.delete_edges(eids_to_remove)

                stage += 1

            if not bcs_exist:
                break

            # if stage == 2:
            #     continue

            # Rejoin unnecessary external fusions
            # np.random.shuffle(ext_fusions)
            # rejoined_bcss = []
            # for i_ext_fusion, ext_fusion in enumerate(ext_fusions):
            #     vname1, vname2, bcs = ext_fusion
            #
            #     try:
            #         if graph.degree(vname1) == 1 or graph.degree(vname2) == 1:
            #             # vids_to_remove.extend([vid1, vid2])
            #
            #             if verbose:
            #                 print(f'rejoined fusion = ({vname1}, {vname2})')
            #                 self.plot_graph(unraveled=True,
            #                                 vertices_to_highlight=[vname1, vname2],
            #                                 vertex_color_highlight='orange')
            #                 plt.show()
            #
            #             new_edges = itertools.product(graph.neighbors(vname1), graph.neighbors(vname2))
            #             rejoined_bcss.append(bcs)
            #             graph.add_edges(new_edges)
            #             graph.delete_vertices([vname1, vname2])
            #     except ValueError:
            #         pass
            #
            # if rejoined_bcss:
            #     unraveled_bcss[f'rejoined_{stage}'] = rejoined_bcss
            #
            # else:
            #     break

        self.unraveled_bcss = unraveled_bcss

        return unraveled_bcss

    # @logging_time
    def unravel_cliques(self, verbose=False):
        if self.unraveled_graph is None:
            self.unraveled_graph = self.graph.copy()

        graph = self.unraveled_graph

        vs_attrs = graph.vs.attributes()
        if 'ext_fusion' not in vs_attrs:
            graph.vs['ext_fusion'] = None
        if 'clifford' not in vs_attrs:
            graph.vs['clifford'] = None

        unraveled_cliques = {}
        stage = 1

        def apply_clifford(v, clifford):
            org_clifford = v['clifford']
            if org_clifford is None:
                new_clifford = clifford
            else:
                new_clifford = '-'.join([clifford, org_clifford])
            v['clifford'] = new_clifford

        while True:
            cliques = find_all_nonoverlapping_cliques(graph, get_name=True)

            if not cliques:
                break

            unraveled_cliques[stage] = cliques

            for clique in cliques:
                if verbose:
                    print('clique to unravel =', clique)
                    vertex_color = ['orange' if vname in clique else 'white' for vname in graph.vs['name']]
                    self.plot_graph(unraveled=True, vertex_color=vertex_color)
                    plt.show()
                clique = list(clique)
                clique_size = len(clique)

                # Choose a vertex to apply LC
                degrees = graph.degree(clique)
                min_degree = min(degrees)
                if min_degree == clique_size - 1:
                    # There exists a vertex in the clique that doesn't have outer edges
                    vname_LC = [vname for vname, deg in zip(clique, degrees) if deg == min_degree]
                    need_to_separate = False
                else:
                    vname_LC = clique
                    need_to_separate = True

                if len(vname_LC) > 1:
                    vname_LC = np.random.choice(vname_LC)
                else:
                    vname_LC = vname_LC[0]
                v_LC = graph.vs.find(name=vname_LC)

                # Separate the edges (E) incident to v_LC outside the clique from v_LC
                eids_to_delete = []
                if need_to_separate:
                    # Vertex connected with E
                    new_v1 = graph.add_vertex(name=str(graph.vcount()),
                                              clifford=v_LC['clifford'],
                                              ext_fusion=v_LC['ext_fusion'])

                    # Vertex having an external fusion with v_LC
                    new_v2 = graph.add_vertex(name=str(graph.vcount()),
                                              clifford=None,
                                              ext_fusion=v_LC['name'])

                    graph.add_edge(new_v1, new_v2)

                    ngh_vids = graph.neighbors(vname_LC)
                    for ngh_vid in ngh_vids:
                        if graph.vs[ngh_vid]['name'] not in clique:
                            graph.add_edge(ngh_vid, new_v1)
                            eids_to_delete.append(graph.get_eid(vname_LC, ngh_vid))

                    vname_org_ext_fusion = v_LC['ext_fusion']
                    if vname_org_ext_fusion is not None:
                        v_org_ext_fusion = graph.vs.find(name=vname_org_ext_fusion)
                        v_org_ext_fusion['ext_fusion'] = new_v1['name']

                    v_LC['ext_fusion'] = new_v2['name']

                # Apply LC
                adj_vnames = set(clique) - {vname_LC}
                apply_clifford(v_LC, 'R_X')
                for adj_vname in adj_vnames:
                    adj_v = graph.vs.find(name=adj_vname)
                    apply_clifford(adj_v, 'R_Z')

                new_eids_to_delete = [graph.get_eid(vname1, vname2) for vname1, vname2
                                      in itertools.combinations(adj_vnames, r=2)]
                eids_to_delete.extend(new_eids_to_delete)
                graph.delete_edges(eids_to_delete)

            stage += 1

        self.unraveled_cliques = unraveled_cliques

        return unraveled_cliques

    # @logging_time
    def build_fusion_network(self,
                             use_unraveled_graph=True,
                             plot=False,
                             verbose=False
                             ):
        graph = self.unraveled_graph if use_unraveled_graph else self.graph
        if graph is None:
            raise ValueError("No unraveled qubit graph created.")

        # Fusion network
        network = ig.Graph()
        self.fusion_network = network
        nodes_corr = {}

        # Links inside star graphs
        for v in graph.vs:
            num_internal_nodes = v.degree() - 1
            vname_xrot = v['name']
            if num_internal_nodes >= 1:

                # Add nodes
                nid_init = network.vcount()
                seed = np.random.randint(0, num_internal_nodes)
                node_names = [vname_xrot if i == 0 else f'{vname_xrot}-{i}'
                              for i in itertools.chain(range(seed, -1, -1), range(seed + 1, num_internal_nodes))]
                attr = {
                    'name': node_names,
                    'seed': [True if i == seed else False for i in range(num_internal_nodes)],
                    'clifford_root': None,  # Clifford gate applied on the root qubit
                    'clifford_leaves': None  # Clifford gate applied on the leaf qubits
                }
                network.add_vertices(num_internal_nodes, attributes=attr)
                nodes_corr[vname_xrot] = network.vs[nid_init:nid_init + num_internal_nodes]

                if num_internal_nodes >= 2:
                    # Connect internal links
                    links = [(nid, nid + 1) for nid in range(nid_init, nid_init + num_internal_nodes - 1)]
                    attr = {
                        # One of RR, RL, and LL, which respectively means that a fusion is performed on two roots,
                        # one root and one leaf, and two leaves of GHZ-3 states.
                        'kind': "RL",
                        # Root node name, if kind == 'RL'
                        'root_node': [node_names[i + 1] if i < seed else node_names[i]
                                      for i in range(num_internal_nodes - 1)],
                        # Whether the fusion is external or internal
                        # 'external': False,
                    }
                    network.add_edges(links, attributes=attr)

        # Links between star graphs
        for e in graph.es:
            vs = [e.source_vertex, e.target_vertex]
            deg_vs = graph.degree(vs)

            if deg_vs[0] > 1 and deg_vs[1] > 1:
                nodes_to_connect = []
                for v in vs:
                    vname = v['name']
                    nodes = nodes_corr[vname].select(lambda node: node.degree() < (2 if node['seed'] else 3))
                    nodes_to_connect.append(np.random.choice(nodes))
                network.add_edge(*nodes_to_connect,
                                 kind='LL',
                                 root_node=None,
                                 # external=False
                                 )

        if verbose:
            print("Fusion network of the unraveled graph:")
            self.plot_fusion_network()
            plt.show()

        # Set of seed node names where root vertices are connected by external fusions.
        root_connected = set()
        is_node_not_full \
            = lambda node: node.degree() < (2 if node['seed'] and (node['name'] not in root_connected) else 3)

        def get_nodes_containing_vertex(vname):
            try:
                node = nodes_corr[vname].find(name=vname)
                root = True
            except KeyError:
                root_name = graph.vs[graph.neighbors(vname)[0]]['name']
                nodes = nodes_corr[root_name].select(is_node_not_full)
                node = np.random.choice(nodes)
                root = False

            return node, root

        # Clifford operations
        vs_attrs = graph.vs.attributes()
        if 'clifford' in vs_attrs:
            for v_cl in graph.vs.select(clifford_ne=None):
                vname = v_cl['name']
                node, root = get_nodes_containing_vertex(vname)
                cl = v_cl['clifford']

                if root:
                    node['clifford_root'] = cl

                else:
                    key = v_cl['ext_fusion']
                    if key is None:
                        key = f'final_{v_cl["name"]}'

                    clifford_leaves = node['clifford_leaves']
                    try:
                        clifford_leaves[key] = cl
                    except TypeError:
                        node['clifford_leaves'] = {key: cl}

        if verbose:
            print("Apply Clifford operations:")
            self.plot_fusion_network()
            plt.show()

        # Links by external fusions
        if 'ext_fusion' in vs_attrs:
            done = set()
            for v1 in graph.vs.select(ext_fusion_ne=None):

                vname1 = v1['name']
                if vname1 in done:
                    continue
                vname2 = v1['ext_fusion']
                done.add(vname2)

                node1, root1 = get_nodes_containing_vertex(vname1)
                node2, root2 = get_nodes_containing_vertex(vname2)
                if root1 and root2:
                    kind = 'RR'
                    root_node = None
                    root_connected.add(node1['name'])
                    root_connected.add(node2['name'])
                elif not root1 and not root2:
                    kind = 'LL'
                    root_node = None
                else:
                    kind = 'RL'
                    root_node = node1['name'] if root1 else node2['name']
                    root_connected.add(root_node)
                network.add_edge(node1, node2, kind=kind, root_node=root_node)

        network.es['name'] = [str(eid) for eid in range(network.ecount())]

        self.fusion_network = network

        if plot or verbose:
            print("Final:")
            self.plot_fusion_network()
            plt.show()

    def contract_edge(self, fusion_network, ename_to_merge):
        ename_to_merge = str(ename_to_merge)

        e_to_merge = fusion_network.es.find(name=ename_to_merge)
        v_merged, v_removed = e_to_merge.source_vertex, e_to_merge.target_vertex
        enames_updated_weight = []

        if v_merged.degree() < v_removed.degree():
            v_merged, v_removed = v_removed, v_merged

        vname_merged = v_merged['name']
        vname_removed = v_removed['name']

        v_merged['weight'] = e_to_merge['weight']
        v_merged['step'] = max(v_merged['step'], v_removed['step']) + 1

        # if vname_merged == vname_removed:
        #     # eids_to_delete = [e_to_merge.index]
        #
        # else:
        assert vname_merged != vname_removed

        eids_to_delete = list(set(fusion_network.incident(v_removed)))
        v_removed['on'] = False

        for eid_connected in eids_to_delete:
            e_connected = fusion_network.es[eid_connected]
            ename_connected = e_connected['name']
            if ename_connected != ename_to_merge:
                enames_updated_weight.append(ename_connected)
                vs_ngh = e_connected.source_vertex, e_connected.target_vertex
                new_edge = [v_merged if v_ngh['name'] == vname_removed else v_ngh for v_ngh in vs_ngh]
                if new_edge[0] == new_edge[1]:  # If a loop is formed
                    v_vrt = fusion_network.add_vertex(name=f'vrt_{fusion_network.vcount()}', weight=0, step=0)
                    new_edge = [v_merged, v_vrt]

                fusion_network.add_edge(*new_edge,
                                        name=ename_connected,
                                        weight=None,
                                        )

        fusion_network.delete_edges(eids_to_delete)

        p_succ = fusion_network['p_succ']
        for eid_connected in list(set(fusion_network.incident(v_merged))):
            e_connected = fusion_network.es[eid_connected]
            v_ngh1, v_ngh2 = e_connected.source_vertex, e_connected.target_vertex
            # if v_ngh1.index == v_ngh2.index:
            #     new_weight = v_ngh1['weight'] / p_succ
            # else:
            assert v_ngh1 != v_ngh2
            new_weight = (v_ngh1['weight'] + v_ngh2['weight']) / p_succ
            e_connected['weight'] = new_weight

        self.fusion_network.es.find(name=ename_to_merge)['step'] = v_merged['step']

        return enames_updated_weight

    # @logging_time
    def calculate_overhead(self,
                           p_succ=0.5,
                           strategy='weight_and_matching',
                           fusion_order=None,
                           get_fusion_order=False,
                           # use_gt=False
                           ):
        """
        Sample one generating scheme for a graph state and calculate the correponding resource overhead.

        Parameters
        ----------
        eta :
        strategy :

        Returns
        -------
        """

        if self.fusion_network is None:
            raise ValueError("No fusion network created")

        # if use_gt and not is_graph_tool_imported():
        #     raise ModuleNotFoundError('Module graph-tool is not installed.')

        # Trivial cases
        node_num = self.fusion_network.vcount()
        if node_num == 0:
            self.data['overhead'] = 0
            self.data['step'] = 0
            return self.data
        elif node_num == 1:
            self.data['overhead'] = 1
            self.data['step'] = 0
            return self.data

        if fusion_order is None:
            fusion_order = []
            is_fusion_order_given = False
        else:
            is_fusion_order_given = True

        self.fusion_network.es['step'] = None

        # Initialize intermediate fusion network
        network = self.fusion_network.copy()
        network['p_succ'] = p_succ
        network.vs['weight'] = 1
        network.vs['step'] = 0
        network.vs['on'] = True
        network.es['weight'] = 2 / p_succ
        del network.es['step']

        turn = 0

        # t_converting = 0
        # t_matching = 0
        # t_postprocess = 0

        # Iterate until no edges remain in the fusion network
        while True:
            if not network.ecount():
                break

            if is_fusion_order_given:
                enames_curr_step = [str(fusion_order[turn])]
                is_parellel = True

            elif strategy == 'weight':
                min_weight = min(network.es['weight'])
                eids_min_weight = network.es.select(weight=min_weight)
                enames_curr_step = eids_min_weight['name']
                is_parellel = len(enames_curr_step) == 1

            elif strategy == 'betweenness':
                eb = np.array(network.edge_betweenness())
                min_eb = np.min(eb)
                eids_curr_step = np.nonzero(eb == min_eb)[0]
                enames_curr_step = [network.es[eid]['name'] for eid in eids_curr_step]
                is_parellel = len(enames_curr_step) == 1

            elif strategy == 'weight_and_betweenness':
                min_weight = min(network.es['weight'])
                eids_min_weight = network.es.select(weight=min_weight)
                es_min_ovh = eids_min_weight

                eb = network.edge_betweenness()
                ebs_min_ovh = np.array([eb[e.index] for e in es_min_ovh])
                min_eb = np.min(ebs_min_ovh)
                enames_curr_step = [es_min_ovh[i]['name'] for i in np.nonzero(ebs_min_ovh == min_eb)[0]]
                is_parellel = len(enames_curr_step) == 1

            elif 'matching' in strategy:
                if strategy == 'weight_and_matching':
                    min_weight = min(network.es['weight'])
                    es_min_weight = network.es.select(weight=min_weight)
                    subnetwork = network.subgraph_edges(es_min_weight)
                else:
                    subnetwork = network

                # for eid, loop in enumerate(subnetwork.is_loop()):
                #     if loop:
                #         if subnetwork is network:
                #             subnetwork = network.copy()
                #
                #         edge = subnetwork.es[eid]
                #         v_connected = edge.source_vertex
                #         v_temp = subnetwork.add_vertex()
                #         subnetwork.add_edge(v_connected, v_temp, name=edge['name'])

                # if use_gt:
                #     # t0 = time.time()
                #     subnetwork_gt = subnetwork.to_graph_tool(edge_attributes={'name': 'string'})
                #     # t_converting += time.time() - t0
                #     # t0 = time.time()
                #     matching = gt.max_cardinality_matching(subnetwork_gt, edges=True)
                #     # t_matching += time.time() - t0
                #     # t0 = time.time()
                #     es_curr_step = gt.find_edge(subnetwork_gt, matching, 1)
                #     enames = subnetwork_gt.ep['name']
                #     enames_curr_step = [enames[e] for e in es_curr_step]
                #     # t_postprocess += time.time() - t0

                # else:
                # t0 = time.time()
                subnetwork_nx = subnetwork.to_networkx()
                # t_converting += time.time() - t0
                # t0 = time.time()
                matching = nx.max_weight_matching(subnetwork_nx, weight=None)
                # t_matching += time.time() - t0
                # t0 = time.time()
                enames_curr_step = subnetwork.es[subnetwork.get_eids(matching)]['name']
                # t_postprocess += time.time() - t0

                is_parellel = True

            elif strategy == 'random':
                enames_curr_step = [np.random.choice(network.es['name'])]
                is_parellel = True

            else:
                raise ValueError

            recalculated_enames = []

            if get_fusion_order and not is_fusion_order_given:
                fusion_order_curr_step = set()
                fusion_order.append(fusion_order_curr_step)

            while True:
                if not is_parellel:
                    for rec_ename in recalculated_enames:
                        try:
                            enames_curr_step.remove(rec_ename)
                        except ValueError:
                            pass

                if not enames_curr_step:
                    break

                recalculated_enames.clear()

                if is_parellel:
                    ename_to_merge = enames_curr_step.pop()
                else:
                    ename_to_merge = np.random.choice(enames_curr_step)
                    enames_curr_step.remove(ename_to_merge)

                if get_fusion_order:
                    e_to_merge = self.fusion_network.es.find(name=ename_to_merge)
                    v1, v2 = e_to_merge.source_vertex, e_to_merge.target_vertex
                    fusion_order_curr_step.add((v1['name'], v2['name'], ename_to_merge))

                enames_updated_weight = self.contract_edge(network, ename_to_merge)

                recalculated_enames.extend(enames_updated_weight)

        v_final = network.vs.select(on=True)
        overhead = sum(v_final['weight'])
        step = max(v_final['step'])

        results = {
            'overhead': overhead,
            'step': step,
            'fusion_order_strategy': strategy,
            'p_succ': p_succ
        }

        if get_fusion_order:
            results['fusion_order'] = fusion_order

        self.data.update(results)

        # print(t_converting, '\n', t_matching, '\n', t_postprocess)

        return self.data

    # def _calculate_overhead_gt(self,
    #                            p_succ=0.5,
    #                            strategy='overhead_and_matching',
    #                            get_fusion_order=False,
    #                            ):
    #     network = self.fusion_network.to_graph_tool(vertex_attributes={'name': 'string'},
    #                                                 edge_attributes={'name': 'string'})
    #     network.vp['overhead'] = network.new_vp('double', val=1)
    #     network.vp['step'] = network.new_vp('int', val=0)
    #     network.ep['overhead'] = network.new_ep('double', val=2 / p_succ)
    #
    #     fusion_order = []
    #
    #     vids = network.vertex_index
    #     vnames = network.vp['name']
    #     enames = network.ep['name']
    #     v_overheads = network.vp['overhead']
    #     e_overheads = network.ep['overhead']
    #     steps = network.vp['step']
    #
    #     network.set_fast_edge_removal()
    #
    #     # Iterate until no edges remain in the fusion network
    #     while True:
    #         if not network.num_edges():
    #             break
    #
    #         # print("All edges:")
    #         # for e in network.edges():
    #         #     print(e, network.edge_index[e], e_overheads[e])
    #
    #         # Determine edges to merge parallelly in the current step
    #         if strategy == 'overhead_and_matching':
    #             overheads = e_overheads.fa
    #             # print(f"overheads = {overheads}")
    #             # print('num_edges =', network.num_edges())
    #             # for e in network.edges():
    #             #     print(e, network.edge_index[e], e_overheads[e])
    #             # print(e_overheads.a)
    #             # print(e_overheads.fa)
    #             # print()
    #             subnetwork = gt.GraphView(network, efilt=(overheads == overheads.min()))
    #             subnetwork = subnetwork.copy()
    #             subnetwork.purge_edges()
    #         else:
    #             subnetwork = network
    #
    #         # if network.num_edges() == 3:
    #         #     return subnetwork
    #
    #         eprop_matching = gt.max_cardinality_matching(subnetwork, edges=True)
    #         es_curr_step_subnetwork = gt.find_edge(subnetwork, eprop_matching, 1)
    #         # es_curr_step = [gt.find_edge(network, network.edge_index, subnetwork.edge_index[e_sn])
    #         #                 for e_sn in es_curr_step_subnetwork]
    #         es_curr_step = []
    #         for e_sn in es_curr_step_subnetwork:
    #             eid = subnetwork.edge_index[e_sn]
    #             e = gt.find_edge(network, network.edge_index, eid)[0]
    #             es_curr_step.append(e)
    #         #     vid1 = subnetwork.vertex_index[e_sn.source()]
    #         #     vid2 = subnetwork.vertex_index[e_sn.target()]
    #         #     es = network.edge(vid1, vid2, all_edges=True)
    #         #     if len(es) == 1:
    #         #         e = es[0]
    #         #     else:
    #         #         ename = subnetwork.ep['name'][e_sn]
    #         #         for e_ in es:
    #         #             if enames[e_] == ename:
    #         #                 e = e_
    #         #                 break
    #         #     es_curr_step.append(e)
    #
    #         # print("Subnetwork edges:")
    #         # for e in subnetwork.edges():
    #         #     print(e, subnetwork.edge_index[e], e_overheads[e])
    #         #
    #         # print(f"eprop = {eprop_matching.a}")
    #         #
    #         # print("Selected edges:")
    #         # for e in es_curr_step:
    #         #     print(e, network.edge_index[e])
    #         # print()
    #         #
    #         # time.sleep(1)
    #
    #         if get_fusion_order:
    #             fusion_order_curr_step = set()
    #             fusion_order.append(fusion_order_curr_step)
    #
    #         vs_to_remove = []
    #         for e_to_merge in es_curr_step:
    #             ename_to_merge = enames[e_to_merge]
    #
    #             v_merged, v_removed = e_to_merge.source(), e_to_merge.target()
    #             assert v_merged != v_removed
    #             deg_merged, deg_removed = network.get_out_degrees([v_merged, v_removed])
    #             if deg_merged < deg_removed:
    #                 v_merged, v_removed = v_removed, v_merged
    #             vid_merged = vids[v_merged]
    #             vid_removed = vids[v_removed]
    #             vname_merged = vnames[v_merged]
    #             vname_removed = vnames[v_removed]
    #
    #             vs_to_remove.append(v_removed)
    #
    #             if get_fusion_order:
    #                 fusion_order_curr_step.add((vname_merged, vname_removed, ename_to_merge))
    #
    #             # Contract the edge
    #             v_overheads[v_merged] = e_overheads[e_to_merge]
    #             updated_step = max(steps[v_merged], steps[v_removed]) + 1
    #             steps[v_merged] = updated_step
    #             self.fusion_network.es.find(name=ename_to_merge)['step'] = updated_step
    #
    #             for e_to_remove in v_removed.out_edges():
    #                 if e_to_remove != e_to_merge:
    #                     v1, v2 = e_to_remove.source(), e_to_remove.target()
    #                     vids_new_edge = [vid_merged if vid_ == vid_removed else vid_
    #                                      for vid_ in [vids[v1], vids[v2]]]
    #
    #                     if vids_new_edge[0] == vids_new_edge[1]:
    #                         v_vrt = network.add_vertex()
    #                         new_edge = network.add_edge(v_merged, v_vrt)
    #                         v_overheads[v_vrt] = 0
    #                     else:
    #                         new_edge = network.add_edge(min(vids_new_edge), max(vids_new_edge))
    #
    #                     enames[new_edge] = enames[e_to_remove]
    #
    #                 network.remove_edge(e_to_remove)
    #
    #             for e_to_update in v_merged.out_edges():
    #                 v1, v2 = e_to_update.source(), e_to_update.target()
    #                 e_overheads[e_to_update] = (v_overheads[v1] + v_overheads[v2]) / p_succ
    #
    #         network.remove_vertex(vs_to_remove)
    #
    #     overhead = float(v_overheads.a.sum())
    #     step = int(steps.a.max())
    #
    #     results = {
    #         'overhead': overhead,
    #         'step': step,
    #         'fusion_order_strategy': strategy,
    #         'p_succ': p_succ
    #     }
    #
    #     if get_fusion_order:
    #         results['fusion_order'] = fusion_order
    #
    #     self.data.update(results)
    #
    #     return self.data

    def simulate(self,
                 n_samples=1,
                 p_succ=0.5,
                 mp=False,
                 n_procs=None,
                 get_all_data=False,
                 get_all_graphs=False,
                 get_all_fusion_networks=False,
                 unravel=True,
                 unravel_bcs_first='random',
                 fusion_order_strategy='weight_and_matching',
                 seed='keep',
                 verbose=False,
                 pbar=False,
                 **kwargs
                 ):

        t0 = time.time()

        if seed != 'keep':
            np.random.seed(seed)

        if mp:
            if n_procs is None:
                n_procs = os.cpu_count()
            mp = mp and n_samples >= n_procs

        if not mp:
            if verbose:
                print("No multiprocessing.")
                print(f"Calculating for n_samples = {n_samples}... ", end='')

            if n_samples == 1 and seed is not None and seed != 'keep':
                seeds_samples = [seed]
            else:
                seeds_samples = np.random.randint(0, _max_seed(), size=n_samples)

            overheads = [] if get_all_data else None
            steps = [] if get_all_data else None
            seeds = [] if get_all_data else None
            unravalled_graphs = [] if get_all_graphs else None
            fusion_networks = [] if get_all_fusion_networks else None

            best_sample = None
            lowest_overhead = None
            for i_sample in range(n_samples):
                seed_sample = seeds_samples[i_sample]
                np.random.seed(seed_sample)

                if unravel:
                    self.unraveled_graph = None
                    try:
                        self.unravel_graph(unravel_bcs_first=unravel_bcs_first)
                    except:
                        print('Error occurs during unraveling')
                        print('seed =', seed_sample)
                        raise ValueError

                try:
                    self.build_fusion_network(use_unraveled_graph=unravel)
                except:
                    print('Error occurs during building fusion network')
                    print('seed =', seed_sample)
                    raise ValueError

                try:
                    data_now = self.calculate_overhead(p_succ=p_succ,
                                                       strategy=fusion_order_strategy,
                                                       **kwargs)
                except:
                    print('Error occurs during calculating overhead')
                    print('seed =', seed_sample)
                    raise ValueError

                overhead_now = data_now['overhead']
                step_now = data_now['step']
                self.data['seed'] = seed_sample

                if lowest_overhead is None or overhead_now < lowest_overhead:
                    best_sample = i_sample
                    lowest_overhead = overhead_now
                    best_ogs = self.copy()

                if get_all_data:
                    overheads.append(overhead_now)
                    steps.append(step_now)
                    seeds.append(seed_sample)

                if get_all_graphs:
                    unravalled_graphs.append(self.unraveled_graph)

                if get_all_fusion_networks:
                    fusion_networks.append(self.fusion_network)

            # res = best_ogs.data
            # res['n_samples'] = n_samples
            # res['unravel'] = unravel

            res = {
                'best_overhead': best_ogs.data['overhead'],
                'best_step': best_ogs.data['step'],
                'best_seed': best_ogs.data['seed'],
                # 'best_unraveled_graph': best_ogs.unraveled_graph,
                # 'best_fusion_network': best_ogs.fusion_network,
                'n_samples': n_samples,
                'p_succ': p_succ,
                'unravel': unravel,
                'fusion_order_strategy': fusion_order_strategy,
            }

            if unravel:
                res['unravel_bcs_first'] = best_ogs.data['unravel_bcs_first']

            if get_all_data or get_all_graphs or get_all_fusion_networks:
                res['best_sample'] = best_sample

                if get_all_data:
                    res['overheads'] = overheads
                    res['steps'] = steps
                    res['seeds'] = seeds

                if get_all_graphs:
                    res['unraveled_graphs'] = unravalled_graphs

                if get_all_fusion_networks:
                    res['fusion_networks'] = fusion_networks

        else:
            if 'parmap' not in sys.modules:
                raise ModuleNotFoundError("Package parmap is not installed.")

            if verbose:
                print(f"Use multiprocessing (n_procs = {n_procs})")
                print(f"Calculating for n_samples = {n_samples}, n_procs = {n_procs}... ", end='')

            additional_keys = []
            if get_all_data:
                additional_keys.extend(['overheads', 'steps', 'seeds'])
            if get_all_graphs:
                additional_keys.append('unraveled_graphs')
            if get_all_fusion_networks:
                additional_keys.append('fusion_networks')

            left = n_samples % n_procs
            ns_samples = [n_samples // n_procs] * n_procs
            for i in range(left):
                ns_samples[i] += 1

            seeds = np.random.randint(0, _max_seed(), size=n_procs)

            res_procs = parmap.starmap(_simulate_single,
                                       list(zip(ns_samples, seeds)),
                                       self.graph,
                                       p_succ=p_succ,
                                       get_all_data=get_all_data,
                                       get_all_graphs=get_all_graphs,
                                       get_all_fusion_networks=get_all_fusion_networks,
                                       unravel=unravel,
                                       unravel_bcs_first=unravel_bcs_first,
                                       fusion_order_strategy=fusion_order_strategy,
                                       pm_pbar=pbar,
                                       **kwargs)
            best_overheads = [res_each['best_overhead'] for res_each in res_procs]
            best_proc = np.argmin(best_overheads)
            res = res_procs[best_proc]
            res['n_samples'] = n_samples
            best_ogs = res['best_ogs']
            del res['best_ogs']

            if additional_keys:
                res['best_sample'] += sum(ns_samples[:best_proc])

            for key in additional_keys:
                vals = [res_each[key] for res_each in res_procs]
                res[key] = list(itertools.chain(*vals))

        if verbose:
            print(f"Done. Best: {res['best_overhead']:.2f} ({time.time() - t0:.2f} s)")

        self.unraveled_graph = best_ogs.unraveled_graph
        self.fusion_network = best_ogs.fusion_network
        self.unraveled_bcss = best_ogs.unraveled_bcss
        self.unraveled_cliques = best_ogs.unraveled_cliques
        self.data = best_ogs.data

        return res

    def simulate_adaptive(self,
                          init_n_samples,
                          p_succ=0.5,
                          mul=2,
                          mp=False,
                          n_procs=None,
                          get_all_data=False,
                          get_all_graphs=False,
                          get_all_fusion_networks=False,
                          unravel=True,
                          unravel_bcs_first='random',
                          fusion_order_strategy='weight_and_matching',
                          seed='keep',
                          verbose=True,
                          pbar=True,
                          **kwargs
                          ):
        if mp and n_procs is None:
            n_procs = os.cpu_count()

        if seed != 'keep':
            np.random.seed(seed)

        additional_keys = []
        if get_all_data:
            additional_keys.extend(['overheads', 'steps'])
        if get_all_graphs:
            additional_keys.append('unraveled_graphs')
        if get_all_fusion_networks:
            additional_keys.append('fusion_networks')

        if verbose:
            if mp:
                print(f"Multiprocessing (n_procs = {n_procs})")
            else:
                print("No multiprocessing")

        n_samples_history = []
        n_samples_now = init_n_samples
        res = None
        while True:
            if verbose:
                print(f"Calculating for n_samples = {n_samples_now}... ", end='')
            t0 = time.time()

            n_samples_history.append(n_samples_now)
            res_now = self.simulate(n_samples=n_samples_now,
                                    p_succ=p_succ,
                                    mp=mp,
                                    n_procs=n_procs,
                                    get_all_data=get_all_data,
                                    get_all_graphs=get_all_graphs,
                                    get_all_fusion_networks=get_all_fusion_networks,
                                    unravel=unravel,
                                    unravel_bcs_first=unravel_bcs_first,
                                    fusion_order_strategy=fusion_order_strategy,
                                    verbose=False,
                                    pbar=pbar,
                                    **kwargs)

            if res is None:
                res = res_now
                best_ogs = self.copy()
                n_samples_now *= mul

            else:
                for key in additional_keys:
                    res[key].extend(res_now[key])

                if res_now['best_overhead'] < res['best_overhead']:
                    for key in additional_keys:
                        res_now[key] = res[key]
                    res = res_now
                    best_ogs = self.copy()

                    n_samples_now *= mul

                else:
                    if verbose:
                        print(f"Done. Best: {res['best_overhead']:.2f} ({time.time() - t0:.2f} s)")
                    break

            if verbose:
                print(f"Done. Best: {res['best_overhead']:.2f} ({time.time() - t0:.2f} s)")

        res['n_samples'] = sum(n_samples_history)

        if additional_keys:
            res['best_sample'] += res['best_sample']

        self.unraveled_graph = best_ogs.unraveled_graph
        self.fusion_network = best_ogs.fusion_network
        self.unraveled_bcss = best_ogs.unraveled_bcss
        self.unraveled_cliques = best_ogs.unraveled_cliques
        self.data = best_ogs.data

        return res

    def plot_graph(self,
                   unraveled=False,
                   **kwargs):

        graph = self.unraveled_graph if unraveled else self.graph
        if graph is None:
            raise ValueError("No unraveled graph created.")

        fig, ax = plot_graph(graph, **kwargs)

        return fig, ax

    def plot_fusion_network(self,
                            inter=False,
                            **kwargs):
        network = self.inter_fusion_network if inter else self.fusion_network
        if network is None:
            if inter:
                raise ValueError('No intermediate fusion network created.')
            else:
                raise ValueError('No fusion network created.')

        fig, ax = plot_fusion_network(network, **kwargs)

        return fig, ax

    def copy(self):
        ogs = OptGraphState(graph=self.graph,
                            unraveled_graph=self.unraveled_graph,
                            fusion_network=self.fusion_network)
        # ogs.inter_fusion_network = self.inter_fusion_network
        ogs.unraveled_bcss = self.unraveled_bcss
        ogs.unraveled_cliques = self.unraveled_cliques
        ogs.data = self.data.copy()
        if self.graph_info is not None:
            ogs.graph_info = self.graph_info.copy()

        return ogs


def _simulate_single(n_samples,
                     seed,
                     graph,
                     **kwargs):
    # To ensure that different processes have different random seeds.
    ogs = OptGraphState(graph=graph)
    res = ogs.simulate(n_samples=n_samples, seed=seed, **kwargs)
    res['best_ogs'] = ogs

    return res
