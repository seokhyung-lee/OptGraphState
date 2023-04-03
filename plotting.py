import igraph as ig
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch


def _plot_base(g: ig.Graph, ax=None, layout='auto', figsize=(10, 10), **visual_style):
    if not isinstance(g, ig.Graph):
        raise TypeError("Parameter 'g' should be an igraph.Graph instance.")

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    ig.plot(g, target=ax, layout=layout, **visual_style)

    return ax


def plot_graph(graph,
               unravel=False,
               ax=None,
               layout='auto',
               figsize=(7, 7),
               vertex_color_normal='white',
               vertex_color_clifford='orange',
               vertices_to_highlight=None,
               vertex_color_highlight='purple',
               edge_color_normal='black',
               edge_color_fusion='red',
               **kwargs):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if isinstance(graph, dict):
        key = 'unraveled_graph' if unravel else 'graph'
        graph = graph[key]

    if not isinstance(graph, ig.Graph):
        raise TypeError("Parameter 'graph' is not ig.Graph.")

    unraveled = 'ext_fusion' in graph.vs.attributes()

    if vertices_to_highlight is None:
        vertices_to_highlight = []

    if unraveled:
        graph = graph.copy()
        org_ecount = graph.ecount()

        done = set()
        vs_with_ext_fusion = graph.vs.select(ext_fusion_ne=None)
        for v in vs_with_ext_fusion:
            if v['name'] not in done:
                vname_fusion = v['ext_fusion']
                graph.add_edge(v, vname_fusion)
                done.add(vname_fusion)
        # unlabelled_qubits = set()
        # ext_fusions = self._ext_fusions_from_cliques.copy()
        # ext_fusions.update(self._ext_fusions_from_quads)
        # graph.add_edges(ext_fusions.items())
        # for vname1, vname2 in ext_fusions.items():
        #     graph.add_edge(vname1, vname2)
        # unlabelled_qubits.update([graph.vs.find(name=vname1).index, graph.vs.find(name=vname2).index])

        vertex_colors = []
        for v in graph.vs:
            if v['name'] in vertices_to_highlight or v.index in vertices_to_highlight:
                color = vertex_color_highlight
            elif v['clifford'] is not None:
                color = vertex_color_clifford
            else:
                color = vertex_color_normal
            vertex_colors.append(color)

        visual_style = {
            # 'vertex_label': ['' if v.index in unlabelled_qubits else v['name'] for v in graph.vs],
            'vertex_label': graph.vs['name'],
            'vertex_color': vertex_colors,
            'edge_color': [edge_color_normal] * org_ecount + [
                edge_color_fusion] * round(len(vs_with_ext_fusion) / 2),
            # 'edge_width': edge_width,
        }

        visual_style.update(kwargs)
        _plot_base(graph, ax=ax, layout=layout, **visual_style)

        children = ax.get_children()
        vcount = graph.vcount()
        ecount = graph.ecount()
        lines = children[2 * vcount:2 * vcount + ecount]
        for eid in range(org_ecount, ecount):
            lines[eid].set(linestyle='--')

    else:
        graph = graph
        visual_style = {
            'vertex_label': list(range(graph.vcount())),
            'vertex_color': vertex_color_normal,
            'edge_color': edge_color_normal,
            # 'edge_width': edge_width,
        }

        visual_style.update(kwargs)
        _plot_base(graph, ax=ax, layout=layout, **visual_style)

    return fig, ax


def plot_fusion_network(network,
                        ax=None,
                        layout='auto',
                        figsize=(7, 7),
                        show_vertex_overhead=False,
                        show_edge_overhead=False,
                        show_fusion_order=True,
                        show_edge_labels=False,
                        vertex_color_normal='white',
                        vertex_color_clifford='orange',
                        uniform_edge_style=False,
                        edge_color_LL='black',
                        edge_color_RL='blue',
                        edge_color_RR='red',
                        edge_style_LL='-',
                        edge_style_RL='-',
                        edge_style_RR='--',
                        **kwargs):
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if isinstance(network, dict):
        network = network['fusion_network']

    if not isinstance(network, ig.Graph):
        raise TypeError("Parameter 'network' is not ig.Graph.")

    show_vertex_overhead = show_vertex_overhead and 'overhead' in network.vs.attributes()
    show_edge_overhead = show_edge_overhead and 'overhead' in network.es.attributes()
    show_fusion_order = show_fusion_order and 'step' in network.es.attributes()

    # if ignore_isolated_vertices:
    #     isolated_vs = network.vs.select(_degree=0)
    #     if isolated_vs:
    #         network = network.copy()
    #         network.delete_vertices(isolated_vs)

    # vcount = network.vcount()
    # ecount = network.ecount()

    if not uniform_edge_style and network.ecount() and 'RL' in network.es['kind']:
        network_directed = network.copy()
        network_directed.to_directed()

        es_to_delete = []
        for edge in network.es:
            v1, v2 = edge.source_vertex, edge.target_vertex
            vname1, vname2 = v1['name'], v2['name']
            if edge['kind'] == 'RL':
                root_node = edge['root_node']
                e_to_delete = (vname1, vname2) if root_node == vname2 else (vname2, vname1)
            else:
                e_to_delete = (vname2, vname1)
            es_to_delete.append(e_to_delete)

        network_directed.delete_edges(es_to_delete)
        network = network_directed

    edge_color_dict = {'RR': edge_color_RR, 'RL': edge_color_RL, 'LL': edge_color_LL}
    cliffords = [v['clifford_root'] is not None or v['clifford_leaves'] is not None
                 for v in network.vs]
    if show_vertex_overhead:
        vertex_label = []
        for v in network.vs:
            name = v['name']
            overhead = v['overhead']
            if abs(overhead % 1) < 1e-6:
                overhead = round(overhead)
            vertex_label.append(f'{name}:{overhead}')
    else:
        vertex_label = network.vs['name'] if network.vcount() else []

    visual_style = {
        'vertex_color': [vertex_color_clifford if clifford else vertex_color_normal for clifford in cliffords],
        'vertex_shape': ['square' if clifford else 'circle' for clifford in cliffords],
        'vertex_label': vertex_label,
        'edge_align_label': True,
    }

    if network.ecount() and not uniform_edge_style:
        visual_style['vertex_size']: 0.5
        visual_style['edge_arrow_size'] = [0.02 if e['kind'] == 'RL' else 0 for e in network.es]
        visual_style['edge_color'] = [edge_color_dict[e['kind']] for e in network.es]

    if network.ecount() and (show_fusion_order or show_edge_overhead or show_edge_labels):
        visual_style['edge_background'] = 'white'

        if sum([show_fusion_order, show_edge_overhead, show_edge_labels]) == 1:
            if show_fusion_order:
                key = 'step'
            elif show_edge_overhead:
                key = 'overhead'
            else:
                key = 'name'
            edge_label = network.es[key]

        else:
            edge_label = []
            for e in network.es:
                label = []
                if show_edge_labels:
                    label.append(f"L{e['name']}")
                if show_fusion_order:
                    label.append(f"O{e['step']}")
                if show_edge_overhead:
                    overhead = e['overhead']
                    if abs(overhead % 1) < 1e-6:
                        overhead = round(overhead)
                    label.append(f"W{overhead}")
                label = '_'.join(label)
                edge_label.append(label)

        visual_style['edge_label'] = edge_label

    visual_style.update(kwargs)

    _plot_base(network, ax=ax, layout=layout, **visual_style)

    if not uniform_edge_style and network.ecount():
        children = ax.get_children()
        lines = [child for child in children if isinstance(child, PathPatch)]

        # if network.is_directed():
        #     lines = children[vcount:vcount + 2 * ecount:2]
        # else:
        #     lines = children[2 * vcount + ecount:2 * vcount + 2 * ecount]

        edge_style = {'RR': edge_style_RR, 'RL': edge_style_RL, 'LL': edge_style_LL}
        for eid, line in enumerate(lines):
            line.set(linestyle=edge_style[network.es[eid]['kind']])

    return fig, ax
