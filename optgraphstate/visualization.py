import igraph as ig
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch


def _plot_base(g: ig.Graph,
               ax=None,
               layout='auto',
               figsize=(10, 10),
               **visual_style):
    if not isinstance(g, ig.Graph):
        raise TypeError("Parameter 'g' should be an igraph.Graph instance.")

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    ig.plot(g, target=ax, layout=layout, **visual_style)

    return ax


def plot_graph(graph,
               ax=None,
               layout='auto',
               figsize=(5, 5),
               save=None,
               show_vertex_names=True,
               vertex_color_normal='white',
               vertex_color_clifford='orange',
               vertices_to_highlight=None,
               vertex_color_highlight='purple',
               edge_color_normal='black',
               edge_color_fusion='red',
               edge_style_fusion='--',
               **kwargs):
    """
    Plot a graph or an unraveled graph.

    See the description of `optgraphstate.GraphState.plot_graph()` for details.

    Parameters
    ----------
    graph : igraph.Graph
        Graph or unraveled graph to plot.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes object.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if not isinstance(graph, ig.Graph):
        raise TypeError("Parameter 'graph' is not ig.Graph.")

    unraveled = 'ext_fusion' in graph.vs.attributes()
    clifford_exists = 'clifford' in graph.vs.attributes()

    if vertices_to_highlight is None:
        vertices_to_highlight = []

    vertices_to_highlight = [str(vname) for vname in vertices_to_highlight]
    vertex_colors = []
    for v in graph.vs:
        if v['name'] in vertices_to_highlight:
            color = vertex_color_highlight
        elif clifford_exists and v['clifford'] is not None:
            color = vertex_color_clifford
        else:
            color = vertex_color_normal
        vertex_colors.append(color)

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

        visual_style = {
            'vertex_color': vertex_colors,
            'edge_color': [edge_color_normal] * org_ecount + [
                edge_color_fusion] * round(len(vs_with_ext_fusion) / 2), }

        if show_vertex_names:
            visual_style['vertex_label'] = graph.vs['name']

        visual_style.update(kwargs)
        _plot_base(graph, ax=ax, layout=layout, **visual_style)

        children = ax.get_children()
        ecount = graph.ecount()
        lines = [child for child in children if isinstance(child, PathPatch)]
        assert len(lines) == ecount
        for eid in range(org_ecount, ecount):
            lines[eid].set(linestyle=edge_style_fusion)

    else:
        graph = graph
        visual_style = {
            'vertex_color': vertex_colors,
            'edge_color': edge_color_normal
        }

        if show_vertex_names:
            visual_style['vertex_label'] = graph.vs['name']

        visual_style.update(kwargs)
        _plot_base(graph, ax=ax, layout=layout, **visual_style)

    plt.tight_layout()

    if save is not None:
        plt.savefig(save, transparent=True)

    return fig, ax


def plot_fusion_network(network,
                        ax=None,
                        layout='auto',
                        figsize=(5, 5),
                        save=None,
                        show_only_structure=False,
                        show_node_names=True,
                        node_color='white',
                        show_link_names=False,
                        show_fusion_orders=True,
                        show_link_cliffords=True,
                        link_color_ll='black',
                        link_color_rl='blue',
                        link_color_rr='red',
                        arrow_size=0.02,
                        **kwargs):
    """
    Plot a fusion network.

    See the description of `optgraphstate.GraphState.plot_fusion_network()`
    for details.

    Parameters
    ----------
    network : igraph.Graph
        Fusion network to plot.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes object.
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    if isinstance(network, dict):
        network = network['fusion_network']

    if not isinstance(network, ig.Graph):
        raise TypeError("Parameter 'network' is not ig.Graph.")

    if show_only_structure:
        show_node_names = False
        show_link_names = False
        show_fusion_orders = False
        show_link_cliffords = False

    # show_vertex_overhead = show_vertex_overhead and 'overhead' in
    # network.vs.attributes()
    # show_edge_overhead = show_edge_overhead and 'overhead' in
    # network.es.attributes()
    show_fusion_orders = show_fusion_orders and 'order' in network.es.attributes()
    rl_exists = 'RL' in network.es['kind']

    if rl_exists:
        network_directed = network.copy()
        network_directed.to_directed()

        es_to_delete = []
        for edge in network.es:
            v1, v2 = edge.source_vertex, edge.target_vertex
            vname1, vname2 = v1['name'], v2['name']
            if edge['kind'] == 'RL':
                root_node = edge['root_node']
                e_to_delete = (vname2, vname1) if root_node == vname2 \
                    else (vname1, vname2)
            else:
                e_to_delete = (vname2, vname1)
            es_to_delete.append(e_to_delete)

        network_directed.delete_edges(es_to_delete)
        network = network_directed

    link_color_dict = {
        'RR': link_color_rr,
        'RL': link_color_rl,
        'LL': link_color_ll
    }

    # cliffords = [
    #     v['clifford_root'] is not None or v['clifford_leaves'] is not None for
    #     v in network.vs]
    # if show_vertex_overhead:
    #     vertex_label = []
    #     for v in network.vs:
    #         name = v['name']
    #         overhead = v['overhead']
    #         if abs(overhead % 1) < 1e-6:
    #             overhead = round(overhead)
    #         vertex_label.append(f'{name}:{overhead}')
    # else:

    visual_style = {
        # 'vertex_color': [node_color_clifford if clifford else node_color
        #                  for clifford in cliffords],
        'vertex_color': node_color,
        # 'vertex_shape': ['square' if clifford else 'circle' for clifford in
        # cliffords],
        'edge_align_label': False,
        'edge_color': [link_color_dict[kind] for kind in network.es['kind']],
    }

    if show_node_names:
        visual_style['vertex_label'] \
            = network.vs['name'] if network.vcount() else []

    if rl_exists:
        visual_style['vertex_size'] = 0.4
        visual_style['edge_arrow_size'] = [
            arrow_size if e['kind'] == 'RL' else 0 for
            e in network.es]
        # visual_style['edge_color'] = [edge_color_dict[e['kind']] for e in
        #                               network.es]

    if network.ecount() and (show_fusion_orders or show_link_names):
        visual_style['edge_background'] = 'white'

        if sum([show_fusion_orders, show_link_names]) == 1:
            if show_fusion_orders:
                key = 'order'
            # elif show_edge_overhead:
            #     key = 'overhead'
            else:
                key = 'name'
            edge_label = network.es[key]

        else:
            edge_label = []
            for e in network.es:
                # label = []
                # if show_edge_labels:
                #     label.append(f"{e['name']}")
                # if show_fusion_order:
                #     label.append(f"{e['order']}")
                # if show_edge_overhead:
                #     overhead = e['overhead']
                #     if abs(overhead % 1) < 1e-6:
                #         overhead = round(overhead)
                #     label.append(f"W{overhead}")
                # label = '_'.join(label)
                label = f"{e['name']}-{e['order']}"
                edge_label.append(label)

        if show_link_cliffords:
            for link in network.es.select(cliffords_ne={}):
                eid = link.index
                edge_label[eid] = str(edge_label[eid]) + 'C'

        visual_style['edge_label'] = edge_label

    visual_style.update(kwargs)

    _plot_base(network, ax=ax, layout=layout, **visual_style)

    if network.ecount():
        children = ax.get_children()
        lines = [child for child in children if isinstance(child, PathPatch)]

        edge_style = {
            'RR': '--',
            'RL': '-',
            'LL': '-'
        }
        for eid, line in enumerate(lines):
            link = network.es[eid]
            line.set(linestyle=edge_style[link['kind']])

    plt.tight_layout()

    if save is not None:
        plt.savefig(save, transparent=True)

    return fig, ax
