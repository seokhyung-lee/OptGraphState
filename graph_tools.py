import itertools

import numpy as np
import igraph as ig


def get_graph(edges):
    return ig.Graph(edges=edges)


def get_sample_graph(shape, *prms):
    """
    Generate a predefined graph with a given shape and parameters.
    > 'complete': Full graph where prms is the number of vertices.
    > 'star': Star graph where prms is the number of leaf vertices (total prms + 1 vertices).
    > 'tree': Tree graph where all vertices at distance d from the root have prms[d] children.
    > 'rhg': Raussendorf-Harrington-Goyal lattice where
        - prms[0]: (int) size of the lattice along the x-axis in the unit of a cell,
        - prms[1]: (bool) whether the boundaries perpendicular to the x-axis is primal or not,
        - prms[2] (int) and prms[3] (bool) determine the size and boundaries for the y-axis,
        - prms[4] (int) and prms[5] (bool) determine the size and boundaries for the z-axis.

    Parameters
    ----------
    shape : string
        Determine the shape of the graph. One of {'complete', 'star', 'tree', 'tree_regular', 'rhg'}.
    prms : number or tuple of numbers
        Parameters for generating the graph.

    Returns
    -------
    graph : igraph.Graph.
        Generated graph.
    """
    graph_info = {'shape': shape}

    if shape == 'random':
        g = ig.Graph.Erdos_Renyi(n=prms[0], m=prms[1])

    elif shape == 'complete':
        g = ig.Graph.Full(n=prms[0])

    elif shape == 'star':
        g = ig.Graph.Star(prms[0])

    elif shape == 'linear':
        g = ig.Graph.Lattice(dim=prms, circular=False)

    elif shape == 'cycle':
        g = ig.Graph.Ring(prms[0])

    elif shape == 'lattice':
        g = ig.Graph.Lattice(dim=prms, circular=False)
        graph_info['size'] = tuple(prms)

    elif shape == 'tree':
        graph_info['num_children'] = tuple(prms)
        g: ig.Graph = ig.Graph.Star(prms[0] + 1)
        parents_start = 1
        n_parents = prms[0]
        for n_children in prms[1:]:
            for parent in range(parents_start, parents_start + n_parents):
                g.add_vertices(n_children)
                vcount = g.vcount()
                g.add_edges([(parent, child) for child in range(vcount - n_children, vcount)])

            parents_start += n_parents
            n_parents *= n_children

    elif shape == 'rhg':
        graph_info['size'] = tuple(prms)
        g = _get_rhg_lattice(*prms)

    elif shape == 'repeater':
        graph_info['m'] = prms[0]
        g, _ = get_sample_graph('complete', 2 * prms[0])
        vcount = g.vcount()
        for v in range(vcount):
            new_v = g.add_vertex()
            g.add_edge(v, new_v)

    elif shape == 'parity_encoding':
        graph_info['logical_graph'] = prms[0]
        graph_info['n'] = prms[1]
        graph_info['m'] = prms[2]
        g = _get_parity_encoded_graph(*prms)

    elif shape == 'ptqc':
        graph_info['n'] = prms[0]
        graph_info['m'] = prms[1]
        graph_info['HIC'] = prms[2]
        graph_info['central'] = prms[3]
        g = _get_ptqc_graph(*prms)

    else:
        raise ValueError("Unsupported shape")

    graph_info['num_vertices'] = g.vcount()
    graph_info['num_edges'] = g.ecount()

    return g, graph_info


def greedy_coloring(g: ig.Graph, vertex_selection_method='largest_first'):
    vcount = g.vcount()
    if not vcount:
        return {}

    if vertex_selection_method == 'largest_first':
        vids_sorted = np.argsort(g.vs.degree())[::-1]
    elif vertex_selection_method == 'random':
        vids_sorted = np.arange(vcount)
        np.random.shuffle(vids_sorted)
    else:
        raise ValueError

    colors = np.full(vcount, np.nan, dtype='int32')
    for vid in vids_sorted:
        colors_ngh = colors[g.neighbors(vid)]
        for color in itertools.count():
            if color not in colors_ngh:
                break
        colors[vid] = color

    return colors


# Find all quadrangles (cycle with length 4) in the graph that do not share any edges.
def find_separated_quadrangles(g: ig.Graph, get_name=False):
    quads = []

    # Edges that are gauranteed to be either contained in a quadrangle or not contained in any quadrangles.
    edges_marked = set()

    vids = np.arange(g.vcount())
    np.random.shuffle(vids)

    for vid in vids:
        ngh_vids = g.neighbors(vid)
        np.random.shuffle(ngh_vids)
        for ngh_vid1, ngh_vid2 in itertools.combinations(ngh_vids, r=2):
            if (vid, ngh_vid1) in edges_marked or (vid, ngh_vid2) in edges_marked:
                continue

            common_ngh_ngh_vids = set(g.neighbors(ngh_vid1)) & set(g.neighbors(ngh_vid2)) - {vid}
            cond = lambda vid_: (ngh_vid1, vid_) not in edges_marked and (ngh_vid2, vid_) not in edges_marked
            common_ngh_ngh_vids = [vid_ for vid_ in common_ngh_ngh_vids if cond(vid_)]

            num_common_ngh_ngh_vs = len(common_ngh_ngh_vids)
            if num_common_ngh_ngh_vs == 1:
                ngh_ngh_vid = common_ngh_ngh_vids[0]
            elif num_common_ngh_ngh_vs > 1:
                ngh_ngh_vid = np.random.choice(common_ngh_ngh_vids)
            else:
                continue

            edges_quadrangle = [(vid, ngh_vid1), (ngh_vid1, ngh_ngh_vid), (ngh_ngh_vid, ngh_vid2),
                                (ngh_vid2, vid)]
            edges_marked.update(edges_quadrangle + [t[::-1] for t in edges_quadrangle])

            quad = (vid, ngh_vid1, ngh_ngh_vid, ngh_vid2)
            if get_name:
                quad = tuple([g.vs[vid]['name'] for vid in quad])
            quads.append(quad)

        for ngh_vid in ngh_vids:
            edges_marked.update([(vid, ngh_vid), (ngh_vid, vid)])

    return quads


def find_all_nonoverlapping_bcs(g: ig.Graph, get_name=False):
    """
    Find bipartitely-complete subgraphs that do not share any vertices.

    Parameters
    ----------
    g : ig.Graph
        Input graph
    get_name : bool
        Whether to get the names or indices of vertices

    Return
    ------
    bcss : List of list of lists of integers or strings.
    List of bipartitely-complete subgraphs expressed as
    [[[vertices in one part of the first bcs], [vertices in the other part of the first bcs]],
      [vertices in one part of the second bcs], [vertices in the other part of the second bcs]],
      ...]
    """
    bcss = []

    vids_in_bcs = set()  # Vertices that are already contained in a bcs
    edges_checked = set()  # Edges that are not contained in any bcs

    vids = np.arange(g.vcount())
    np.random.shuffle(vids)

    for vid in vids:
        if vid in vids_in_bcs:
            continue

        ngh_vids = g.neighbors(vid)
        np.random.shuffle(ngh_vids)

        for ngh_vid1 in ngh_vids:
            if ngh_vid1 in vids_in_bcs or (vid, ngh_vid1) in edges_checked:
                continue

            for ngh_vid2 in ngh_vids:
                if ngh_vid1 >= ngh_vid2 or ngh_vid2 in vids_in_bcs or (vid, ngh_vid2) in edges_checked:
                    continue

                part1 = set(g.neighbors(ngh_vid1)) & set(g.neighbors(ngh_vid2))

                if len(part1) == 1:
                    continue

                marked = False
                for part1_vid in part1:
                    if part1_vid in vids_in_bcs:
                        marked = True
                        break
                if marked:
                    continue

                part1_neighbors = [set(g.neighbors(vid)) for vid in part1]
                part2 = set.intersection(*part1_neighbors)

                marked = False
                for part2_vid in part2:
                    if part2_vid in vids_in_bcs:
                        marked = True
                        break
                if marked:
                    continue

                if get_name:
                    bcss.append([g.vs[part1]['name'], g.vs[part2]['name']])
                else:
                    bcss.append([part1, part2])

                vids_in_bcs.update(part1)
                vids_in_bcs.update(part2)

        for ngh_vid in ngh_vids:
            edges_checked.update({(vid, ngh_vid), (ngh_vid, vid)})

    return bcss


# def find_complete_bipartite_subgraph(g: ig.Graph, get_name=False):
#     """
#     Find one complete bipartite subgraph.
#
#     Parameters
#     ----------
#     g : ig.Graph
#         Input graph
#     get_name : bool
#         Whether to get the names or indices of vertices
#
#     Return
#     ------
#     part1, part2 : Lists of the names or indices of the vertices in the two parts respectively
#     If no complete bipartite subgraphs are found, (None, None) is returned.
#     """
#
#     # Edges that are gauranteed to be not contained in any bcs.
#     edges_marked = set()
#
#     vids = np.arange(g.vcount())
#     np.random.shuffle(vids)
#
#     for vid in vids:
#         ngh_vids = g.neighbors(vid)
#         np.random.shuffle(ngh_vids)
#
#         for ngh_vid1, ngh_vid2 in itertools.combinations(ngh_vids, r=2):
#             if (vid, ngh_vid1) in edges_marked or (vid, ngh_vid2) in edges_marked:
#                 continue
#
#             part1 = set(g.neighbors(ngh_vid1)) & set(g.neighbors(ngh_vid2)) - {vid}
#             cond = lambda part1_vid: (ngh_vid1, part1_vid) not in edges_marked \
#                                      and (ngh_vid2, part1_vid) not in edges_marked
#             part1 = [part1_vid for part1_vid in part1 if cond(part1_vid)]
#             part1.append(vid)
#
#             if len(part1) == 1:
#                 continue
#
#             part1_neighbors = [set(g.neighbors(vid)) for vid in part1]
#             part2 = set.intersection(*part1_neighbors) - {ngh_vid1, ngh_vid2}
#             cond = lambda part2_vid: all([(part1_vid, part2_vid) not in edges_marked for part1_vid in part1])
#             part2 = [part2_vid for part2_vid in part2 if cond(part2_vid)]
#             part2.extend([ngh_vid1, ngh_vid2])
#
#             if get_name:
#                 return g.vs[part1]['name'], g.vs[part2]['name']
#             else:
#                 return part1, part2
#
#         for ngh_vid in ngh_vids:
#             edges_marked.update({(vid, ngh_vid), (ngh_vid, vid)})
#
#     return None, None


# Find all cliques (fully-connected subgraphs) in the graph where any pair of them share at most one photon.
# If multiple cliques share two or more photons, select randomly among the cliques with the largest size.
def find_all_nonoverlapping_cliques(g: ig.Graph, get_name=False):
    cliques = g.maximal_cliques(min=3)
    cliques = [set(clique) for clique in cliques]
    # num_cliques = len(cliques)

    # Remove overlapping cliques
    np.random.shuffle(cliques)
    nonoverlapping_cliques = []
    all_vids = set()
    for clique in cliques:
        if not (clique & all_vids):
            nonoverlapping_cliques.append(clique)
            all_vids.update(clique)

    cliques = nonoverlapping_cliques
    if get_name:
        cliques = [{g.vs[vid]['name'] for vid in clique} for clique in cliques]

    return cliques


def _connect_vertex_sets(graph, inds1, inds2, **attrs):
    if inds1 and inds2:
        graph.add_edges(itertools.product(inds1, inds2), attributes=attrs)


def _get_ptqc_graph(n, m, hic, center):
    if center:
        graph = ig.Graph(2 * n * m + 1, vertex_attrs={"clifford": None})
    else:
        graph = ig.Graph(3 * n * m, vertex_attrs={"clifford": None})

    def build_logical_qubit(lq, inbetween_edge_each_block):
        if center:
            vid_start = 1 + lq * n * m
        else:
            vid_start = lq * n * m

        if inbetween_edge_each_block:
            vids_all_1 = range(vid_start, vid_start + (n - 1) * m + 1, m)
            if m > 1:
                for vid_x_1 in vids_all_1:
                    vids_x_not1 = range(vid_x_1 + 1, vid_x_1 + m)
                    _connect_vertex_sets(graph, [vid_x_1], vids_x_not1)

            graph.vs[vids_all_1]['clifford'] = 'H'

            return vids_all_1

        else:
            vids_1_all = range(vid_start, vid_start + m)

            if n > 1:
                vids_not1_1 = range(vid_start + m, vid_start + (n - 1) * m + 1, m)
                _connect_vertex_sets(graph, vids_1_all, vids_not1_1)

                if m > 1:
                    for vid_x_1 in vids_not1_1:
                        vids_x_not1 = range(vid_x_1 + 1, vid_x_1 + m)
                        _connect_vertex_sets(graph, [vid_x_1], vids_x_not1)

                graph.vs[vids_not1_1]['clifford'] = 'H'

            return vids_1_all

    if hic:
        inbetween_edge_each_block = [True, True] if center else [True, False, False]

    else:
        inbetween_edge_each_block = [False, False] if center else [True, True, False]

    vids_lqs = [build_logical_qubit(lq, inbetween_edge_each_block[lq]) for lq in range(2 if center else 3)]

    if center:
        _connect_vertex_sets(graph, [0], itertools.chain(*vids_lqs))

    else:
        _connect_vertex_sets(graph, vids_lqs[0], vids_lqs[1])
        _connect_vertex_sets(graph, vids_lqs[1], vids_lqs[2])

    return graph


def _get_parity_encoded_graph(logical_graph: ig.Graph, n, m):
    logical_vcount = logical_graph.vcount()
    vcount = logical_vcount * n * m
    g = ig.Graph(vcount, vertex_attrs={"clifford": None})

    # Internal structure of each logical qubit
    for first_vid in range(0, vcount, n * m):
        # The first block is connected with the first vertex of each block.
        vids_first_block = range(first_vid, first_vid + m)
        first_vids_each_block = range(first_vid + m, first_vid + n * m, m)
        _connect_vertex_sets(g, vids_first_block, first_vids_each_block)

        # For each block besides the first one, the first vertex is connected with the other vertices.
        for first_vid_each_block in range(first_vid + m, first_vid + n * m, m):
            other_vids = range(first_vid_each_block + 1, first_vid_each_block + m)
            _connect_vertex_sets(g, [first_vid_each_block], other_vids)

        # The first vertex of each block has the Hadamard gate.
        g.vs[first_vids_each_block]['clifford'] = 'H'

    # Connection between logical qubits
    for logical_edge in logical_graph.es:
        logical_vids = logical_edge.source, logical_edge.target
        vids_first_blocks = [range(logical_vid * n * m, logical_vid * n * m + m) for logical_vid in logical_vids]
        _connect_vertex_sets(g, *vids_first_blocks)

    return g


def _get_rhg_lattice(Lx, Ly, Lz):
    g = ig.Graph()
    size = (2 * Lx, 2 * Ly, 2 * Lz)

    def _adjacent_coords(*coords):
        for axis, coord in enumerate(coords):
            for diff in [1, -1]:
                adjacent_vertex = list(coords[:])
                adjacent_vertex[axis] += diff
                if 0 <= adjacent_vertex[axis] <= size[axis]:
                    yield tuple(adjacent_vertex)

    def duality_condition(x, y, z, primal):
        if primal:
            return not (x + y + z) % 2 and (x % 2 or y % 2 or z % 2)  # one even, two odds
        else:
            return (x + y + z) % 2 and not (x % 2 and y % 2 and z % 2)  # two evens, one odd

    def add_qubit(x, y, z, primal):
        vertex = g.add_vertex(x=x, y=y, z=z, primal=primal)
        return vertex

    # Add vertices
    primal_qubits = []
    for x in range(size[0] + 1):
        for y in range(size[1] + 1):
            for z in range(size[2] + 1):
                if duality_condition(x, y, z, True):
                    new_qubit = add_qubit(x, y, z, True)
                    primal_qubits.append(new_qubit)
                elif duality_condition(x, y, z, False):
                    add_qubit(x, y, z, False)

    # Add edges
    get_vertex_by_coords = lambda x, y, z: g.vs.find(x=x, y=y, z=z)
    for vertex in primal_qubits:
        edges = [(vertex, get_vertex_by_coords(*adj_coords))
                 for adj_coords in _adjacent_coords(vertex["x"], vertex["y"], vertex["z"])
                 if duality_condition(*adj_coords, False)]
        g.add_edges(edges)

    return g
