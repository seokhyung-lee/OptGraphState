import itertools
import random

import numpy as np
import igraph as ig


def get_graph_from_edges(edges):
    """
    Generate an [`igraph.Graph`](https://python.igraph.org/en/stable/api/igraph.Graph.html) object from the list of edges.

    Parameters
    ----------
    edges : list of 2-tuple of int
        List of edges that form the graph, where each integer indicates a
        vertex label.

    Returns
    -------
    graph : igraph.Graph
        Generated graph.

    """
    return ig.Graph(edges=edges)


def get_sample_graph(shape, *prms):
    """
    Generate a predefined graph with a given shape and parameters.

    See the description of `optgraphstate.GraphState.__init__()`.

    Returns
    -------
    graph : igraph.Graph
        Generated graph.
    """
    graph_info = {'shape': shape}

    if shape == 'random':
        assert len(prms) in [2, 3]
        if len(prms) == 3:
            random.seed(prms[2])
        g = ig.Graph.Erdos_Renyi(n=prms[0], m=prms[1])
        random.seed()

    elif shape == 'complete':
        assert len(prms) == 1
        g = ig.Graph.Full(n=prms[0])

    elif shape == 'star':
        assert len(prms) == 1
        g = ig.Graph.Star(prms[0])

    elif shape == 'linear':
        g = ig.Graph.Lattice(dim=prms, circular=False)

    elif shape == 'cycle':
        g = ig.Graph.Ring(prms[0])

    elif shape == 'lattice':
        g = ig.Graph.Lattice(dim=prms, circular=False)
        # graph_info['size'] = tuple(prms)

    elif shape == 'tree':
        # graph_info['num_children'] = tuple(prms)
        g: ig.Graph = ig.Graph.Star(prms[0] + 1)
        parents_start = 1
        n_parents = prms[0]
        for n_children in prms[1:]:
            for parent in range(parents_start, parents_start + n_parents):
                g.add_vertices(n_children)
                vcount = g.vcount()
                g.add_edges([(parent, child) for child in
                             range(vcount - n_children, vcount)])

            parents_start += n_parents
            n_parents *= n_children

    elif shape == 'rhg':
        assert len(prms) == 3
        g = _get_rhg_lattice(*prms)
        # g['size'] = tuple(prms)

    elif shape == 'repeater':
        assert len(prms) == 1
        g, _ = get_sample_graph('complete', 2 * prms[0])
        vcount = g.vcount()
        for v in range(vcount):
            new_v = g.add_vertex()
            g.add_edge(v, new_v)
        # g['m'] = prms[0]

    elif shape == 'parity_encoding':
        assert len(prms) == 3
        g = _get_parity_encoded_graph(*prms)
        # g['logical_graph'] = prms[0]
        # g['n'] = prms[1]
        # g['m'] = prms[2]

    elif shape == 'ptqc':
        assert len(prms) == 4
        g = _get_ptqc_graph(*prms)
        # g['n'] = prms[0]
        # g['m'] = prms[1]
        # g['HIC'] = prms[2]
        # g['central'] = prms[3]

    else:
        raise ValueError("Unsupported shape")

    return g


def find_nonoverlapping_bcss(g: ig.Graph, get_name=False):
    """
    Find a maximum set of bipartitely-complete subgraphs (BCSs) that do not
    share any vertices.

    A *maximum* set means that it cannot be enlarged by adding another BCS.
    The obtained set may vary each time the function is called since the
    iterations of vertices are randomized by `numpy.random`.

    Parameters
    ----------
    g : igraph.Graph
        Traget graph.

    get_name : bool (default: False)
        Whether to get the names or indices of vertices.

    Returns
    -------
    bcss : list of 2-tuple of list of {int or str}
        Each element corresponds to a BCS found and has the structure of
        ([v1, v2, ...], [u1, u2, ...]), where v1, v2, ... are the
        indices/names of the vertices in one part of the BCS and
        u1, u2, ... are the indices/names of the vertices in another part.
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
                if ngh_vid1 >= ngh_vid2 or ngh_vid2 in vids_in_bcs or (
                        vid, ngh_vid2) in edges_checked:
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
                    bcss.append((g.vs[part1]['name'], g.vs[part2]['name']))
                else:
                    bcss.append((part1, part2))

                vids_in_bcs.update(part1)
                vids_in_bcs.update(part2)

        for ngh_vid in ngh_vids:
            edges_checked.update({(vid, ngh_vid), (ngh_vid, vid)})

    return bcss


def find_nonoverlapping_cliques(g: ig.Graph, get_name=False):
    """
    Find a maximum set of cliques that do not share any vertices.

    A *maximum* set means that it cannot be enlarged by adding another clique.
    The obtained set may vary each time the function is called since the
    iterations of vertices are randomized by `numpy.random`.

    Parameters
    ----------
    g : igraph.Graph
        Target graph.

    get_name : bool (default: False)
        Whether to get the names or indices of vertices.

    Returns
    -------
    cliques : list of set of {int or str}
        Each element is the set of the indices/names of vertices in a clique.
    """
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


# def get_all_vertices(graph):
#     """
#     Get all the vertex names of a graph.
#
#     Parameters
#     ----------
#     graph : igraph.Graph
#         Target graph.
#
#     Returns
#     -------
#     vertices : list of str
#         List of the names of the vertices in `graph`.
#     """
#
#     return graph.vs['name']
#
#
# def get_adjacency(graph):
#     """
#     Get the adjacency matrix of a graph.
#
#     Parameters
#     ----------
#     graph : igraph.Graph
#         Target graph.
#
#     Returns
#     -------
#     adjacency : numpy.ndarray
#         Adjacency matrix of `graph`.
#     """
#     adj = graph.get_adjacency()
#     return np.array(list(adj))
#
#
# def get_vertex_attrs(graph, vertex):
#     """
#     Get the attributes of a vertex of a graph.
#
#     Parameters
#     ----------
#     graph : igraph.Graph
#         Graph that the vertex belongs to.
#
#     vertex : str or int
#         Name of the vertex.
#
#     Returns
#     -------
#     attributes : dict
#         Attributes of the vertex.
#     """
#     vertex = str(vertex)
#     attrs = graph.vs.find(name=vertex).attributes()
#
#     return attrs
#
#
# def get_neighbors(graph, vertex):
#     """
#     Get the neighbors of a vertex in a graph.
#
#     Parameters
#     ----------
#     graph: igraph.Graph
#         Graph that the vertex belongs to.
#
#     vertex: str or int
#         Name of the vertex.
#
#     Returns
#     -------
#     neighbors: list of str
#         List of the names of the neighbors of the vertex.
#     """
#
#     vertex = str(vertex)
#
#     return graph.vs[graph.neighbors(str(vertex))]['name']
#
#
# def get_all_edges(graph):
#     """
#     Get all the edges of a graph in terms of pairs of vertex names.
#
#     Parameters
#     ----------
#     graph : igraph.Graph.
#         Target graph.
#
#     Returns
#     -------
#     edges : list of 2-tuple of str
#         List of the connected pairs of vertex names.
#     """
#
#     edges = []
#     for e in graph.es:
#         edges.append((e.source_vertex['name'], e.target_vertex['name']))
#     return edges
#
#
# def get_edge_attrs(graph, v1, v2):
#     """
#     Get the attributes of an edge in a graph.
#
#     Parameters
#     ----------
#     graph : igraph.Graph.
#         Graph that the edge belongs to.
#     v1, v2 : str or int
#         Names of the vertices connected by the edge.
#
#     Returns
#     -------
#     attributes : dict
#         Dictionary that contains the attributes of the edge.
#     """
#     n1 = str(v1)
#     n2 = str(v2)
#
#     link_attrs = graph.es[graph.get_eid(v1, v2)].attributes()
#     return link_attrs


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
                vids_not1_1 = range(vid_start + m,
                                    vid_start + (n - 1) * m + 1,
                                    m)
                _connect_vertex_sets(graph, vids_1_all, vids_not1_1)

                if m > 1:
                    for vid_x_1 in vids_not1_1:
                        vids_x_not1 = range(vid_x_1 + 1, vid_x_1 + m)
                        _connect_vertex_sets(graph, [vid_x_1], vids_x_not1)

                graph.vs[vids_not1_1]['clifford'] = 'H'

            return vids_1_all

    if hic:
        inbetween_edge_each_block = [True, True] if center else [True, False,
                                                                 False]

    else:
        inbetween_edge_each_block = [False, False] if center else [True, True,
                                                                   False]

    vids_lqs = [build_logical_qubit(lq, inbetween_edge_each_block[lq]) for lq
                in range(2 if center else 3)]

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

        # For each block besides the first one, the first vertex is
        # connected with the other vertices.
        for first_vid_each_block in range(first_vid + m, first_vid + n * m, m):
            other_vids = range(first_vid_each_block + 1,
                               first_vid_each_block + m)
            _connect_vertex_sets(g, [first_vid_each_block], other_vids)

        # The first vertex of each block has the Hadamard gate.
        g.vs[first_vids_each_block]['clifford'] = 'H'

    # Connection between logical qubits
    for logical_edge in logical_graph.es:
        logical_vids = logical_edge.source, logical_edge.target
        vids_first_blocks = [
            range(logical_vid * n * m, logical_vid * n * m + m) for logical_vid
            in logical_vids]
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
            return not (x + y + z) % 2 and (
                    x % 2 or y % 2 or z % 2)  # one even, two odds
        else:
            return (x + y + z) % 2 and not (
                    x % 2 and y % 2 and z % 2)  # two evens, one odd

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
        edges = [(vertex, get_vertex_by_coords(*adj_coords)) for adj_coords in
                 _adjacent_coords(vertex["x"], vertex["y"], vertex["z"]) if
                 duality_condition(*adj_coords, False)]
        g.add_edges(edges)

    return g
