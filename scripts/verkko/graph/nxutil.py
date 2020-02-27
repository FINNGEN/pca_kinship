# Richard Darst, July 2012

"""Module to handle networkx graphs and their interfaces.

"""

import collections
import numpy
import random

import networkx

def from_string(s):
    """Create a simple graph from a string.

    This is designed to be a simple way to create small graphs for
    test scripts, not a serious way of making graphs.

    Examples:

    - `from_string('1-2 2-3 3-1')`

      This makes a simple clique graph with nodes 1,2,3

    - `from_string('1-2-3-1')`

      When there are multiple edges specified without a space, add all
      the edges independently, in this case '1-2 2-3 3-1' are added.

    These are the only operations defined for now.
    """
    edges = s.split()
    g = networkx.Graph()
    for e in edges:
        e = e.split('-')
        for i in range(len(e)-1):
            g.add_edge(eval(e[i]), eval(e[i+1]))
    return g

def edges_between(g, nbunch1, nbunch2):
    """List of edges between two node sets.

    For all information on this function, see documentation on
    ``n_edges_between``.

    Return value: list of tuples ``(node1, node2)``
        List of edges between nbunch1 and nbunch2.  For directed
        graphs, tuples are (from, to).
    """
    import verkko.misc.testutil as testutil ; testutil.warn_untested()
    # switch so nbunch1 is the smaller set.
    if (isinstance(g, networkx.Graph) and not isinstance(g, networkx.DiGraph)
            and len(nbunch1) > len(nbunch2) ):
        nbunch1, nbunch2 = nbunch2, nbunch1
    edges = [ ]
    for n1 in nbunch1:
        for neighbor in g.adj[n1]:  # iterate neighbors (or successors)
            if neighbor in nbunch2:
                edges.append((n1, neighbor))
    return edges
def n_edges_between(g, nbunch1, nbunch2):
    """Number of edges betwen nbunches.

    This returns the edges that begin in nbunch1 and end in nbunch2.
    This interpertation holds true for directed and undirected graphs.
    If nbunch1==nbunch2 and the graph is undirected, returns **double
    the number of edges**.  Any edges with both ends in the
    intersection of nbunch1 and nbunch2 will be included twice for
    directed graphs.

    Return value:
        Set of all edges (tuples ``(n1, n2)``) between two sets of
        nodes nbunch1 and nbunch2.


    Complexity:
        time: ``O(len(nbunch1))``  or  ``O( min(len(nbunch1), len(nbunch2)) )``
            For undirected graphs, nbunch1 and nbunch2 will be can be
            switched for efficiency purposes.

            nbunch2 (or nbunch1, if switched) MUST be sets or support
            ``O(1)`` ``__contains__`` in order to achieve this optimal
            efficiency, otherwise complexity is
            ``O(len(nbunch1)*len(nbunch2))``

        space: O(1)

    See also:
        n_edges_between for a version that returns only the
        counts between two nbunches.
    """
    import verkko.misc.testutil as testutil ; testutil.warn_untested()
    n = 0
    # switch so nbunch1 is the smaller set.
    if (isinstance(g, networkx.Graph) and not isinstance(g, networkx.DiGraph)
            and len(nbunch1) > len(nbunch2) ):
        nbunch1, nbunch2 = nbunch2, nbunch1
    for n1 in nbunch1:
        for neighbor in g.adj[n1]:  # iterate neighbors (or successors)
            if neighbor in nbunch2:
                n += 1
    return n



def graphcolor(g, colormap=None, distance=1):
    """Return a coloring of a graph.

    This function returns a graph coloring.  It is designed to be fast
    and hopefully suitable for visualization, not give the fewest
    number of colors.

    Arguments:

    g: undirected networkx.Graph
        graph to be colored
    colormap: str, matplotlib colormap name.
        If given, should be a string listing a colormap name.  Instead
        of returning node->int map, return node->color tuple map.
    distance: int, default 1
        Nodes have unique colors for the n-th nearest neighbors,
        instead of only nearest neighbors.  Sometimes when
        visualizing, you may want extra colors for clarity.

    Returns:

    dict node->int
        indexes for colors (or color tuples).

    Complexity:
        roughly time: O(N nodes * avg degree)
        space: O(N nodes + avg degree)
    """
    colors = set()
    coloring = { }
    # Go through and make a greedy coloring.  Coloring-satisfying
    # integer for each node.
    for n in g:
        # 1-st nearest neighbors
        neighbors = set(g.neighbors(n))
        # n-th nearest neighbors
        for _ in range(distance-1):
            neighbors = set(y for neigh in neighbors
                            for y in g.neighbors(neigh))
        #
        used_colors = set(coloring[neigh] for neigh in neighbors
                           if neigh in coloring)
        avail_colors = colors - used_colors
        if avail_colors:
            color = random.choice(list(avail_colors))
        else:
            color = len(colors)
            colors.add(color)
        coloring[n] = color
    # This step goes through and re-randomizes choices, given the
    # minimal set of colors we found already.
    for n in g:
        # 1-st nearest neighbors
        neighbors = set(g.neighbors(n))
        # n-th nearest neighbors
        for _ in range(distance-1):
            neighbors = set(y for neigh in neighbors
                            for y in g.neighbors(neigh))

        used_colors = set(coloring[neigh] for neigh in neighbors)
        avail_colors = colors - used_colors
        avail_colors.add(coloring[n])
        color = random.choice(list(avail_colors))
        coloring[n] = color

    # If wanted, do a matplotlib colormap
    if colormap:
        import matplotlib.cm as cm
        import matplotlib.colors as mcolors

        colormap = cm.get_cmap(colormap)
        normmap = mcolors.Normalize(vmin=0, vmax=len(colors))
        coloring = dict((n, colormap(normmap(color)))
                          for n, color in coloring.iteritems())
    return coloring

