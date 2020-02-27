
import warnings

import networkx

class NotTestedWarning(Warning):
    """Emit a warning for untested code.

    Consider using the warn_untested function instead.

    Sample usage::
        import verkko.misc.testutil as testutil
        import warnings
        warnings.warn('untested', testutil.NotTestedWarning)
    """
    pass
def warn_untested(msg="This code is untested"):
    """Emit a warning for untested code.

    To use this, simply call this in an untested function.  When the
    function runs, it will emit a warning like this::

        verkko/graph/nxutil.py:63: NotTestedWarning: This code is untested
          import verkko.misc.testutil as testutil ; testutil.warn_untested()


    Arguments:
        message: str, warning message, default 'This code is untested'
            This is printed as the warning.  You could include extra
            information here.

    """
    warnings.warn(msg, NotTestedWarning, stacklevel=2)


def assert_isomorphic(g1, g2, msg=None):
    """Assertion function for networkx isomorphism.

    This is a thin wrapper for `networkx.is_isomorphic`, with
    additional printing of the difference in node/edge sets if the
    graphs are not isomorphic.  Note that this printing of differences
    only is useful for *equal* graphs, not *isomorphic* graphs."""
    if not networkx.is_isomorphic(g1, g2):
        msg_ = ["%r and %r not isomorphic"%(g1, g2),
                #"A: %s"%networkx.to_dict_of_lists(g1),
                #"B: %s"%networkx.to_dict_of_lists(g2),
                ]
        n1 = set(g1.nodes_iter())
        n2 = set(g2.nodes_iter())
        if n1 != n2:
            msg_.append("Nodes in 1 only: %s"%(n1-n2))
            msg_.append("Nodes in 2 only: %s"%(n2-n1))
        e1 = set(frozenset((a,b)) for a,b in g1.edges_iter())
        e2 = set(frozenset((a,b)) for a,b in g2.edges_iter())
        if e1 != e2:
            msg_.append("Edges in 1 only: %s"%(' '.join('(%s,%s)'%(a,b) for a,b in e1-e2)))
            msg_.append("Edges in 2 only: %s"%(' '.join('(%s,%s)'%(a,b) for a,b in e2-e1)))
        if msg: msg_.insert(0, msg)
        raise AssertionError('\n'.join(msg_))

def assert_edges_equal(g1, g2, msg=None):
    """Assert that two graphs have the same edge sets.

    This is a more strict condition than assert_isomorphic, since the
    node labels must be the same, and less strict than g1 == g2 since
    this does not consider edge attributes (such as weight).

    This function does not currently work on multigraphs, but does
    work on digraphs.  When comparing a digraph to a simple graph,
    *undirected* edges are taken.
    """
    if (isinstance(g1, (networkx.MultiGraph))
        or isinstance(g2, (networkx.MultiGraph))):
        raise RuntimeError("This function doesn't work on multigraphs yet.")

    # For directed graphs, we need ordered edge pairs.
    edgetype = frozenset
    if isinstance (g1, networkx.DiGraph) and isinstance(g2, networkx.DiGraph):
        edgetype = tuple

    # Create sets of all edges in the graphs and compare them.
    e1 = set(edgetype((a,b)) for a,b in g1.edges_iter())
    e2 = set(edgetype((a,b)) for a,b in g2.edges_iter())
    #print e1, e2
    if e1 != e2:
        msg_ = ['%r and %r do not have the same edges'%(g1, g2)]
        if e2 - e1:
            msg_.append("Edges not in the first:")
            for e in e2-e1:
                msg_.append("  %s"%(e, ))
        if e1 - e2:
            msg_.append("Edges not in the second:")
            for e in e1-e2:
                msg_.append("  %s"%(e, ))
        if msg: msg_.insert(0, msg)
        raise AssertionError('\n'.join(msg_))
