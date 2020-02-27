from nose.tools import *
import networkx as nx
from ..misc.testutil import assert_isomorphic
import nxutil

def test_from_string():
    g1 = nxutil.from_string('1-2')
    g2 = nx.Graph()
    g2.add_edge(1, 2)
    assert_isomorphic(g1, g2)

    g1 = nxutil.from_string('1-2-3')
    g2.add_edge(1, 2)
    g2.add_edge(2, 3)
    assert_isomorphic(g1, g2)

    g1 = nxutil.from_string('1-2-3-1')
    g2.add_edge(1, 2)
    g2.add_edge(2, 3)
    g2.add_edge(1, 3)
    assert_isomorphic(g1, g2)

    g1 = nxutil.from_string('1-2 3-4')
    g2 = nx.Graph()
    g2.add_edge(1, 2)
    g2.add_edge(3, 4)
    assert_isomorphic(g1, g2)

    g1 = nxutil.from_string('1-2-3 4-5  5-2')
    g2 = nx.Graph()
    g2.add_edges_from([(1,2), (2,3), (4,5), (5,2) ])
    assert_isomorphic(g1, g2)

