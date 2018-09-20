# *****************************************************************************
#       Copyright (C) 2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
from ..sweep import Sweep

import numpy as np
import networkx as nx


def test_unit_cube(unit_cube):
    s = Sweep(unit_cube.events)
    assert round(s.calculate_volume(), 8) == 1.


def test_simplex(simplex):
    s = Sweep(simplex.events)
    assert round(s.calculate_volume(), 8) == round(1. / np.math.factorial(simplex.dim), 8)


def test_non_regular_overlapping(cube_simplex_overlapping_3d):
    s = Sweep(cube_simplex_overlapping_3d.events)
    assert round(s.calculate_volume(), 8) == 1.


def test_simplex_cube_overlapping_2d(cube_simplex_overlapping_2d):
    s = Sweep(cube_simplex_overlapping_2d.events)
    assert round(s.calculate_volume(), 8) == 1.


def test_unbounded_non_regular(unbounded_non_regular):
    s = Sweep(unbounded_non_regular.events, sweep_plane=np.array([1, 0]))
    assert s.calculate_volume(1) == 1
    assert s.calculate_volume(2) == 4
    assert s.calculate_volume(400) == 400 ** 2


def test_simplex_in_cube(simplex_in_cube):
    s = Sweep(simplex_in_cube.events)
    assert round(s.calculate_volume(), 5) == 4


def test_precision(cube_simplex_overlapping_3d_imprecise):
    s = Sweep(cube_simplex_overlapping_3d_imprecise.events)
    assert round(s.calculate_volume(), 3) == round(1 + 8. / 6. - 0.5 ** 3, 3)


def test_overlapping_simplices(overlapping_simplices_2):
    s = Sweep(overlapping_simplices_2.events)
    assert round(s.calculate_volume(), 3) == round((8 ** 3) / 6 + 0.5, 3)


def test_nodes_and_edges(simplex):
    s = Sweep(simplex.events)
    nodes, edges = s.nodes_and_edges()
    assert len(nodes) == 1
    assert len(edges) == 0


def test_graph_one_connected_component(hole_2d):
    s = Sweep(hole_2d.events)
    graph = s.graph()
    assert len([cc for cc in nx.connected_components(graph.to_undirected())]) == 1


def test_graph_two_connected_components(simplex_cube_disconnected):
    s = Sweep(simplex_cube_disconnected.events)
    graph = s.graph()
    assert len([cc for cc in nx.connected_components(graph.to_undirected())]) == 2
