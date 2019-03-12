# *****************************************************************************
#       Copyright (C) 2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
from sweepvolume.cell_decomposition import Cell_Decomposition
from sweepvolume.sweep import Sweep
from sweepvolume.geometry import Hyperplane

import numpy as np


def test_halfspace_restriction():
    h0 = Hyperplane(np.array([0, 1]), -1)
    h1 = Hyperplane(np.array([0, 1]), 0)
    h2 = Hyperplane(np.array([1, 0]), 0)
    h3 = Hyperplane(np.array([1, 0]), -1)
    h4 = Hyperplane(np.array([1, 1]), 0)

    hyperplanes_final = (h0, h1, h2, h3, h4)
    polytope = {(0, -1), (1, 1), (2, 1)}
    polytope_final = {(0, -1), (1, 1), (2, 1), (3, -1), (4, 1)}
    cd = Cell_Decomposition(hyperplanes_final, [polytope_final])
    cd1 = Cell_Decomposition((h0, h1, h2), [polytope])
    cd1.restrict_to_halfspace(h3, -1)
    cd1.restrict_to_halfspace(h4, 1)
    assert set(cd1.events).issubset(set(cd.events))
    assert set(cd.events).issubset(set(cd1.events))
    assert np.isclose(Sweep(cd.events).calculate_volume(), Sweep(cd1.events).calculate_volume())


def test_halfspace_restriction_2():
    a_0 = Hyperplane(np.array([1, 0, ]), 0)
    a_1 = Hyperplane(np.array([0, 1]), 0)
    a_2 = Hyperplane(np.array([1, 1]), -1)
    a_3 = Hyperplane(np.array([1, 0]), -1)
    a_4 = Hyperplane(np.array([0, 1]), -1)
    hyperplanes_start = [a_0, a_1, a_2, a_3]
    hyperplanes_final = [a_0, a_1, a_2, a_3, a_4]

    simplex = {(0, 1), (1, 1), (2, -1)}
    unit_cube_start = {(0, 1), (1, 1), (3, -1)}
    unit_cube = {(0, 1), (1, 1), (3, -1), (4, -1)}

    cd = Cell_Decomposition(hyperplanes_final, [simplex, unit_cube])
    cd1 = Cell_Decomposition(hyperplanes_start, [simplex, unit_cube_start])
    cd1.restrict_to_halfspace(a_4, -1)

    assert set(cd.events).issubset(set(cd1.events))
    assert set(cd1.events).issubset(set(cd.events))
    assert round(Sweep(cd.events).calculate_volume(), 8) == round(
        Sweep(cd1.events).calculate_volume(), 8)


def test_halfspace_restriction_3(translated_triangles):
    a = Hyperplane(np.array([1., 0., ]), 0.)
    translated_triangles.restrict_to_halfspace(a, 1)
    assert (len(translated_triangles.events)) == 3
    s = Sweep(translated_triangles.events)
    assert round(s.calculate_volume(), 6) == 0.25


def test_halfspace_restriction_4():
    """
    Halfspace restriction can create new events even if hyperplane is already in cd.hyperplanes.
    Case: From two vertex cones that neutralized each other in original cd only one is valid in
    restricted cd.
    Here vertex at (0,1) is not part of original cd but of restricted cd.
    """
    a_0 = Hyperplane(np.array([1, 0, ]), 0)
    a_1 = Hyperplane(np.array([0, 1]), 0)
    a_2 = Hyperplane(np.array([1, 0]), -1)
    a_3 = Hyperplane(np.array([0, 1]), -1)

    b_0 = Hyperplane(np.array([1, 0, ]), 0.5)
    b_1 = Hyperplane(np.array([0, 1]), -0.5)
    b_2 = Hyperplane(np.array([1, 0]), -0.5)
    b_3 = Hyperplane(np.array([0, 1]), -1.5)

    hyperplanes = [a_0, a_1, a_2, a_3, b_0, b_1, b_2, b_3]
    cube_1 = {(0, 1), (1, 1), (2, -1), (3, -1)}
    cube_2 = {(4, 1), (5, 1), (6, -1), (7, -1)}

    cube_2_restricted = {(4, 1), (5, 1), (6, -1), (7, -1), (3, -1)}
    cd = Cell_Decomposition(hyperplanes, [cube_1, cube_2])
    cd.restrict_to_halfspace(a_3, -1)
    s = Sweep(cd.events)
    cd2 = Cell_Decomposition(hyperplanes, [cube_1, cube_2_restricted])
    s2 = Sweep(cd2.events)
    assert round(s2.calculate_volume(), 5) == 1.25
    assert round(s.calculate_volume(), 5) == 1.25


def test_both_halfspace_restriction(simplex_in_cube):
    '''
    test restricts cell decomposition to other halfspace of shared hyperplane of simplex and
    unit cube. All events lay on the hyperplane or outside of the restriction halfspace.
    No cone can be legal for events on the halfspace.
    '''
    cut_plane = Hyperplane(np.array([1, 0]), 1)
    simplex_in_cube.restrict_to_halfspace(cut_plane, -1)
    assert len(simplex_in_cube.events) == 0


def test_as_dict_json_compatible(hole_4d):
    """
    Test just dumps the cell decomposition and succeeds if no error is thrown.
    """
    import json
    print(json.dumps(hole_4d.as_dict(), sort_keys=True, indent=4, separators=(',', ': ')))