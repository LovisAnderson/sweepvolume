# *****************************************************************************
#       Copyright (C) 2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
import pytest
from ..geometry import Hyperplane
import numpy as np
from ..cell_decomposition import Cell_Decomposition


@pytest.fixture
def paper_example():
    f_0 = Hyperplane(np.array([6, 11]), -100)
    f_1 = Hyperplane(np.array([2, 1]), -28)
    f_2 = Hyperplane(np.array([2, -3]), 4)
    f_3 = Hyperplane(np.array([-2, 7]), -52)

    hyperplanes = [f_0, f_1, f_2, f_3]

    p_0 = set([(0, 1), (1, -1), (3, -1)])
    p_1 = set([(0, 1), (2, -1), (3, -1)])

    return Cell_Decomposition(hyperplanes, [p_0, p_1])


@pytest.fixture
def simplex():
    a_0 = Hyperplane(np.array([1, 0, 0]), 0)
    a_1 = Hyperplane(np.array([0, 1, 0]), 0)
    a_2 = Hyperplane(np.array([0, 0, 1]), 0)
    a_3 = Hyperplane(np.array([-1, -1, -1]), 1)

    hyperplanes = [a_0, a_1, a_2, a_3]

    simplex = {(0, 1), (1, 1), (2, 1), (3, 1)}
    return Cell_Decomposition(hyperplanes, [simplex])


@pytest.fixture
def unit_cube():
    a_0 = Hyperplane(np.array([1, 0, 0]), 0)
    a_1 = Hyperplane(np.array([0, 1, 0]), 0)
    a_2 = Hyperplane(np.array([0, 0, 1]), 0)
    a_3 = Hyperplane(np.array([1, 0, 0]), -1)
    a_4 = Hyperplane(np.array([0, 1, 0]), -1)
    a_5 = Hyperplane(np.array([0, 0, 1]), -1)

    hyperplanes = [a_0, a_1, a_2, a_3, a_4, a_5]

    unit_cube = {(0, 1), (1, 1), (2, 1), (3, -1), (4, -1), (5, -1)}
    return Cell_Decomposition(hyperplanes, [unit_cube])


@pytest.fixture
def simplex_in_cube_3d():
    a_0 = Hyperplane(np.array([1, 0, 0]), 0)
    a_1 = Hyperplane(np.array([0, 1, 0]), 0)
    a_2 = Hyperplane(np.array([0, 0, 1]), 0)
    a_3 = Hyperplane(np.array([-1, -1, -1]), 1)
    a_4 = Hyperplane(np.array([1, 0, 0]), -1)
    a_5 = Hyperplane(np.array([0, 1, 0]), -1)
    a_6 = Hyperplane(np.array([0, 0, 1]), -1)

    hyperplanes = [a_0, a_1, a_2, a_3, a_4, a_5, a_6]

    simplex = {(0, 1), (1, 1), (2, 1), (3, 1)}
    unit_cube = {(0, 1), (1, 1), (2, 1), (4, -1), (5, -1), (6, -1)}
    return Cell_Decomposition(hyperplanes, [unit_cube, simplex])


@pytest.fixture
def cube_simplex_overlapping_3d():
    c_0 = Hyperplane(np.array([1, 0, 0]), 0)
    c_1 = Hyperplane(np.array([0, 1, 0]), 0)
    c_2 = Hyperplane(np.array([0, 0, 1]), 0)
    c_3 = Hyperplane(np.array([1, 0, 0]), -1)
    c_4 = Hyperplane(np.array([0, 1, 0]), -1)
    c_5 = Hyperplane(np.array([0, 0, 1]), -1)

    s_0 = Hyperplane(np.array([1, 1, 1]), -2.)
    s_1 = Hyperplane(np.array([1, 0, 0]), -0.5)
    s_2 = Hyperplane(np.array([0, 1, 0]), -0.5)
    s_3 = Hyperplane(np.array([0, 0, 1]), -0.5)

    hyperplanes = [c_0, c_1, c_2, c_3, c_4, c_5, s_0, s_1, s_2, s_3]

    simplex = {(6, -1), (7, 1), (8, 1), (9, 1)}
    unit_cube = {(0, 1), (1, 1), (2, 1), (3, -1), (4, -1), (5, -1)}
    return Cell_Decomposition(hyperplanes, [unit_cube, simplex])


@pytest.fixture
def translated_triangles():
    # geometry: <\<\ two congruent triangles that meet in 0
    h0 = Hyperplane(np.array([0, 1]), 0)
    h1 = Hyperplane(np.array([-1, 1]), -1)
    h2 = Hyperplane(np.array([1, 1]), 0)
    h3 = Hyperplane(np.array([-1, 1]), 0)
    h4 = Hyperplane(np.array([1, 1]), -1)

    hyperplanes = [h0, h1, h2, h3, h4]
    poly1 = {(0, 1), (1, -1), (2, -1)}
    poly2 = {(0, 1), (3, -1), (4, -1)}
    return Cell_Decomposition(hyperplanes, [poly1, poly2])


@pytest.fixture
def cube_simplex_overlapping_3d_imprecise():
    k = 1e+12
    c_0 = Hyperplane(np.array([k, 0, 0]), 0)
    c_1 = Hyperplane(np.array([0, k, 0]), 0)
    c_2 = Hyperplane(np.array([0, 0, k]), 0)
    c_3 = Hyperplane(np.array([k, 0, 0]), -k)
    c_4 = Hyperplane(np.array([0, k, 0]), -k)
    c_5 = Hyperplane(np.array([0, 0, k]), -k)

    s_0 = Hyperplane(np.array([k, k, k]), -3.5 * k)
    s_1 = Hyperplane(np.array([k, 0, 0]), -0.5 * k)
    s_2 = Hyperplane(np.array([0, k, 0]), -0.5 * k)
    s_3 = Hyperplane(np.array([0, 0, k]), -0.5 * k)

    hyperplanes = [c_0, c_1, c_2, c_3, c_4, c_5, s_0, s_1, s_2, s_3]

    simplex = {(6, -1), (7, 1), (8, 1), (9, 1)}
    unit_cube = {(0, 1), (1, 1), (2, 1), (3, -1), (4, -1), (5, -1)}
    return Cell_Decomposition(hyperplanes, [unit_cube, simplex])


@pytest.fixture
def cube_simplex_overlapping_2d():
    c_0 = Hyperplane(np.array([1, 0]), 0)
    c_1 = Hyperplane(np.array([0, 1]), 0)
    c_2 = Hyperplane(np.array([1, 0]), -1)
    c_3 = Hyperplane(np.array([0, 1]), -1)

    s_0 = Hyperplane(np.array([1, 0]), -0.5)
    s_1 = Hyperplane(np.array([0, 1]), -0.5)
    s_2 = Hyperplane(np.array([1, 1]), -1.5)
    hyperplanes = [c_0, c_1, c_2, c_3, s_0, s_1, s_2]

    simplex = {(4, 1), (5, 1), (6, -1)}
    unit_cube = {(0, 1), (1, 1), (2, -1), (3, -1)}
    return Cell_Decomposition(hyperplanes, [unit_cube, simplex])


@pytest.fixture
def unbounded_non_regular():
    a_0 = Hyperplane(np.array([-1, 1]), 0)
    a_1 = Hyperplane(np.array([1, 1]), 0)
    a_2 = Hyperplane(np.array([-(1 / 2), -1]), 0)
    a_3 = Hyperplane(np.array([1 / 2, -1]), 0)

    hyperplanes = [a_0, a_1, a_2, a_3]

    left_poly = {(1, -1), (2, -1)}
    right_poly = {(0, -1), (3, -1)}
    return Cell_Decomposition(hyperplanes, [left_poly, right_poly])


@pytest.fixture
def simplex_in_cube():
    s_0 = Hyperplane(np.array([1, 0]), 0)
    s_1 = Hyperplane(np.array([0, 1]), 0)
    s_2 = Hyperplane(np.array([-1, -1]), 1)

    c_0 = Hyperplane(np.array([1, 0]), 1)
    c_1 = Hyperplane(np.array([0, 1]), 1)
    c_2 = Hyperplane(np.array([1, 0]), -1)
    c_3 = Hyperplane(np.array([0, 1]), -1)
    hyperplanes = [s_0, s_1, s_2, c_0, c_1, c_2, c_3]
    simplex = {(0, 1), (1, 1), (2, 1)}
    cube = {(3, 1), (4, 1), (5, -1), (6, -1)}
    return Cell_Decomposition(hyperplanes, [simplex, cube])


@pytest.fixture
def non_regular_pyramid_3d():
    f_0 = Hyperplane(np.array([-1, 0, -1]), 0)
    f_1 = Hyperplane(np.array([0, -1, -1]), 0)
    f_2 = Hyperplane(np.array([0, 1, -1]), 0)
    f_3 = Hyperplane(np.array([1, 0, -1]), 0)
    f_4 = Hyperplane(np.array([0, 0, 1]), -1)

    hyperplanes = [f_0, f_1, f_2, f_3, f_4]
    p_0 = {(0, -1), (1, -1), (2, -1), (3, -1), (4, -1)}
    return Cell_Decomposition(hyperplanes, [p_0])


@pytest.fixture
def overlapping_simplices():
    s_0 = Hyperplane(np.array([1, 0, 0]), 0)
    s_1 = Hyperplane(np.array([0, 1, 0]), 0)
    s_2 = Hyperplane(np.array([0, 0, 1]), 0)
    s_3 = Hyperplane(np.array([1, 1, 1]), -4)

    s_4 = Hyperplane(np.array([0, 0, 1]), -1)
    s_5 = Hyperplane(np.array([0, 1, 0]), -1)
    s_6 = Hyperplane(np.array([1, 0, 0]), -1)
    s_7 = Hyperplane(np.array([1, 1, -1]), 0)

    hyperplanes = [s_0, s_1, s_2, s_3, s_4, s_5, s_6, s_7]
    simplex_1 = {(0, 1), (1, 1), (2, 1), (3, -1)}
    simplex_2 = {(4, -1), (5, 1), (6, 1), (7, -1)}

    return Cell_Decomposition(hyperplanes, [simplex_1, simplex_2])


@pytest.fixture
def overlapping_simplices_2():
    s_0 = Hyperplane(np.array([1, 0, 0]), 0)
    s_1 = Hyperplane(np.array([0, 1, 0]), 0)
    s_2 = Hyperplane(np.array([0, 0, 1]), 0)
    s_3 = Hyperplane(np.array([1, 1, 1]), -8)

    s_4 = Hyperplane(np.array([0, 0, 1]), -1)
    s_5 = Hyperplane(np.array([0, 1, 0]), -1)
    s_6 = Hyperplane(np.array([1, 0, 0]), -1)
    s_7 = Hyperplane(np.array([1, 1, -1]), -3)

    hyperplanes = [s_0, s_1, s_2, s_3, s_4, s_5, s_6, s_7]
    simplex_1 = {(0, 1), (1, 1), (2, 1), (3, -1)}
    simplex_2 = {(4, -1), (5, 1), (6, 1), (7, -1)}
    polytopes = [simplex_1, simplex_2]
    return Cell_Decomposition(hyperplanes, polytopes)


@pytest.fixture
def hole_2d():
    a1 = np.array([1, 0])
    a2 = np.array([-1, 0])

    a3 = np.array([0, 1])
    a4 = np.array([0, -1])

    p1 = Hyperplane(a1, 0)
    p2 = Hyperplane(a2, -1)
    p3 = Hyperplane(a3, -1)
    p4 = Hyperplane(a4, 0)
    P1 = set(zip(range(4), [-1] * 4))

    q1 = Hyperplane(a1, -1)
    q2 = Hyperplane(a2, 0)
    q3 = Hyperplane(a3, 0)
    q4 = Hyperplane(a4, -1)
    P2 = set(zip(range(4, 8), [-1] * 4))

    r1 = Hyperplane(a1, -2)
    r2 = Hyperplane(a2, 1)
    r3 = Hyperplane(a3, -1)
    r4 = Hyperplane(a4, 0)
    P3 = set(zip(range(8, 12), [-1] * 4))

    s1 = Hyperplane(a1, -1)
    s2 = Hyperplane(a2, 0)
    s3 = Hyperplane(a3, -2)
    s4 = Hyperplane(a4, 1)
    P4 = set(zip(range(12, 16), [-1] * 4))

    hyperplanes = [p1, p2, p3, p4, q1, q2, q3, q4, r1, r2, r3, r4, s1, s2, s3, s4]
    return (Cell_Decomposition(hyperplanes, [P1, P2, P3, P4]))


@pytest.fixture
def hole_3d():
    a1 = np.array([1, 0, 0])
    a2 = np.array([-1, 0, 0])

    a3 = np.array([0, 1, 0])
    a4 = np.array([0, -1, 0])

    a5 = np.array([0, 0, 1])
    a6 = np.array([0, 0, -1])

    p0 = Hyperplane(a5, -1)
    p1 = Hyperplane(a6, 0)

    p2 = Hyperplane(a1, 0)
    p3 = Hyperplane(a2, -1)
    p4 = Hyperplane(a3, -1)
    p5 = Hyperplane(a4, 0)
    P1 = set(zip(range(6), [-1] * 6))

    q1 = Hyperplane(a1, -1)
    q2 = Hyperplane(a2, 0)
    q3 = Hyperplane(a3, 0)
    q4 = Hyperplane(a4, -1)
    P2 = set(zip(range(2), [-1] * 2) + zip(range(6, 10), [-1] * 4))

    r1 = Hyperplane(a1, -2)
    r2 = Hyperplane(a2, 1)
    r3 = Hyperplane(a3, -1)
    r4 = Hyperplane(a4, 0)
    P3 = set(zip(range(2), [-1] * 2) + zip(range(10, 14), [-1] * 4))

    s1 = Hyperplane(a1, -1)
    s2 = Hyperplane(a2, 0)
    s3 = Hyperplane(a3, -2)
    s4 = Hyperplane(a4, 1)
    P4 = set(zip(range(2), [-1] * 2) + zip(range(14, 18), [-1] * 4))
    hyperplanes = [p0, p1, p2, p3, p4, p5, q1, q2, q3, q4, r1, r2, r3, r4, s1, s2, s3, s4]

    return Cell_Decomposition(hyperplanes, [P1, P2, P3, P4])


@pytest.fixture
def hole_4d():
    a1 = np.array([1, 0, 0, 0])
    a2 = np.array([-1, 0, 0, 0])

    a3 = np.array([0, 1, 0, 0])
    a4 = np.array([0, -1, 0, 0])

    a5 = np.array([0, 0, 1, 0])
    a6 = np.array([0, 0, -1, 0])

    a7 = np.array([0, 0, 0, 1])
    a8 = np.array([0, 0, 0, -1])

    p0 = Hyperplane(a5, -1)
    p1 = Hyperplane(a6, 0)
    p2 = Hyperplane(a7, -1)
    p3 = Hyperplane(a8, 0)

    p4 = Hyperplane(a1, 0)
    p5 = Hyperplane(a2, -1)
    p6 = Hyperplane(a3, -1)
    p7 = Hyperplane(a4, 0)
    P1 = set(zip(range(8), [-1] * 8))

    q1 = Hyperplane(a1, -1)
    q2 = Hyperplane(a2, 0)
    q3 = Hyperplane(a3, 0)
    q4 = Hyperplane(a4, -1)
    P2 = set(zip(range(4), [-1] * 4) + zip(range(8, 12), [-1] * 4))

    r1 = Hyperplane(a1, -2)
    r2 = Hyperplane(a2, 1)
    r3 = Hyperplane(a3, -1)
    r4 = Hyperplane(a4, 0)
    P3 = set(zip(range(4), [-1] * 4) + zip(range(12, 16), [-1] * 4))

    s1 = Hyperplane(a1, -1)
    s2 = Hyperplane(a2, 0)
    s3 = Hyperplane(a3, -2)
    s4 = Hyperplane(a4, 1)
    P4 = set(zip(range(4), [-1] * 4) + zip(range(16, 20), [-1] * 4))
    hyperplanes = [p0, p1, p2, p3, p4, p5, p6, p7, q1, q2, q3, q4, r1, r2, r3, r4, s1, s2, s3, s4]

    return Cell_Decomposition(hyperplanes, [P1, P2, P3, P4])


@pytest.fixture
def simplex_cube_disconnected():
    s_0 = Hyperplane(np.array([1, 0]), 0)
    s_1 = Hyperplane(np.array([0, 1]), 0)
    s_2 = Hyperplane(np.array([1, 1]), -1)

    c_0 = Hyperplane(np.array([1, 0]), -2)
    c_1 = Hyperplane(np.array([0, 1]), -2)
    c_2 = Hyperplane(np.array([1, 0]), -3)
    c_3 = Hyperplane(np.array([0, 1]), -3)
    hyperplanes = [s_0, s_1, s_2, c_0, c_1, c_2, c_3]
    simplex = {(0, 1), (1, 1), (2, -1)}
    cube = {(3, 1), (4, 1), (5, -1), (6, -1)}
    return Cell_Decomposition(hyperplanes, [simplex, cube])
