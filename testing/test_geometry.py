# *****************************************************************************
#       Copyright (C) 2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
import numpy as np


def test_position():
    from ..geometry import Vertex, Hyperplane
    v1 = Vertex.vertex_from_coordinates(np.array([1., 1., 1.]) / 2.)
    v2 = Vertex.vertex_from_coordinates(np.array([1., 1., 0.]) / 2.)
    v3 = Vertex.vertex_from_coordinates(np.array([1., 1., 1.]) / 4.)
    h = Hyperplane(np.array([1., 1., 1.]), -1.)
    assert v1.position(h) == 1
    assert v2.position(h) == 0
    assert v3.position(h) == -1


def test_vector_distance():
    from ..geometry import vector_distance
    v1 = np.array([1, 0])
    v2 = np.array([0, 1])
    v3 = np.array([1, 1])
    assert vector_distance(v1, v2) == 1
    assert vector_distance(v1, v1) == 0
    assert np.isclose(vector_distance(v1, v3), 1 - np.reciprocal(np.sqrt(2)))
