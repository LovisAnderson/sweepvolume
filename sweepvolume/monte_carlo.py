#*****************************************************************************
#       Copyright (C) 2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
import numpy as np
from .geometry import Polytope, Vertex


class MonteCarloVolume(object):
    def __init__(self, polytopes, iterations=10000):
        """
        Method to compute the volume with a monte carlo algorithm.
        :param polytopes: List of geometry.Polytope objects
        :param iterations: number of points drawn in monte carlo process
        """

        self.dim = polytopes[0].dim
        self.polytopes = polytopes
        self.bounding_box = self._bounding_box()
        self.iterations = iterations
        self.hits = self.random_drawing()
        self.volume = self._volume()

    def _bounding_box(self):
        min_vector = np.array([1e+14] * self.dim)
        max_vector = np.array([-1e+14] * self.dim)
        polytope_vertices = set().union(*[poly.vertices for poly in self.polytopes])
        from geometry import Polytope
        convex_hull = Polytope(vertices=polytope_vertices)
        for vertex in convex_hull.vertices:
            for i in range(self.dim):
                if vertex.coordinates[i] < min_vector[i]:
                    min_vector[i] = vertex.coordinates[i]
                if vertex.coordinates[i] > max_vector[i]:
                    max_vector[i] = vertex.coordinates[i]
        return np.transpose(np.array([min_vector, max_vector]))

    def _is_legal_point(self, point):
        for polytope in self.polytopes:
            is_in_polytope = polytope.point_in_poly(point)
            if is_in_polytope:
                return True
        return False

    def random_drawing(self):
        inside_body = 0
        for i in range(self.iterations):
            point = np.array([np.random.uniform(coordinate[0], coordinate[1])
                              for coordinate in self.bounding_box])
            v = Vertex.vertex_from_coordinates(point)
            if self._is_legal_point(v):
                inside_body += 1
        return inside_body

    def _volume(self):
        bb_volume = 1
        for coordinate in self.bounding_box:
            bb_volume = bb_volume * (abs(coordinate[1] - coordinate[0]))
        volume = bb_volume * (float(self.hits) / float(self.iterations))
        return volume


def test_monte_carlo():
    from geometry import Hyperplane, Polytope
    a1 = np.array([1, 0])
    a2 = np.array([-1, 0])

    a3 = np.array([0, 1])
    a4 = np.array([0, -1])

    p1 = Hyperplane(a1, 0)
    p2 = Hyperplane(a2, -1)
    p3 = Hyperplane(a3, -1)
    p4 = Hyperplane(a4, 0)
    unit_cube_1 = Polytope(halfspaces=zip([p1, p2, p3, p4], [-1]*4))

    q1 = Hyperplane(a1, -1)
    q2 = Hyperplane(a2, 0)
    q3 = Hyperplane(a3, 0)
    q4 = Hyperplane(a4, -1)
    unit_cube_2 = Polytope(halfspaces=zip([q1, q2, q3, q4], [-1]*4))
    m = MonteCarloVolume([unit_cube_1, unit_cube_2], iterations=100000)
    error = np.abs(m.volume - 2)
    assert error < 0.1
