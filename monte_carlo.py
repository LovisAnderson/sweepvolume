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


class MonteCarloVolume(object):
    def __init__(self, cell_decomposition, iterations=10000):
        self.vertices = cell_decomposition.events
        self.dimension = cell_decomposition.dimension
        self.polytopes = cell_decomposition.polytope_vectors
        self.hyperplanes = cell_decomposition.hyperplanes
        self.bounding_box = self._bounding_box()
        self.iterations = iterations
        self.hits = self.random_drawing()
        self.volume = self._volume()

    def _bounding_box(self):
        min_vector = np.array([10000]*self.dimension)
        max_vector = np.zeros(self.dimension)
        for vertex in self.vertices:
            for i in range(self.dimension):
                if len(vertex.cones) > 0:
                    if vertex.coordinates[i] < min_vector[i]:
                        min_vector[i] = vertex.coordinates[i]
                    if vertex.coordinates[i] > max_vector[i]:
                        max_vector[i] = vertex.coordinates[i]
        return np.transpose(np.array([min_vector, max_vector]))

    def _is_legal_point(self, point):
        for polytope in self.polytopes:
            is_in_polytope = True
            for halfspace in polytope:
                scalar_product = np.dot(self.hyperplanes[halfspace[0]].a, point)
                orient = -1 if scalar_product + self.hyperplanes[halfspace[0]].b < 0 else 1
                if not orient == halfspace[1]:
                    is_in_polytope = False
            if is_in_polytope:
                return True
        return False

    def random_drawing(self):
        inside_body = 0
        for i in range(self.iterations):
            point = np.array([np.random.uniform(coordinate[0], coordinate[1]) for coordinate in self.bounding_box])
            if self._is_legal_point(point):
                inside_body += 1
        return inside_body

    def _volume(self):
        bb_volume = 1
        for coordinate in self.bounding_box:
            bb_volume = bb_volume * (abs(coordinate[1] - coordinate[0]))
        volume = bb_volume * (float(self.hits) / float(self.iterations))
        return volume