# *****************************************************************************
#       Copyright (C) 2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************

import logging
from copy import deepcopy

import numpy as np


class Event(object):
    def __init__(self, vertex, incidences=None):

        self.vertex = vertex
        self.incidences = incidences if incidences else set()

    def set_position_vector(self, hyperplanes):
        self.position_vector = np.zeros(len(hyperplanes))
        for index in range(len(hyperplanes)):
            if index in self.incidences:
                continue
            self.position_vector[index] = self.vertex.position(hyperplanes[index])
            if self.position_vector[index] == 0:
                self.incidences.add(index)

    def update_position_vector(self, hyperplanes):
        """
        :param hyperplanes: list of hyperplanes for which
         the updated position vector is to be calculated

        """

        for hyperplane in hyperplanes[len(self.position_vector):]:
            position = self.vertex.position(hyperplane)
            self.position_vector = np.append(self.position_vector, position)

    def polytope_incidences(self, polytope_vectors):
        incidences = []
        for polytope_index in range(len(polytope_vectors)):
            if self.vector_matches_position_vector(polytope_vectors[polytope_index]):
                incidences.append(polytope_index)
        return incidences

    def update_polytope_incidences(self, polytope_vectors):
        """
        Method recomputes polytope incidences if polytope vectors are expanded. Only checks
        if incidences are still true. For polytope vectors that have changed through halfspace
        restriction no new incidences can arise but old incidences can vanish.

        """
        nr_incidences_before = len(self.incident_polytopes)
        incidences = []
        for incidence in self.incident_polytopes:
            if self.vector_matches_position_vector(polytope_vectors[incidence]):
                incidences.append(incidence)
        self.incident_polytopes = incidences
        if len(self.incident_polytopes) > nr_incidences_before:
            logging.WARNING('number of incident polytopes higher than before')

    def vector_matches_position_vector(self, vector):
        """
        tests if input pos vector matches vertex position vector. For input vector not all
        position have to be specified.

        :param vector: list of tuples (hyperplane_index, orientation)
        :return: matches: Boolean that indicates if self.position_vector matches
                          position vector match.
        """
        matches = True
        for hyperplane_index, orientation in vector:
            if self.position_vector[hyperplane_index] == 0:
                continue
            elif self.position_vector[hyperplane_index] != orientation:
                return False
        return matches

    def __str__(self):
        return '[' + ','.join(['{:.3f}'.format(c) for c in self.vertex.coordinates]) + ']'

    def __hash__(self):
        return hash('<{}/{}>'.format(self.vertex.vector, self.vertex.denominator))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __deepcopy__(self, memo):
        # We override the copy function to make a halfway deepcopy
        # Some attributes are fixed and do not need a deepcopy
        fixed_attributes = ['vertex']
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            if k in fixed_attributes:
                setattr(result, k, v)
            else:
                setattr(result, k, deepcopy(v, memo))
        return result
