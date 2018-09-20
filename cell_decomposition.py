# *****************************************************************************
#       Copyright (C) 2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
import copy
import itertools
import logging

import json
import os

import numpy as np
import sympy

from geometry import Hyperplane, Polytope, Cone

from event import Event


class Cell_Decomposition(object):
    def __init__(self,
                 hyperplanes,
                 polytope_vectors,
                 bounding_box=None):
        """
         Constructor for cell_decomposition object
        :param hyperplanes: List of hyperplane objects
        :param polytope_vectors: List of list of tuples (index,position), whereby index specifies
                                 index of hyperplane in hyperplanes and position specifies
                                 halfspace in which polytope lies.
        :param bounding_box: Bounding box in which all polytopes are included.
                             Given as list (length dimension) of lists (length 2) which give lower
                             bound upper bound for each dimension.
        """
        """
       
        :param hyperplanes:
        :param polytope_vectors: 
        """
        self.polytope_vectors = polytope_vectors
        self.hyperplanes = np.array(hyperplanes)
        self.nr_of_hyperplanes = len(hyperplanes)
        self.bbox = bounding_box
        if len(hyperplanes) > 0:
            self.dimension = hyperplanes[0].dimension
            self.possible_events = set(self._find_all_possible_events())
            self.events = self.find_events()
        else:
            self.possible_events = {}
            self.events = []
            self.dimension = 0

    def _find_all_possible_events(self):

        """
        Main method for calculating vertices and their incidences in the decomposition.
        :return: vertices, a list of vertices with incidences saved
        """

        events = dict()
        for combination in itertools.combinations(range(self.nr_of_hyperplanes), self.dimension):
            vertex = self.solve_les(combination)
            if vertex is None:
                continue
            if not self.inside_bbox(vertex):
                continue
            if events.get(hash(vertex)):
                e = events[hash(vertex)]
                e.incidences.update(combination)
            else:
                e = Event(vertex)
                e.incidences.update(combination)
                e.set_position_vector(self.hyperplanes)
                events[hash(vertex)] = e
        return events.values()

    def inside_bbox(self, vertex):
        if not self.bbox:
            return True
        eps = 1e-2
        for i, coordinate in enumerate(vertex.coordinates):
            if not (self.bbox[0][i] - eps < coordinate < self.bbox[1][i] + eps):
                return False
        return True

    def solve_les(self, hyperplane_indices):
        hyperplanes = self.hyperplanes[np.array(hyperplane_indices)]
        p = Polytope(zip(hyperplanes, [0] * self.dimension))
        if len(p.vertices) == 1:
            return p.vertices[0]
        else:
            return None

    def find_events(self):

        events = set()
        for event in self.possible_events:
            event.incident_polytopes = event.polytope_incidences(self.polytope_vectors)
            event.cones = self.vertex_cones(event) if len(event.incident_polytopes) > 0 else []
            if len(event.cones) > 0:
                events.add(event)
        return events

    def is_legal_cone(self, cone_position_vector, event):
        """

        :param cone_position_vector: quasi position vector of the cone. if it is equal to the event position vector for
            every index at which the vpv is not equal zero and equal to the cone psv at the other entries
        :param event: event at which cone is supported
        :return: boolean indicating if cone is in conic decomposition
        """

        def position_vector_match(polytope_position_vector, cone_position_vector):
            for index, position in polytope_position_vector:
                if cone_position_vector[index] != position:
                    return False
            return True

        for poly in event.incident_polytopes:
            if position_vector_match(self.polytope_vectors[poly], cone_position_vector):
                return True
        return False

    def vertex_cones(self, event):
        """
        Method computes the active cones for a vertex in regularized form.
        :param event: Vertex object
        :return: List of cones in regularized form
        """
        cones = []
        incidences = list(event.incidences)
        cone_pos_vector = copy.deepcopy(event.position_vector)
        for combination in itertools.product([-1, 1], repeat=len(event.incidences)):
            cone = []
            incidence_dict = {}
            for i in range(len(combination)):
                incidence_dict[incidences[i]] = combination[i]
                cone_pos_vector[incidences[i]] = combination[i]
                cone.append((incidences[i], combination[i]))
            if self.is_legal_cone(cone_pos_vector, event):
                if len(event.incidences) > self.dimension:
                    # Test if cone is empty, if so we can skip it
                    halfspaces = [(self.hyperplanes[index], orientation)
                                  for index, orientation in cone]
                    poly = Polytope(halfspaces=halfspaces)
                    if len(poly.poly.minimized_generators()) <= 1:
                        continue
                    cones += self.regularize_cone(cone, event.vertex)
                else:
                    cone = [(self.hyperplanes[halfspace[0]], halfspace[1]) for halfspace in cone]
                    cones.append(Cone(cone, incidences=incidence_dict))
        # remove neutralizing vertex cones only when all cones are regular.
        if self.dimension == len(event.incidences):
            cones = self.remove_neutralizing_vertex_cones(cones)
        return cones

    @staticmethod
    def remove_neutralizing_vertex_cones(cones):
        """
        If a vertex has two cones whose position vectors differ at exactly one entry their
        contribution to the volume is together 0. We therefor remove them from the list of cones.
        Typical case in 2d is  -> _|_ <-where the vertex is on the crossing and the two cones are
        left and right.
        :param cones:
        :return: list without cancelled out vertex cones
        """
        contributing_cones = []
        while cones:
            cone1 = cones.pop()
            contributing = True
            for cone2 in cones:
                if cone1.is_neutralizing_cone(cone2):
                    cones.remove(cone2)
                    contributing = False
                    break
            if contributing:
                contributing_cones.append(cone1)
        return contributing_cones

    def regularize_cone(self, cone, vertex):

        """"
        Method calculates vertex cones in regularized form.
        :param cone: pos. vec. of cone
        :param vertex:
        :return:
        """

        ordered_cone = self.first_rows_independent(cone)
        hyperplanes, pos_vector = self._construct_shift_polytope(ordered_cone, vertex)
        cd = Cell_Decomposition(hyperplanes, [pos_vector])
        cones = [c for event in cd.events for c in event.cones]
        return cones

    def _construct_shift_polytope(self, cone, vertex):

        """
        Method calculates regularization polytope
        :param cone: cone as pos.vec. s.t. first d hyperplanes have l.i. normal vectors
        :param vertex: vertex of cone
        :return: hyperplanes, position vector corrsponding to shift polytope
        """

        polytope_vertices = [vertex]
        polytope_halfspaces = [(self.hyperplanes[index], orientation)
                               for index, orientation in cone[:self.dimension]]
        for hyperplane_index, orientation in cone[self.dimension:]:
            halfspace = self.compute_halfspace_shift(
                [self.hyperplanes[hyperplane_index], orientation], polytope_vertices)
            polytope_halfspaces.append(halfspace)
            polytope_vertices = Polytope(polytope_halfspaces).vertices
        hyperplanes = [hyperplane for hyperplane, orientation in polytope_halfspaces]
        orientations = [orientation for hyperplane, orientation in polytope_halfspaces]
        return hyperplanes, zip(np.arange(len(hyperplanes)), orientations)

    def first_rows_independent(self, halfspaces):
        """
        Method orders halfspaces s.t. the first d have independent normals
        :param halfspaces: halfspaces that correspond to cone as list of tuples
                          (halfspace_index, orientation)
        :return: halfspaces ordered s.t. the first d have independent normals
        """
        first = []
        last = []
        mat = np.array([self.hyperplanes[hyperplane_ind].a for hyperplane_ind, _ in halfspaces])
        _, inds = sympy.Matrix(mat).T.rref()
        if len(inds) != self.dimension:
            logging.warning('Not d many linear independent normals in regularization step')
        for index in range(len(halfspaces)):
            if index in inds:
                first.append(halfspaces[index])
            else:
                last.append(halfspaces[index])
        return first + last

    @staticmethod
    def compute_halfspace_shift(halfspace, vertices):
        """
        Method computes shifted hyperplane for regularization. New halfspaces contains all vertices
         in its interior, but has the same normal vector as input halfspace
        :param halfspace: A tuple (Hyperplane, orientation)
        :param vertices: list of vertices to be contained in open halfspace
        :return: Tuple of (Hyperplane, orientation)
        """
        shifting_distance = round(np.linalg.norm(halfspace[0].a))
        a = halfspace[0].a
        min_b = max_b = - halfspace[0].b
        for vertex in vertices:
            shift = np.dot(a, vertex.coordinates)
            min_b = min(min_b, shift)
            max_b = max(max_b, shift)
        if halfspace[1] == -1:
            b = round(max_b + shifting_distance)
        else:
            b = round(min_b - shifting_distance)
        hyperplane = Hyperplane(a, -b)
        hyperplane.integer = True
        return (hyperplane, halfspace[1])

    def restrict_to_halfspace(self, hyperplane, orientation):
        """
        Method is called with halfspace and restricts cell decomposition to that halfspace.
        :param hyperplane: Hyperplane object
        :param orientation:
        :return:
        """
        index = np.where(self.hyperplanes == hyperplane)
        hyperplane_is_new = True
        if len(index[0]) == 0:
            self.hyperplanes = np.append(self.hyperplanes, hyperplane)
            index = len(self.hyperplanes) - 1
            self.nr_of_hyperplanes += 1
        else:
            index = index[0][0]
            hyperplane_is_new = False
        self._update_polytopes((index, orientation), hyperplane_is_new)
        for event in list(self.events):
            self._restrict_event_to_halfspace(event, (index, orientation), hyperplane_is_new)
        if hyperplane_is_new:
            self._find_new_events((index, orientation))
        else:
            self._update_vertex_cones((index, orientation))

    def _update_polytopes(self, halfspace, hyperplane_is_new):
        # todo: We might want to remove polytope if it is not full dimensional anymore
        for poly in self.polytope_vectors:
            poly.add(halfspace)

    def _update_vertex_cones(self, halfspace):
        # events might not have been relevant before due to neutralizing vertex cones
        # for all possible events vertex cones have to be recalculated
        for event in self.possible_events:
            if halfspace[0] not in event.incidences:
                continue
            else:
                event.update_position_vector(self.hyperplanes)
                event.incident_polytopes = event.polytope_incidences(self.polytope_vectors)
                if len(event.incident_polytopes) > 0:
                    event.cones = self.remove_neutralizing_vertex_cones(self.vertex_cones(event))
                    if event.cones:
                        self.events.add(event)

    def _find_new_events(self, halfspace):
        """
        Method computes new events that emerge from halfspace intersection.
        :param halfspace: Tuple of (hyperplane_index, orientation)
        :return:
        """
        # iterate over all combinations of d hyperplanes that contain the given hyperplane
        hyperplane_indices = set(range(self.nr_of_hyperplanes))
        hyperplane_indices.remove(halfspace[0])
        for combination in itertools.combinations(hyperplane_indices, self.dimension - 1):
            combination = np.append(np.array(combination), halfspace[0])
            x = self.solve_les(combination)
            if x is None or not self.inside_bbox(x):
                continue
            e = Event(x)
            if e not in self.events:

                e.incidences.update(combination)
                e.set_position_vector(self.hyperplanes)
                e.incident_polytopes = e.polytope_incidences(self.polytope_vectors)
                if len(e.incident_polytopes) > 0:
                    e.cones = self.remove_neutralizing_vertex_cones(self.vertex_cones(e))
                    if e.cones:
                        self.events.add(e)
                    self.possible_events.add(e)
        return

    def _restrict_event_to_halfspace(self, event, halfspace, hyperplane_is_new):
        """
        Method is refreshing event for new halfspace. Event is removed if it lies on the other side
        of the hyperplane. Incidences and cones are actualized if it lies on the hyperplane and
        hyperplane is new.
        :param event: An Event
        :param halfspace: New halfspace as tuple (hyperplane_index, orientation)
        :param hyperplane_is_new: Boolean
        :return:
        """
        hyperplane = self.hyperplanes[halfspace[0]]
        position = event.vertex.position(hyperplane)
        if hyperplane_is_new:
            event.position_vector = np.append(event.position_vector, position)
        event.update_polytope_incidences(self.polytope_vectors)
        if halfspace[1] != position and position != 0:
            self.events.remove(event)
            return
        if not hyperplane_is_new and position == 0:
            event.update_position_vector(self.hyperplanes)
            event.cones = self.vertex_cones(event)
            if len(event.cones) == 0:
                self.events.remove(event)
            return
        if position == 0:
            event.incidences.update([halfspace[0]])
            event.cones = self.vertex_cones(event)
            if len(event.cones) == 0:
                self.events.remove(event)
        return

    def as_dict(self, initialized=True):
        hyperplanes = map(lambda hyp: repr(np.append(hyp.a, hyp.b)), self.hyperplanes)
        events = [str(event) for event in self.events] if initialized else []
        cd_dict = {
            'hyperplanes': hyperplanes,
            'events': events,
            'polytope_vectors': (map(lambda p: list(p), self.polytope_vectors)),
            'dimension': self.dimension
        }
        return cd_dict

    def dump_json(self, output_dir='../../data/cds', fn='cell_decomposition.json'):

        with open(os.path.join(output_dir, fn), 'w') as file:
            json.dump(self.as_dict(), file, indent=4)
        return

    def remove_outside_vertices(self):
        self.possible_events = set(filter(
            self.inside_bbox,
            self.possible_events)
        )

    def __deepcopy__(self, memo):
        # We override the copy function to make a halfway deepcopy
        # Some attributes are fixed and do not need a deepcopy
        from copy import copy, deepcopy
        fixed_attributes = ['dimension', 'bbox']
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            if k in fixed_attributes:
                setattr(result, k, v)
            else:
                setattr(result, k, deepcopy(v, memo))
        return result


def cd_from_json(path):
    with open(path) as file:
        cd_dict = json.load(file)
    hyperplanes = np.array(map(lambda hyp: Hyperplane(np.array(hyp[0]), float(hyp[1])),
                               cd_dict['hyperplanes']))
    polytopes = cd_dict['polytope_vectors']
    return Cell_Decomposition(hyperplanes, polytopes)
