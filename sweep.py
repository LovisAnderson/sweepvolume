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
from event import Event

import networkx as nx
from graphviz import Digraph

import logging


class Sweep(object):

    def __init__(self, events, sweep_plane=None):
        if not events:
            logging.warning("--Sweep called with no vertices--")
            self.dimension = 0
            self.sorted_events = []
        else:
            e = events.pop()
            events.add(e)
            self.dimension = e.vertex.dim
            if sweep_plane is not None:
                self.sweep_plane = sweep_plane
            else:
                self.sweep_plane = self.random_sweep_plane(self.dimension)
            # tuple of (event, lambda)
            self.sorted_events = self.sort_events_by_lambda(events)
        self.gammas = []
        self.sweep()

    @staticmethod
    def random_sweep_plane(dimension):
        np.random.seed(0)
        random_orientation = np.random.choice([-1, 1], [dimension])
        return np.random.rand(dimension) * random_orientation

    def sort_events_by_lambda(self, events):
        dtype = [('event', Event), ('lambda', float)]
        events_with_lambda = list()

        for event in events:
            lam = np.dot(self.sweep_plane, event.vertex.coordinates)
            events_with_lambda.append((event, lam))
        events_with_lambda = np.array(events_with_lambda, dtype=dtype)
        sorted_by_lambda = np.sort(events_with_lambda, order='lambda')

        return sorted_by_lambda

    def sweep(self):
        """
        Computes all volume contributios (gammas) for events
        :return:
        """
        for event, lam in self.sorted_events:
            gamma = 0
            for cone in event.cones:
                gamma += self.calculate_gamma(cone)
            self.gammas.append(gamma)

    def calculate_gamma(self, cone):
        product = np.prod(np.dot(cone.rays.astype(float), self.sweep_plane))
        gamma = np.reciprocal(product) * abs(cone.determinant)
        return gamma

    def calculate_volume(self, lam=None):
        """"
        Calculates the sweep volume
        :param lam: gives together with the sweep-direction the halfspace up too which the volume is calculated
        :return: the volume
        """
        if len(self.sorted_events) == 0:
            return 0
        if lam is None:
            lam = self.sorted_events[-1][1]
        vol = 0
        for i in range(len(self.sorted_events)):
            lam_i = self.sorted_events[i][1]
            if lam < lam_i:
                break
            vol += self.gammas[i]* (lam - lam_i) ** self.dimension
        return np.reciprocal(float(np.math.factorial(self.dimension))) *vol

    def calculate_volumes(self, lambdas):
        f = []
        for i in range(len(self.sorted_events)):
            lam_i = self.sorted_events[i][1]
            f.append(np.piecewise(lambdas,
                                  [lambdas > lam_i],
                                  [lambda lambdas: self.gammas[i] * np.power((lambdas - lam_i), self.dimension)]))
        return np.reciprocal(float(np.math.factorial(self.dimension))) * np.sum(f, axis=0)

    def graphviz_graph(self, contract=True):

        dot = Digraph(comment='Sweep incidence Tree')
        nodes, edges = self.nodes_and_edges(contract=contract)
        for event in nodes:
            id = str(hash(event))
            node = dot.node(id, '{}'.format(event.incident_polytopes))
        for edge in edges:
            dot.edge(*[str(hash(edge[i])) for i in [0,1]])
        return dot

    def graph(self, contract=True):

        g = nx.DiGraph()
        nodes, edges = self.nodes_and_edges(contract=contract)
        g.add_nodes_from(nodes)
        g.add_edges_from(edges)

        return g

    def nodes_and_edges(self, contract=True):
        """
        Method creates lists of nodes, edges that build sweep-plane graph
        :param contract: Boolean that if sets to true contracts nodes that have same incidences as predecessor
        :return: nodes, edges: nodes is a list of vertices, edges a list of tuples (vertex, vertex)
        """
        polytope_last_event = {}

        nodes = set()
        edges = set()
        for event, lam in self.sorted_events:
            event_edges = []
            for incidence in event.incident_polytopes:
                last_node = polytope_last_event.get(incidence)
                if last_node:
                    if contract and set(event.incident_polytopes) == set(last_node.incident_polytopes):
                        continue
                    event_edges.append((last_node, event))
                nodes.add(event)
                polytope_last_event[incidence] = event
            edges = edges.union(event_edges)
        return nodes, edges