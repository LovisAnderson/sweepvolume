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

import numpy as np
from ppl import (Variable, point,
                 Constraint_System, Generator_System,
                 C_Polyhedron, Poly_Con_Relation, Poly_Gen_Relation)
from gmpy2 import mpz, gcd, divexact
from scipy.spatial import ConvexHull


class Hyperplane(object):
    """
    Class for hyperplane object.
    Representation: a1*x1 + ... + ad*xd + b = 0
    """

    def __init__(self, a, b):
        assert isinstance(a, np.ndarray)
        self.dim = a.shape[0]
        self.a = a
        self.b = b
        if not (
                all(isinstance(a_i, mpz) for a_i in np.append(a, b))
        ):
            self.convert_to_integer()
        self.minimial_representation()
        # norm of the composite vector (a,b)
        self.a_norm = np.linalg.norm([float(a_i) for a_i in self.a])
        self.matrix = self.as_matrix()

    def is_active(self, vertex):
        if np.dot(self.a, vertex) + self.b == 0:
            return True
        else:
            return False

    def normalize(self):
        """
        norms hyperplane such that b>=0.
        :return: True if hyperplane had to be flipped and False elsewise
        """
        if self.b < 0:
            self.a, self.b = -self.a, -self.b
            self.matrix = self.as_matrix()
            return True
        else:
            return False

    def norm(self):
        return np.linalg.norm(np.append(self.a.astype(float), float(self.b)))

    def convert_to_integer(self):
        if all([float(a_i).is_integer() for a_i in self.a]) and float(self.b).is_integer():
            self.a = np.array([mpz(a_i) for a_i in self.a])
            self.b = mpz(self.b)
            self.integer = True
            return
        k = 1e9
        norm = self.norm()
        self.a = self.a / norm
        self.a = np.array([mpz(k * a) for a in self.a])
        self.b = mpz(k * (self.b / norm))
        self.integer = True

    def as_matrix(self):
        import mpmath
        return mpmath.matrix(self.a)

    def pertubate(self):
        sigma = np.array([1e-6] * self.dim)
        covariance = np.diag(sigma ** 2)
        a = np.random.multivariate_normal(self.a / self.a_norm, covariance)
        return Hyperplane(a, self.b / self.a_norm)

    def minimial_representation(self):
        gcd = self.gcd()
        self.divide_by(gcd)

    def gcd(self):
        hyperplane_gcd = self.a[0]
        for a_i in np.append(self.a[1:], self.b):
            hyperplane_gcd = gcd(hyperplane_gcd, a_i)
        return hyperplane_gcd

    def divide_by(self, div):
        self.a = np.array([divexact(a_i, div) for a_i in self.a])
        self.b = divexact(self.b, div)
        self.matrix = self.as_matrix()
        self.a_norm = np.linalg.norm([float(a_i) for a_i in self.a])

    def __str__(self):
        """
        Output method.
        """

        cx = ['{}*x{}'.format(c, i) for i, c in enumerate(self.a, 1)]

        outputStr = " + ".join(cx)
        outputStr += " + {} == 0".format(self.b)

        return outputStr

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __deepcopy__(self, memo):
        # We make a shallow copy here to be more memory efficient
        # and still be able to run deepcopy on cell decomposition
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result


class Vertex(object):
    def __init__(self, ppl_point):
        assert ppl_point.is_point()
        self.dim = ppl_point.space_dimension()
        self.ppl = ppl_point
        self.vector = ppl_point.coefficients()
        self.denominator = ppl_point.divisor()

        self.coordinates = np.array(self.vector) / float(self.denominator)

    def position(self, hyperplane):
        x = [Variable(i) for i in range(self.dim)]

        p = C_Polyhedron(Generator_System(self.ppl))
        s_rel = p.relation_with(
            sum(hyperplane.a[i] * x[i] for i in range(self.dim)) + hyperplane.b < 0)
        if s_rel.implies(Poly_Con_Relation.is_included()):
            return -1
        else:
            b_rel = p.relation_with(
                sum(hyperplane.a[i] * x[i] for i in range(self.dim)) + hyperplane.b > 0)
            if b_rel.implies(Poly_Con_Relation.is_included()):
                return 1
            else:
                return 0

    def __hash__(self):
        return hash('<{}/{}>'.format(self.vector, self.denominator))

    def __eq__(self, other):
        return self.__hash__() == other.__hash__()

    def __str__(self):
        return '[' + ','.join(['{:.3f}'.format(c) for c in self.coordinates]) + ']'

    @staticmethod
    def vertex_from_coordinates(coordinates):
        x = [Variable(i) for i in range(len(coordinates))]
        denom = 1e+10
        p = point(sum(x[i] * int(coordinates[i] * denom) for i in range(len(coordinates))), denom)
        return Vertex(p)


class Cone(object):
    def __init__(self, halfspaces, save_halfspaces=False, incidences=None):
        if save_halfspaces:
            self.halfspaces = halfspaces
        else:
            self.halfspaces = None
        self.dim = halfspaces[0][0].dim
        self.vertex = None
        self.rays = self.rays_from_constraints(halfspaces)
        self.determinant = self.compute_determinant()
        self.incidences = incidences

    def compute_determinant(self):
        return np.linalg.det([ray.astype(float) for ray in self.rays])

    def rays_from_constraints(self, halfspaces):
        # Create PPL variables.
        x = [Variable(i) for i in range(self.dim)]

        # Init constraint system.
        constraints = Constraint_System()

        # Add polytope facets to constraint systems.
        # Format of ConvexHull.equations: [a1 a2 .. b] => a1*x1 + a2*x2 + .. + b <= 0.
        for hyp, orient in halfspaces:

            # PPL works on integral values, so all
            # coordinates are scaled and truncated.

            if orient == -1:
                constraints.insert(sum(hyp.a[i] * x[i] for i in range(self.dim)) + hyp.b <= 0)
            elif orient == 1:
                constraints.insert(sum(hyp.a[i] * x[i] for i in range(self.dim)) + hyp.b >= 0)
            elif orient == 0:
                constraints.insert(sum(hyp.a[i] * x[i] for i in range(self.dim)) + hyp.b == 0)

        # Build PPL polyhedra.
        poly = C_Polyhedron(constraints)
        # Get vertices of polytope resulting from intersection. only rays
        rays = []
        for gen in poly.minimized_generators():
            if gen.is_ray():
                ray = np.array(gen.coefficients()).astype(float)
                ray = ray / np.linalg.norm(ray)
                rays.append(ray)
        return np.array(rays)

    def is_neutralizing_cone(self, cone):
        """
        Method checks if cones are identical except for one flipped halfspace
        :param cone: another Cone object
        :return: True/False
        """
        shared_incidences = [[], []]
        for inc, orient in self.incidences.items():
            orient2 = cone.incidences.get(inc)
            if orient2:
                shared_incidences[0].append(orient)
                shared_incidences[1].append(orient2)
        if len(shared_incidences[0]) != self.dim:
            return False
        # if the vertices differ at exactly 1 entry, the dot product is exactly d-2
        # that is because all entries are in {-1,1}
        if np.dot(shared_incidences[0], shared_incidences[1]) == self.dim - 2:
            return True
        return False

    def __str__(self):
        rays = 'rays: ' + str(self.rays)

        if self.halfspaces:
            halfspaces = 'halfspace: ' + str(self.halfspaces)
            return rays + halfspaces
        else:
            return rays

    def __deepcopy__(self, memo):
        # We make a shallow copy here to be more memory efficient
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result


class Polytope(object):

    def __init__(self, halfspaces=None, vertices=None):

        if halfspaces is not None:
            self.halfspaces = halfspaces
            self.dim = halfspaces[0][0].dim
            self.poly = self.poly_from_constraints(halfspaces)
            self.vertices = self.vertices_from_poly()

        elif vertices is not None:
            v = vertices.pop()
            self.dim = v.dim
            vertices.add(v)
            self.poly = self.poly_from_vertices(vertices)
            self.vertices = self.vertices_from_poly()
            self.hyperplanes = self.hyperplanes_from_poly()
            self.halfspaces = zip(self.hyperplanes, [1] * len(self.hyperplanes))
        else:
            self.hyperplanes = []
            self.halfspaces = []
            self.vertices = []
            self.poly = C_Polyhedron(0)

    def as_convex_hull(self):
        vertices = np.array([vertex.coordinates for vertex in self.vertices])
        convex_hull = ConvexHull(vertices, qhull_options='Pp')
        return convex_hull

    def get_vertex_coordinates(self):
        return np.array([v.coordinates for v in self.vertices])

    def vertices_from_poly(self):
        # Get vertices of polytope resulting from intersection. only points not rays
        return np.array([Vertex(pt) for pt in self.poly.minimized_generators() if pt.is_point()])

    def poly_from_constraints(self, halfspaces):
        # Create PPL variables.
        x = [Variable(i) for i in range(self.dim)]

        # Init constraint system.
        constraints = Constraint_System()

        # Add polytope facets to constraint systems.
        # Format of ConvexHull.equations: [a1 a2 .. b] => a1*x1 + a2*x2 + .. + b <= 0.
        for hyp, orient in halfspaces:

            if orient == -1:
                constraints.insert(sum(hyp.a[i] * x[i] for i in range(self.dim)) + hyp.b <= 0)
            elif orient == 1:
                constraints.insert(sum(hyp.a[i] * x[i] for i in range(self.dim)) + hyp.b >= 0)
            elif orient == 0:
                constraints.insert(sum(hyp.a[i] * x[i] for i in range(self.dim)) + hyp.b == 0)

        # Build PPL polyhedra.
        return C_Polyhedron(constraints)

    def poly_from_vertices(self, vertices):
        # Create PPL variables.
        x = [Variable(i) for i in range(self.dim)]

        # Init constraint system.
        generators = Generator_System()

        # ppl needs coordinates in form of point(sum(a_i *x_i), denom) with a_i integers
        for vertex in vertices:
            generators.insert(vertex.ppl)
        # Build PPL polyhedra.
        return C_Polyhedron(generators)

    def hyperplanes_from_poly(self):

        logging.debug('<nr generators: {}>, <nr constraints: {}>'
                      .format(len(self.poly.generators()), len(self.poly.constraints())))

        hyperplanes = []

        # constraints in ppl are saved as of the form ax + b >= 0

        for constraint in self.poly.minimized_constraints():
            a = np.array(constraint.coefficients())
            b = constraint.inhomogeneous_term()
            hyperplane = Hyperplane(a, b)
            hyperplanes.append(hyperplane)
        return hyperplanes

    def add_constraint(self, halfspace):
        # Create PPL variables.
        x = [Variable(i) for i in range(self.dim)]
        if halfspace[1] == -1:
            self.poly.add_constraint(
                sum(halfspace[0].a[i] * x[i] for i in range(self.dim)) + halfspace[0].b <= 0)
        elif halfspace[1] == 1:
            self.poly.add_constraint(
                sum(halfspace[0].a[i] * x[i] for i in range(self.dim)) + halfspace[0].b >= 0)
        elif halfspace[1] == 0:
            self.poly.add_constraint(
                sum(halfspace[0].a[i] * x[i] for i in range(self.dim)) + halfspace[0].b == 0)
        self.hyperplanes = self.hyperplanes_from_poly()
        self.halfspaces.append(halfspace)
        self.vertices = self.vertices_from_poly()

    def point_in_poly(self, point):
        return self.poly.relation_with(point.ppl).implies(Poly_Gen_Relation.subsumes())


def vector_distance(vector_1, vector_2, absolute_value=True):
    """
    Function calculates the distance between two vectors as 1 - <a1, a2>
    :param vector_1: numpy array
    :param vector_2: numpy array
    :param absolute_value: Boolean if absolute value should be returned.
    :return: vector distance of vector_1 and vector_2
    """
    a1_normed = vector_1 / np.linalg.norm(vector_1)
    a2_normed = vector_2 / np.linalg.norm(vector_2)
    if absolute_value:
        return 1 - abs(np.dot(a1_normed, a2_normed))
    else:
        return 1 - np.dot(a1_normed, a2_normed)


def norm_halfspace(halfspace):
    """
    Method norms halfspace in a way that offset is non-negative
    :param halfspace: a tuple of (hyperplane, orientation) whereby -1 represents <= and +1 >= for
    the orientation.
    :return: Normed halfspace with positive offset
    """
    hyperplane, orientation = halfspace[0], halfspace[1]
    if hyperplane.b < 0:
        hyperplane.a = -hyperplane.a
        hyperplane.b = -hyperplane.b
        hyperplane.matrix = hyperplane.as_matrix()
        return hyperplane, -orientation
    else:
        return halfspace