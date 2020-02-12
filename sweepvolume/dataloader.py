# *****************************************************************************
#       Copyright (C) 2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
import json
import logging
from sweepvolume.geometry import Hyperplane, Polytope, norm_halfspace
from sweepvolume.cell_decomposition import Cell_Decomposition
from ordered_set import OrderedSet
import numpy as np


def load_json(filepath):
    # Open and read specified input file.
    #  Input format has to be {disj: {polyID1: [...], ..., polyIDk: [...]}}

    with open(filepath, 'r') as f:
        disj = json.load(f)
        logging.debug('Loaded disjunction from file: {}'.format(filepath))
    if len(disj) > 1:
        logging.warning("JSON input file contains multiple disjunctions -"
                        " only first will be considered!")
    return disj


def polytopes_from_json(filepath):
    """
    Method to read polytope points from JSON file.
     Input format is {disj: {polyID1: [...], ..., polyIDk: [...]}}
     The elements polyID: [[a1, ..., ad, b],...,[...]]
     correspond to the halfspace a1x1 + ... + adxd + b <= 0
    """
    polys = []
    disj = load_json(filepath)
    for disjID, _ply in disj.items():
        logging.info("Reading polytopes of disjunction {}".format(disjID))

        for i, (plyID, ply) in enumerate(_ply.items()):
            halfspaces = [(Hyperplane(np.array(hyp_vec[:-1]), hyp_vec[-1]), -1)
                          for hyp_vec in ply]
            polys.append(Polytope(halfspaces=halfspaces))

    return polys


# Method computes cell decompositions for union of polytopes
def get_cell_decomposition(filepath):
    polytopes = polytopes_from_json(filepath)
    hyperplanes, polytope_vectors = hyperplanes_and_polytope_vectors(polytopes)

    logging.info("Computing cell decomposition for union of polytopes"
                 " from {} hyperplanes.".format(len(hyperplanes)))
    union_cd = Cell_Decomposition(hyperplanes,
                                  polytope_vectors)
    return union_cd


def hyperplanes_and_polytope_vectors(polytopes):
    """
    Method computes hyperplanes and position vectors for given polytopes
    :param polytopes:a list of list of tuples: (sweepvolume.geometry.Hyperplane, orientation)
    :return: hyperplanes: list of hyperplanes,
     polytopes as list of tuples (hyperplane_index, orientation)
    """
    assert all([isinstance(poly, Polytope) for poly in polytopes])
    hyperplanes = OrderedSet()
    polytope_vectors = []
    for poly in polytopes:
        ply = set()
        for halfspace in poly.halfspaces:
            halfspace = norm_halfspace(halfspace)
            idx = hyperplanes.add(halfspace[0])
            ply.add((idx, halfspace[1]))
        polytope_vectors.append(ply)
    return hyperplanes, polytope_vectors

