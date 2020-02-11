# *****************************************************************************
#       Copyright (C) 2017-2018 Lovis Anderson  <lovisanderson@gmail.com>
#                     2017-2018 Benjamin Hiller <hiller@zib.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 3 of
#  the License, or (at youroption) any later version.
#                  http://www.gnu.org/licenses/
# *****************************************************************************
from scipy.spatial import ConvexHull

from sweepvolume.geometry import Polytope

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as mpl3D
import numpy as np
import pylab as pl

# well distinguishable colors
FACECOLORS = np.array([[0.90196078, 0.09803922, 0.29411765, 0.5],
                       [0.23529412, 0.70588235, 0.29411765, 0.5],
                       [1., 0.88235294, 0.09803922, 0.5],
                       [0., 0.50980392, 0.78431373, 0.5],
                       [0.96078431, 0.50980392, 0.18823529, 0.5],
                       [0.56862745, 0.11764706, 0.70588235, 0.5],
                       [0.2745098, 0.94117647, 0.94117647, 0.5],
                       [0.94117647, 0.19607843, 0.90196078, 0.5],
                       [0.82352941, 0.96078431, 0.23529412, 0.5],
                       [0.98039216, 0.74509804, 0.74509804, 0.5],
                       [0., 0.50196078, 0.50196078, 0.5],
                       [0.90196078, 0.74509804, 1., 0.5],
                       [0.66666667, 0.43137255, 0.15686275, 0.5],
                       [1., 0.98039216, 0.78431373, 0.5],
                       [0.50196078, 0., 0., 0.5],
                       [0.66666667, 1., 0.76470588, 0.5],
                       [0.50196078, 0.50196078, 0., 0.5],
                       [1., 0.84313725, 0.70588235, 0.5],
                       [0., 0., 0.50196078, 0.5],
                       [0.50196078, 0.50196078, 0.50196078, 0.5],
                       # From This line on colors are just randomly generated
                       [0.95329486, 0.66512501, 0.01362486, 0.5],
                       [0.62528863, 0.67630143, 0.97575655, 0.5],
                       [0.88163737, 0.5116229, 0.05593318, 0.5],
                       [0.45292847, 0.02006605, 0.44344312, 0.5],
                       [0.98342825, 0.36085405, 0.48277939, 0.5],
                       [0.69136181, 0.88392874, 0.92183639, 0.5],
                       [0.21767242, 0.56740529, 0.86849512, 0.5],
                       [0.51096492, 0.92031795, 0.92476999, 0.5],
                       [0.08343842, 0.27880765, 0.0093934, 0.5],
                       [0.84564538, 0.64971208, 0.84468567, 0.5]
                       ])


def plot_convex_hull_cell_decompositions(cds, alpha=1.):
    assert all(len(cd.polytope_vectors) == 1 for cd in cds)
    polytopes_vertices = []
    convex_hulls = []
    facecolors = FACECOLORS
    for i, cd in enumerate(cds):
        polytope = cd.polytope_vectors[0]
        halfspaces = [[cd.hyperplanes[j], orientation] for (j, orientation) in polytope]
        p = Polytope(halfspaces)
        # try except because scipy.ConvexHull cant handle non-full dimensional polytopes
        try:
            h = ConvexHull([vertex.coordinates for vertex in p.vertices])
            convex_hulls.append(h)
            polytopes_vertices.append(h.points[h.vertices])
        except:
            # workaround s.t. facecolor for lower dimensional polytope is not used
            #  -> consistent facecolors in ACD
            facecolors = np.delete(facecolors, i, 0)
            continue
    facecolors[:, -1] = alpha
    if cds[0].dim == 2:
        polycollection = mpl.collections.PolyCollection(np.array(polytopes_vertices),
                                                        facecolors=facecolors)
        fig, ax = plt.subplots(figsize=(7, 7))

        axes = plt.gca()

        ax.add_collection(polycollection)
        limits = get_limits(polytopes_vertices)

        axes.set_xlim(limits[0])
        axes.set_ylim(limits[1])

    elif cd.dim == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for i, poly in enumerate(polytopes_vertices):
            plot3d_polytope(poly, ax, FACECOLORS[i])

        ax.w_xaxis.set_ticklabels([])
        ax.w_yaxis.set_ticklabels([])
        ax.w_zaxis.set_ticklabels([])

    plt.show()


def plot_polytopes_and_cut(cd,
                           cuts=None,
                           alpha=0.5,
                           box=None,
                           plot_convex_hull=True, ax=None):
    polytopes_vertices = []
    convex_hulls = []
    facecolors = FACECOLORS
    for i, polytope in enumerate(cd.polytope_vectors):
        halfspaces = [[cd.hyperplanes[j], orientation] for (j, orientation) in polytope]
        p = Polytope(halfspaces)
        # try except because scipy.ConvexHull cant handle non-full dimensional polytopes
        try:
            h = ConvexHull([vertex.coordinates for vertex in p.vertices])
            convex_hulls.append(h)
            polytopes_vertices.append(h.points[h.vertices])
        except:
            # workaround s.t. facecolor for lower dimensional polytope is not used
            #  -> consistent facecolors in ACD
            facecolors = np.delete(facecolors, i, 0)
            continue

    facecolors[:, -1] = alpha

    if plot_convex_hull:
        convex_hull = ConvexHull(np.concatenate(polytopes_vertices))
        convex_hull_vertices = convex_hull.points[convex_hull.vertices]
        polytopes_vertices.append(convex_hull_vertices)
        FACECOLORS[len(polytopes_vertices) - 1] = (0, 0, 0, 0.2)

    axis_given = not (ax is None)

    if cd.dim == 2:
        polycollection = mpl.collections.PolyCollection(np.array(polytopes_vertices),
                                                        facecolors=facecolors)
        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 7))
        axes = plt.gca()

        if not box:
            box = get_limits(polytopes_vertices)

        axes.set_xlim(box[0])
        axes.set_ylim(box[1])
        ax.add_collection(polycollection)
        if cuts:
            if not isinstance(cuts, list):
                cuts = [cuts]
            for cut in cuts:
                if cut.a[1] != 0:
                    x = np.arange(box[0][0], box[0][1], 0.01)
                    y = (-cut.b - cut.a[0] * x) / cut.a[1]
                else:
                    y = np.arange(box[1][0], box[1][1], 0.01)
                    x = (-cut.b - cut.a[1] * y) / cut.a[0]
                plt.plot(x, y, color='k', linewidth=4.)

    elif cd.dim == 3:

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for i, poly in enumerate(polytopes_vertices):
            plot3d_polytope(poly, ax, FACECOLORS[i])

        ax.w_xaxis.set_ticklabels([])
        ax.w_yaxis.set_ticklabels([])
        ax.w_zaxis.set_ticklabels([])

    if axis_given:
        return ax
    else:
        plt.show()


def plot3d_polytope(points, ax, plot_color):
    """
    Plotting method for 3d polytope
    :param points: points of polytope.
    :param ax: matplotlib axis.
    :param plot_color: color, either color string or rgb tuple (between 0 and 1).
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        show = True
    else:
        show = False

    # Plot points.
    ax.scatter(points[:, 0],
               points[:, 1],
               points[:, 2],
               marker='')
    c = ConvexHull(points)
    # Plot convex hull simplices.
    hull3d = [points[spx] for spx in c.simplices]
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    hullColl = Poly3DCollection(hull3d, linewidths=0.2, edgecolors=plot_color, alpha=plot_color[-1])
    hullColl.set_facecolor(plot_color)
    ax.add_collection3d(hullColl)
    if show:
        plt.show()


def get_limits(polytopes_vertices):
    dimension_flattend = [[coord[i] for polytope in polytopes_vertices for coord in polytope]
                          for i in range(len(polytopes_vertices[0][0]))]
    limits = []
    for i, flat in enumerate(dimension_flattend):
        frame = (max(flat) - min(flat)) / 10.
        limits.append([min(flat) - frame,
                       max(flat) + frame])

    return limits


def plot_sweep(sweep, points=10000):
    start_sweep = sweep.sorted_events[0][1]
    end_sweep = sweep.sorted_events[-1][1]
    sweep_length = end_sweep - start_sweep
    step_length = sweep_length / points
    start_plot = start_sweep - 50*step_length
    end_plot = end_sweep + 50*step_length
    x = np.arange(start_plot, end_plot, step_length)
    y = sweep.calculate_volumes(x)
    fig, ax = plt.subplots()
    ax.plot(x, y)

    ax.set(xlabel='Sweep-Time', ylabel='Volume',
           title=f'Sweep graph for direction {sweep.sweep_plane}')
    plt.show()
