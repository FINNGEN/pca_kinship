import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.colors import ColorConverter
import matplotlib as mpl


def plot_alluvial(ax, ribbon_size_matrix, ribbon_label_matrix=None,
                  module_colors_1=None, module_colors_2=None,
                  ribbon_bglim=20, stable_ribbon_sizes_1=None,
                  stable_ribbon_sizes_2=None, ribbon_label_size=10,
                  ribbon_label_HA="center"):
    """
    Plot a two-sided alluvial diagram to the ax given as parameter.
    See :py:mod:`test_alluvial` for examples.

    Parameters
    ----------
    ax : a matplotlib.axes object
    ribbon_size_matrix : 2D matrix
        a 2D matrix, which describes how the groups change
        element ``ribbon_size_matrix[i,j]``, should correspond to the number
        of common nodes for the modules i and j.
    ribbon_label_matrix : 2D-iterable
        matrix, which contains the ribbon labels
    module_colors_1 : iterable, optional
        colors for the first (left) modules (defaulting to gray)
    module_colors_2 : iterable, optional
        colors for the second (right) modules (defaulting to gray)
    ribbon_bglim : float, optional
        if ribbon contains less than this number of nodes (or other 'mass')
        it is drawn on the background to avoid visual clutter
    stable_ribbon_sizes_1 : 2D np.array, optional
        a matrix describing the number of shaded nodes of each module for the
        first modules (defaulting to no shading)
    stable_ribbon_sizes_2 : 2D np.array, optional
        a matrix describing the number of shaded nodes of each module for the
        second modules (defaulting to no shading)
    ribbon_label_size : number, optional
        the font size for the ribbon labels
    ribbon_label_HA : {"center", "left", "right"}, optional
        the horizontal alignment of the ribbon labels (defaulting to center)

    """
    module_sizes_1 = np.sum(ribbon_size_matrix, axis=1)
    module_sizes_2 = np.sum(ribbon_size_matrix, axis=0)
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_xticks([])
    ax.set_yticks([])

    for loc, spine in ax.spines.items():
        spine.set_color('none')  # don't draw spine

    if module_colors_1 is None:
        module_colors_1 = ["gray"] * len(module_sizes_1)
    if module_colors_2 is None:
        module_colors_2 = ["gray"] * len(module_sizes_2)

    assert len(module_sizes_2) == len(module_colors_2)
    assert len(module_sizes_1) == len(module_colors_1)

    # assumes identical total sizes of the modules
    #(i.e. same nodes on both sides)
    n_modules_1 = len(module_sizes_1)
    n_modules_2 = len(module_sizes_2)
    max_n_modules = np.max([n_modules_1, n_modules_2])
    module_size_sum = np.maximum(np.sum(module_sizes_1),
                                 np.sum(module_sizes_2))
    vertical_pad_btw_modules = 0.01  # percent of
    vertical_pad_up_and_below = 0.00
    horizontal_pad_lr = 0.05
    individual_node_size = (1 - 2 * vertical_pad_up_and_below -
                            vertical_pad_btw_modules *
                            (max_n_modules - 1)) / module_size_sum
    module_width = 0.2  # should be in range [0,0.5-horizontal_pad_lr
    mw = module_width
    # mwnsf: module width non shaded fraction
    #(to be able to differ from fully shaded to non-shaded)
    mwnsf = 0.1
    # plot first modules
    iteration_list = [
        [module_sizes_1, module_colors_1, stable_ribbon_sizes_1],
        [module_sizes_2, module_colors_2, stable_ribbon_sizes_2]
    ]

    # for storing the start y coordinates
    module_y_starts_1_and_2 = [[], []]
    # plot modules
    for i, iteration_data in enumerate(iteration_list):
        module_sizes, module_colors, stable_ribbon_sizes = iteration_data
        module_y_starts = module_y_starts_1_and_2[i]
        current_y = vertical_pad_up_and_below
        if i is 0:
            rect_x_start = horizontal_pad_lr
        else:
            rect_x_start = 1 - horizontal_pad_lr - mw
        for j in range(len(module_sizes)):
            module_size = module_sizes[j]
            color = module_colors[j]
            module_y_starts.append(current_y)
            module_height = individual_node_size * module_size

            rect = Rectangle((rect_x_start, current_y), module_width,
                             module_height, fc=color, ec="0.85")
            ax.add_patch(rect)
            if stable_ribbon_sizes_1 is not None:
                rect = Rectangle(
                    (rect_x_start + i * mw * (1 - mwnsf), current_y),
                    module_width * mwnsf, module_height,
                    fc=color, ec="0.85")
                ax.add_patch(rect)
            current_y += module_height + vertical_pad_btw_modules

    modules_current_fill_starts_1 = np.array(module_y_starts_1_and_2[0])
    modules_current_fill_starts_2 = np.array(module_y_starts_1_and_2[1])
    curvature_param = 0.6
    # plot ribbons in order of biggest modules first?
    zorder = 0

    for i in range(len(ribbon_size_matrix)):
        for j in range(len(ribbon_size_matrix[i])):
            ribbon_size = ribbon_size_matrix[i][j]
            if ribbon_size == 0:
                continue
            ystart1 = modules_current_fill_starts_1[i]
            yend1 = ystart1 + ribbon_size * individual_node_size
            modules_current_fill_starts_1[i] = yend1
            ystart2 = modules_current_fill_starts_2[j]
            yend2 = ystart2 + ribbon_size * individual_node_size
            modules_current_fill_starts_2[j] = yend2
            # the points needed for the bezier
            bezier_verts1 = [
                (horizontal_pad_lr + module_width, ystart1),  # P0
                (curvature_param, ystart1),  # P1
                (1 - curvature_param, ystart2),  # P2
                (1 - horizontal_pad_lr - module_width, ystart2),  # P3
            ]
            bezier_verts2 = [
                (horizontal_pad_lr + module_width, yend1),  # P0
                (curvature_param, yend1),  # P1
                (1 - curvature_param, yend2),  # P2
                (1 - horizontal_pad_lr - module_width, yend2),  # P3
            ]
            if ribbon_size < ribbon_bglim:
                use_zorder = -10000 - j
            else:
                use_zorder = zorder
            _plot_ribbon_using_bezier(ax, use_zorder, bezier_verts1,
                                      bezier_verts2, module_colors_1[i],
                                      module_colors_2[j])
            # plotting node stabilities
            if stable_ribbon_sizes_1 is not None and \
                    stable_ribbon_sizes_2 is not None:
                ribbon_vertical_size_1 = yend1 - ystart1
                ribbon_vertical_size_2 = yend2 - ystart2
                y1_low = ystart1 + (stable_ribbon_sizes_1[i][j] /
                                    float(ribbon_size) *
                                    ribbon_vertical_size_1)
                y1_high = yend1
                y2_low = ystart2 + (stable_ribbon_sizes_2[i][j] /
                                    float(ribbon_size) *
                                    ribbon_vertical_size_2)
                y2_high = yend2
                bezier_verts1_stable = [
                    (horizontal_pad_lr + module_width, y1_low),  # P0
                    (curvature_param, y1_low),  # P1
                    (1 - curvature_param, y2_low),  # P2
                    (1 - horizontal_pad_lr - module_width, y2_low),  # P3
                ]
                linewidth = 0.25
                unstable_mask_color = (1, 1, 1, 0.3)
                _plot_ribbon_using_bezier(ax, use_zorder, bezier_verts1_stable,
                                          bezier_verts2, unstable_mask_color,
                                          unstable_mask_color, lw=linewidth)
                rect = Rectangle((horizontal_pad_lr + mwnsf * mw, y1_low),
                                 mw * (1 - mwnsf), y1_high - y1_low,
                                 fc=unstable_mask_color,
                                 ec=unstable_mask_color, lw=linewidth)
                ax.add_patch(rect)
                rect = Rectangle(
                    (1 - horizontal_pad_lr - module_width, y2_low),
                    mw * (1 - mwnsf), y2_high - y2_low,
                    fc=unstable_mask_color,
                    ec=unstable_mask_color, lw=linewidth)
                ax.add_patch(rect)
            zorder = zorder - 1
            if ribbon_label_matrix is not None:
                s = ribbon_label_matrix[i][j]
                if s != "":
                    ax.text(
                        horizontal_pad_lr +
                        module_width, (ystart1 + yend1) / 2.,
                        s, color="0.1", va='center', ha=ribbon_label_HA,
                        fontsize=ribbon_label_size, family="arial",
                        fontweight='bold')


def _plot_ribbon_using_bezier(ax, zorder, points1, points2, color1="gray",
                              color2="gray", lw=1):
    """ Draw ribbon for alluvial diagram (see plot_alluvial)

    Parameters
    ----------
    ax : a matplotlib.axes object
    zorder : float
        the zorder for the ribbon
    points1 : iterable of float tuples
        the points, which determine the first line of the Bezier ribbon
    points2 : iterable of float tuples
        the points, which determine the second line of the Bezier ribbon
    color1 : a matplotlib compliant color definition
        color for the left side of the ribbon
    color1 : a matplotlib compliant color definition
        color for the right side of the ribbon
    lw : float
        linewidth for the bezier borders
    """
    cc = ColorConverter()
    color1 = np.array(cc.to_rgba(color1))
    color2 = np.array(cc.to_rgba(color2))
    tRange = np.linspace(0, 1, 100)
    xpointsList = []
    ypointsList = []
    for points in [points1, points2]:
        points = np.array(points)
        p1 = points[0]
        p2 = points[1]
        p3 = points[2]
        p4 = points[3]
        allPoints = (p1[:, np.newaxis] * (1 - tRange) ** 3 + p2[:, np.newaxis]
                     * (3 * (1 - tRange) ** 2 * tRange) + p3[:, np.newaxis] *
                     (3 * (1 - tRange) * tRange ** 2) + p4[:, np.newaxis] *
                     tRange ** 3)
        xpoints = allPoints[0]
        xpointsList.append(xpoints)
        ypoints = allPoints[1]
        ypointsList.append(ypoints)
        ax.plot(xpoints, ypoints, "0.85", lw=lw, zorder=zorder + 0.5)
    xpoints = xpointsList[0]
    if (mpl.colors.colorConverter.to_rgba_array(color1) ==
            mpl.colors.colorConverter.to_rgba_array(color2)).all():
        ax.fill_between(xpoints, ypointsList[0], ypointsList[1], lw=lw,
                        facecolor=color1, edgecolor=color1, zorder=zorder)
    else:
        for i in range(len(tRange) - 1):
            #mean = (tRange[i]+tRange[i+1])*0.5
            xnow = np.mean(xpoints[i:i + 2])
            norm_mean = (xnow - xpoints[0]) / (xpoints[-1] - xpoints[0])
            color = color1 * (1 - norm_mean) + color2 * norm_mean
            ax.fill_between(xpoints[i:i + 2], ypointsList[0][i:i + 2],
                            ypointsList[1][i:i + 2], lw=lw, facecolor=color,
                            edgecolor=color, zorder=zorder)
