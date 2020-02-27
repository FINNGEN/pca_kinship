"""
A module for testing/demonstrating the use of the py:mod:alluvial module.
Have a look at the source.
"""
import os

from alluvial import plot_alluvial


def test():
    import matplotlib as mpl
    mpl.use('Agg')
    from matplotlib import pyplot as plt

    ribbon_size_matrix = [
        [5, 2, 1],
        [3, 4, 1],
        [1, 0, 3]
    ]

    ribbon_label_matrix = [
        ["a", "b", "c"],
        ["short", "medium", "long"],
        ["3-1", "3-2", "3-3"]
    ]

    stable_ribbon_sizes_1 = [
        [2, 1, 1],
        [3, 4, 1],
        [1, 0, 2]
    ]

    stable_ribbon_sizes_2 = [
        [5, 2, 1],
        [3, 2, 1],
        [1, 0, 3]
    ]

    module_colors_1 = ["red", "green", "blue"]
    module_colors_2 = ["green", "red", "blue"]
    ribbon_bglim = 1.5

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    # the current full power:
    plot_alluvial(
        ax1, ribbon_size_matrix, ribbon_label_matrix, module_colors_1,
        module_colors_2, ribbon_bglim, stable_ribbon_sizes_1,
        stable_ribbon_sizes_2, ribbon_label_size=10, ribbon_label_HA="center"
    )

    #  minimal setup
    plot_alluvial(ax2, ribbon_size_matrix)

    fname = os.path.join(os.path.dirname(__file__),
                         "gallery/alluvial_example.pdf")

    # save fig to gallery
    plt.savefig(fname, format='pdf')
