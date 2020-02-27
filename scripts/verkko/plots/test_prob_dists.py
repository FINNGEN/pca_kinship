import os
import numpy as np
import matplotlib.pyplot as plt

import prob_dists


def test():
    n_points = 10000
    xvals = np.exp(np.random.randn(n_points))
    fig = plt.figure()
    ax = fig.add_subplot(311)
    fig, ax = prob_dists.plot_pdf(xvals, xscale='log', ax=ax)
    ax = fig.add_subplot(312)
    fig, ax = prob_dists.plot_pdf(xvals, xscale='lin', yscale='log', ax=ax)
    ax = fig.add_subplot(313)
    prob_dists.plot_pdf(xvals, xscale='lin', ax=ax)
    fname = os.path.join(os.path.dirname(__file__), "gallery/PDF_example.pdf")
    fig.savefig(fname)

    fig = plt.figure()
    ax = fig.add_subplot(311)
    fig, ax = prob_dists.plot_ccdf(xvals, xscale='lin', yscale='log', ax=ax)
    ax = fig.add_subplot(312)
    fig, ax = prob_dists.plot_ccdf(xvals, xscale='log', yscale='log', ax=ax)
    ax = fig.add_subplot(313)
    fig, ax = prob_dists.plot_ccdf(xvals, xscale='lin', yscale='lin', ax=ax)
    fname = os.path.join(os.path.dirname(__file__), "gallery/cCDF_example.pdf")
    fig.savefig(fname)
