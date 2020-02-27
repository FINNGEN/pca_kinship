"""Calculate and plot statistics about a node cover.
Usage:
    python nodeCover_stats.py COMMFILE RESNAME [NETFILE]
"""

import sys
from netpython import netio, communities, transforms
import binner
import pylab
import numpy as np
import fig_utils

# Some general parameters for the plots.
fig_width_cm = 14 # The width of figures in cm.
font_sizes = (7,8,7,7,7)
print_titles = True
save_formats = ['png']

def print_statistics(nc, ncName, net=None):
    """Print statistics about the node cover and the network."""
    with open(ncName + "_stats.txt", 'w') as f:
        f.write("Node cover file:\n")
        f.write(sys.argv[1] + "\n\n")
        if net:
            f.write("Network file:\n")
            f.write(sys.argv[3] + "\n\n")

        f.write("Number of communities:   %d\n" % len(nc))
        f.write("Largest community size:  %d nodes\n" % nc.getGiantSize())
        f.write("Smallest community size: %d nodes\n" % nc.getCommunitySizes()[-1])

def plot_avgCommWeight(nc, ncName, net):
    """Plot average weight of communities."""

    def w_gen(comm, net):
        w_avg_total = np.mean(list(net.weights))
        for c in comm:
            subnet = transforms.getSubnet(net,c)
            if len(subnet):
                w = np.sum(list(subnet.weights))/(len(subnet)*w_avg_total)
                yield len(subnet), w

    bins = binner.Bins(int, 1, nc.getGiantSize(), 'linlog', 1.5)
    perc = (0.1, 0.25, 0.5, 0.75, 0.9)
    binned_data = bins.bin_percentile(w_gen(nc, net), perc)

    # Set plotting properties
    l, b, r, t = 0.1, 0.1, 0.98, (0.93 if print_titles else 0.98)
    axes_rect = [l,b,r-l,t-b] #[left, bottom, width, height]
    pylab.rcParams.update(fig_utils.get_rcParams(float(fig_width_cm)),
                          fig_ratio=0.8,
                          font_sizes=font_sizes)
    
    fig = pylab.figure()
    ax = fig.add_axes(axes_rect)

    plot_styles = ('b:', 'b--', 'r-', 'b--', 'b:')
    labels = (r"%d$^{\mathrm{th}}$ percentile" % int(100*perc[0]),
              r"%d$^{\mathrm{th}}$ percentile" % int(100*perc[1]),
              r"Median",
              r"%d$^{\mathrm{th}}$ percentile" % int(100*perc[3]),
              r"%d$^{\mathrm{th}}$ percentile" % int(100*perc[4]))
    for p, data, sty, lbl in zip(perc, binned_data, plot_styles, labels):
        ax.loglog(bins.centers, data, sty, label=lbl)

    ax.set_xlabel(r"Community size")
    ax.set_ylabel(r"$\langle w_{\mathrm{comm}} \rangle /"
                  r"\langle w_{\mathrm{total}} \rangle$")
    if print_titles:
        ax.set_title(r"Ratio of avg community weigth to avg total weight (binned)")

    ax.legend(loc='best')

    # Adjust axis
    V = ax.axis()
    ax.axis((1, bins.bin_limits[-1], V[2], V[3]))

    # Save figure.
    fig_utils.savefig(fig, ncName + "_weightDist", save_formats)



def plot_commSizeDist(nc, ncName):
    """Plot community size distribution."""
    bins = binner.Bins(int, 1, nc.getGiantSize(), 'linlog', 1.5)
    binned_data = bins.bin_count_divide(nc.getCommunitySizes())

    # Set plotting properties
    l, b, r, t = 0.1, 0.1, 0.98, (0.93 if print_titles else 0.98)
    axes_rect = [l,b,r-l,t-b] #[left, bottom, width, height]
    pylab.rcParams.update(fig_utils.get_rcParams(float(fig_width_cm)),
                          fig_ratio=0.8,
                          font_sizes=font_sizes)
    
    fig = pylab.figure()
    ax = fig.add_axes(axes_rect)
    ax.loglog(bins.centers, binned_data, 'b.-')

    ax.set_xlabel(r"Community size")
    ax.set_ylabel(r"Number of communities")
    if print_titles:
        ax.set_title(r"Community size distribution (binned)")

    # Adjust axis
    V = ax.axis()
    ax.axis((1, bins.bin_limits[-1], V[2], V[3]))

    # Save figure.
    fig_utils.savefig(fig, ncName + "_sizeDist", save_formats)

def plot_all(nc, ncName, net=None):
    """Create all plots of node cover.

    If `net` is given, the plots that require the underlying network
    will also be done.

    Parameters
    ----------
    nc : netpython.communities.NodeCover
        The node cover.
    ncName : str
        The name of the node cover. This will be used as the beginning
        of all file names.
    net : pynet.SymmNet
        The network under `nc`.
    """
    print_statistics(nc, ncName, net)
    plot_commSizeDist(nc, ncName)
    if net is not None:
        plot_avgCommWeight(nc, ncName, net)

if __name__ == '__main__':
    # Read in the community structure.
    try:
        commFileName = sys.argv[1]
        ncName = sys.argv[2]
        netFileName = ""
        if len(sys.argv) > 3:
            netFileName = sys.argv[3]
    except:
        sys.stderr.write(__doc__)

    # Read in the community structure.
    with open(commFileName, 'r') as f:
        nc = communities.NodeCover(inputFile=f)

    # Read in the network (if given)
    net = None
    if netFileName:
        net = netio.loadNet(netFileName)

    # Plot everything.
    plot_all(nc, ncName, net)

