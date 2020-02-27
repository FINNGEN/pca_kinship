"""
This script takes in a file created by the com_struct_topol_analysis.py script with the format:
N_NODES N_COMMUNITIES EDGE_DENSITY EDGE_WEIGHT INTENSITY COHERENCE SHORTEST_PATH_LENGTH DIAMETER

and plots the wanted values. Plots WILL NOT BE SAVED unless you give an output directory as argument.

Under the heading: YOUR VARIABLES, you can flag which plots you want to be done.
Knock your socks off, by modifying the plotting scheme according to your wish, and save the file for your favourite net of choice. Now several values are plotted with both scatter and binned plots.

Note: input file and output dir. need to be with absolute path or be in working dir. (No ../ stuff)

Examples of use:
1. Just show plots on screen:
python plot_com_struct_topol_analysis.py n_sms_and_calls_as_w_average_community_analysis.txt

2. Show and save plots as .eps-files with the following name format:
{plotted analysis}{input file name}.eps
python plot_com_struct_topol_analysis.py n_sms_and_calls_as_w_average_community_analysis.txt /home/eiesland/networks/Jon/communities/python_code/plots/black_metal_is_the_best_music/

Jon Eiesland
May 2009
"""

from pylab import *
from binner import * # lce-made binning module for plotting also found in the Python_aux in the cvs
from numpy import *
from sys import argv

def main():
    #argv and read inputfile
    if (len(argv) < 3):
       print '\n PLOTS WILL NOT BE SAVED. If you want them saved, give absolute path (and ending with \"/\") of output directory as argument.'
       saveplots = 0
    analysis_file = argv[1]
    #if dir provided, flag to save plots
    if (len(argv) == 3):
       out_dir = argv[2]
       saveplots = 1
    f=open(analysis_file)
    #read header and discard
    f.readline()
    #read rest to memory
    s=f.readlines()
    f.close()

    save_plot_name = analysis_file.split('/')[-1]
    save_plot_name = analysis_file.split('.')[0]

    #make plot vectors and insert from file
    n_nodes = []
    n_com = []
    edge_dens = []
    edge_weight = []
    intensity = []
    coherence = []
    shortest_pl = []
    diameter = []
    for line in s:
        foo = line.split()
        n_nodes.append(int(foo[0]))
        n_com.append(int(foo[1]))
        edge_dens.append(float(foo[2]))
        edge_weight.append(float(foo[3]))
        intensity.append(float(foo[4]))
        coherence.append(float(foo[5]))
        shortest_pl.append(float(foo[6]))
        diameter.append(float(foo[7]))

    #Make all the figures
    #YOUR VARIABLES
    do_show = 1 # show plots on screen or not
    do_all = 1 # plots all if flagged, regardles of following flags
    do_comdist = 0
    do_edgedens = 0
    do_edgeweight = 0
    do_intensity = 1
    do_coherence = 1
    do_sh_pl = 1
    do_diam = 1

    clique_number = 4 # percolation number used in the com analysis. for x-axis label
    xlabel_text = '# of nodes in ' + str(clique_number) + '-clique community'
    weight_unit = '(# of calls and SMS)' #for weight dist y-axis label


    if (do_comdist or do_all):
        figure() #com dist
        xlabel(xlabel_text)
        ylabel('# of communities')
        loglog(n_nodes, n_com, '*')
        if (saveplots):
            savefig(out_dir + 'com_dist_' + save_plot_name + '.eps')

    if (do_edgedens or do_all):
        figure() #edge density
        xlabel(xlabel_text)
        ylabel('Average edge density')
        loglog(n_nodes, edge_dens, '*')
        if (saveplots):
            savefig(out_dir + 'edge_dens_' + save_plot_name + '.eps')

    if (do_edgeweight or do_all):
        bins = Bins(int, min(n_nodes), max(n_nodes), 'log', 1.6)
        seq = zip(n_nodes, edge_weight)
        bin_av = bins.bin_average(seq)
        figure() #edge weight
        subplot(121)
        loglog(n_nodes, edge_weight, '*')
        xlabel(xlabel_text)
        ylabel('Average edge weight ' + weight_unit)
        ymin, ymax = ylim()
        subplot(122)
        loglog(bins.centers, bin_av, '-*')
        ylim(ymin, ymax)
        xlabel(xlabel_text)
        ylabel('Average edge weight ' + weight_unit)
        if (saveplots):
            savefig(out_dir + 'edge_weight_' + save_plot_name + '.eps')

    if (do_intensity or do_all):
        bins = Bins(int, min(n_nodes), max(n_nodes), 'log', 1.6)
        seq = zip(n_nodes, intensity)
        bin_av = bins.bin_average(seq)
        figure() #intensity
        subplot(121)
        loglog(n_nodes, intensity, '*')
        xlabel(xlabel_text)
        ylabel('Average intensity')
        ymin, ymax = ylim()
        subplot(122)
        loglog(bins.centers, bin_av, '-*')
        ylim(ymin, ymax)
        xlabel(xlabel_text)
        ylabel('Average intensity')
        if (saveplots):
            savefig(out_dir + 'edge_weight_' + save_plot_name + '.eps')

    if (do_coherence or do_all):
        bins = Bins(int, min(n_nodes), max(n_nodes), 'log', 1.6)
        seq = zip(n_nodes, coherence)
        bin_av = bins.bin_average(seq)
        figure() #coherence
        subplot(121)
        loglog(n_nodes, coherence, '*')
        xlabel(xlabel_text)
        ylabel('Average coherence')
        ymin, ymax = ylim()
        subplot(122)
        loglog(bins.centers, bin_av, '-*')
        ylim(ymin, ymax)
        xlabel(xlabel_text)
        ylabel('Average coherence')
        if (saveplots):
            savefig(out_dir + 'coherence_' + save_plot_name + '.eps')

    if (do_sh_pl or do_all):
        bins = Bins(int, min(n_nodes), max(n_nodes), 'log', 1.6)
        seq = zip(n_nodes, shortest_pl)
        bin_av = bins.bin_average(seq)
        figure() #shortest path length
        subplot(121)
        semilogx(n_nodes, shortest_pl, '*')
        xlabel(xlabel_text)
        ylabel('Average shortest path length')
        ymin, ymax = ylim()
        subplot(122)
        semilogx(bins.centers, bin_av, '-*')
        ylim(ymin, ymax)
        xlabel(xlabel_text)
        ylabel('Average shortest path length')
        if (saveplots):
            savefig(out_dir + 'shortest_pl_' + save_plot_name + '.eps')

    if (do_diam or do_all):
        bins = Bins(int, min(n_nodes), max(n_nodes), 'log', 1.6)
        seq = zip(n_nodes, diameter)
        bin_av = bins.bin_average(seq)
        figure() #diameter
        subplot(121)
        semilogx(n_nodes, diameter, '*')
        xlabel(xlabel_text)
        ylabel('Average diameter')
        ylim(ceil(min(diameter) -  1), floor(max(diameter) + 1))
        ymin, ymax = ylim()
        subplot(122)
        semilogx(bins.centers, bin_av, '-*')
        ylim(ymin, ymax)
        xlabel(xlabel_text)
        ylabel('Average diameter')
        if (saveplots):
            savefig(out_dir + 'diameter_' + save_plot_name + '.eps')

    if (do_show):
        show()

if __name__ == "__main__":
   main()
