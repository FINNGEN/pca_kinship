"""
This script takes two input files (the full net, and the communities)
and one output directory as arguments, and gives out a file with
average values of different properties of communities as a function of
com. size in the following format:
N_NODES N_COMMUNITIES EDGE_DENSITY EDGE_WEIGHT INTENSITY COHERENCE
SHORTEST_PATH_LENGTH DIAMETER

To get a description of the input arguments on screen, run the script
without args or see below:

python com_struct_topol_analysis.py

NOTES:
All input files and dirs need absolute path names or be in the
working directory.
Runtime on the mobile phone network, N_nodes > 10^6 was 15 mins.
8 of them to read in the full net.
Output file name is on the format:
{full net's edge file name}_average_community_analysis.txt

EXAMPLE OF USE:
python com_struct_topol_analysis.py 4-3_cliques_links_nonunit_w.edg
   4-3_cliques_links.txt_cliques_3_0_0_communities_0
   /home/eiesland/networks/Jon/communities/python_code/4-3_testdir/

Jon Eiesland
May 2009
"""

from netpython import * # lce-made net library found in
                        # networks/EDEN/nwtool/netpython

from sys import argv
import numpy as np

def main():
    if len(argv) < 4:
        print "Please give arguments. \n \
        1 - full net's edgefile (source dest weight)\n \
        2 - file with nodes in each community from Jussi Kumpula's clique percolation code (all the nodes in one community on one line) \n \
        3 - directory for analysis output (ending with '/') \n "
        exit()

    # Args
    full_net_edg = argv[1]
    comp_nodes_list = argv[2]
    out_dir = argv[3]

    # Open full net, arg 1
    print 'Beginning to read full net (Phone net, N_nodes 10^6 took 8 mins) ...'
    net=netio.loadNet(full_net_edg)
    print 'Full net read. Reading communities, doing analysis and writing to file...'

    # Read component nodes list, arg 2
    f = open(comp_nodes_list)
    communities = f.readlines()
    f.close()

    # make: n_nodes <> n_comommunities from arg 2
    n_nodessort = [] #temp. containers
    n_comsort = []
    for line in communities:
        com_nodes = map(int, line.split())
        N = len(com_nodes)
        # Insert into corresponding indices
        if (N not in n_nodessort):
            n_nodessort.append(N)
            n_comsort.append(1)
        else:
            n_comsort[n_nodessort.index(N)] += 1

    # Sort and update
    n_nodes = sorted(n_nodesort) # Final containers
    n_com = [1 for e in n_comsort]

    for es in n_nodes:
        for ens in n_nodessort:
            if (es == ens):
                i_s = n_nodes.index(es)
                i_ns = n_nodessort.index(ens)
                n_com[i_s] = n_comsort[i_ns]
                break


    # For each community, get the subnet and do the analysis
    cont_length = len(n_nodes) # Number of different com. sizes
    cont_rho = np.zeros(cont_length, float) # Containers of the different averages for the com. sizes
    cont_w = np.zeros(cont_length, float)
    cont_q = np.zeros(cont_length, float)
    cont_I = np.zeros(cont_length, float)
    cont_pl = np.zeros(cont_length, float)
    cont_d = np.zeros(cont_length, float)

    for line in communities:
        com_nodes = map(int, line.split()) #list of the nodes in this community
        N = len(com_nodes)     #number of nodes in com                
        com_net = netext.getSubnet(net, com_nodes) #get an instance of type 'net' of the community
        com_net_edges = com_net.edges #get the edges inside the community
        E = len(com_net_edges) #number of edges in com
        index = n_nodes.index(N) #get the index in container to input averages for this com. also used to get the number of communities with N nodes from n_com[index]

        # Edge density
        rho = 2. * E / N / (N-1)
        cont_rho[index] += rho / n_com[index]

        # Edge weight
        w = 0
        for e in com_net_edges:
            w += e[2] #source dest weight
        cont_w[index] += w / E / n_com[index]

        # Intensity
        I = 1
        for e in com_net_edges:
            I *= pow(e[2], 1./E)
        cont_I[index] += I / n_com[index]

        # Shortest path length and diameter
        node_shortest_pl = [] #container for shortest path length for each node
        diameter = 0.
        for e in com_nodes: #do for each node in community
            pl_dict = netext.getPathLengths(com_net, e) #for all other nodes in com get dictionary: {other node: path length from active node}
            pl_list = [pl_dict[x] for x in pl_dict] #list of the path lengths
            av_shortest = float(sum(pl_list)) / (N - 1) # N-1 paths for each node
            node_shortest_pl.append(av_shortest)
            if (diameter < max(pl_list)): #check if longest shortest pl increases
                diameter = max(pl_list)
        av_sh_pl = sum(node_shortest_pl) / N
        cont_pl[index] += av_sh_pl / n_com[index]
        cont_d[index] += float(diameter) / n_com[index]

    # Coherence
    cont_q = cont_I/cont_w

    #Print to file
    out_file_name = full_net_edg.split('/')[-1]
    out_file_name = out_file_name.split('.')[0]
    all_outfile = out_dir + out_file_name + '_average_community_analysis.txt'
    o = open(all_outfile, 'w')
    #header
    headers = ('N_NODES', 'N_COMMUNITIES', 'EDGE_DENSITY'.rjust(14),
               'EDGE_WEIGHT'.rjust(14), 'INTENSITY'.rjust(14),
               'COHERENCE'.rjust(14), 'SHORTEST_PL'.rjust(14), 'DIAMETER')
    o.write("t".join(headers) + "\n")

    for i in range(cont_length):
        data = (str(n_nodes[i]).rjust(7), str(n_com[i]).rjust(13),
                str(cont_rho[i]).rjust(14), str(cont_w[i]).rjust(14),
                str(cont_I[i]).rjust(14), str(cont_q[i]).rjust(14),
                str(cont_pl[i]).rjust(14), str(cont_d[i]).rjust(8))
        o.write("\t".join(data) + "\n")


    #For getting what's printed to file on screen, uncomment ->
    ## for i in range(len(n_nodes)):
    ##     print n_nodes[i], n_com[i], cont_rho[i], cont_w[i], cont_I[i], cont_q[i], cont_pl[i], cont_d[i]

if __name__ == '__main__':
    main()
