from collections import defaultdict as dd
from collections import Counter
from verkko.binner import binner
import numpy as np
import pandas as pd
import os
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(palette='Set2')
import matplotlib.ticker as ticker
import networkx as nx
from pca_scripts.color_dict import color_dict
from utils import identify_separator,return_header,basic_iterator,mapcount
from matplotlib.lines import Line2D

inf_dict = {}
inf_types = ['Dup/MZ','PO','FS','2nd','3rd']
for i,key in enumerate(inf_types):inf_dict[key] = i

def plot_degree_dist(args,filter_batches = None):
    """
    Plots degree distribution for each inftype
    """
    
    args.degree_fig = os.path.join(args.out_path,args.prefix +'_degree_distribution.pdf')
    if os.path.isfile(args.degree_fig) and not args.force:
        print(f"{args.degree_fig} already generated")
        return

   

    # bins. np.inf added to keep length static across
    bins = [1.5 + i for i in range(10)] 
    xlabels = [str(int(elem)) for elem in bins]
    bins += [np.inf]
    xlabels = np.array(['0'] +  xlabels + ['11+'])
    
    save_count_data =  os.path.join(args.misc_path,args.prefix +'_degree_count.npy')
    save_deg_data =  os.path.join(args.misc_path,args.prefix +'_degree_data.npy')

    if os.path.isfile(save_count_data) and os.path.isfile(save_deg_data):
        deg_count = np.load(save_count_data)
        deg_data = np.load(save_deg_data)
    else:
        G  = read_network(args)

        # TEST FILTERING OUT LARGE FAMILIES
        if filter_batches:
            sample_batch_dict = read_batch_data(args)
            for batch in filter_batches:
                nodes = [sample for sample in sample_batch_dict if sample_batch_dict[sample] == batch]
                print(len(nodes),batch)
                G.remove_nodes_from(nodes)   

        #CREATE ARRAY OF COUNTS
        deg_count = np.zeros(len(bins)+1)
        g_deg = [d for n,d in G.degree()]
        total_samples = mapcount(args.kinship_bed.replace('.bed','.fam'))# total samples in plink data
        deg_count[0] =total_samples - G.number_of_nodes() 
        # assign degree count to bins adding an empty bin to fix indexes
        count = Counter(np.digitize(g_deg,[0.5] + bins))
        for x in count.keys():deg_count[x] = count[x]
        np.save(save_count_data,deg_count)

        # DEGREE EDGE DATA
        deg_data =[]
        # loop through relevant bins 
        for i,bin in enumerate(bins):
            print(i,bin)
            # array that contains counts for each edge weight
            tmp_data = np.zeros(len(inf_types))
            # nodes of degree in bin i
            deg_nodes = [node for node in G.nodes if np.digitize(G.degree[node],bins) == i]
            print(len(deg_nodes),deg_count[i+1]) # should match
            deg_weights = []
            for node in deg_nodes:
                for n2 in G[node]:
                    deg_weights.append(G[node][n2]['weight'])
            count = Counter(deg_weights)
            for key in count:tmp_data[key] = count[key]
            print(tmp_data,np.sum(tmp_data))
            deg_data.append(tmp_data)

        deg_data = np.array(deg_data)
        np.save(save_deg_data, deg_data)

    x_data = np.arange(len(deg_count))
    print(x_data,deg_count.shape,deg_data.shape)  

    # mask missing data
    data_mask = (deg_count > 0)
    x_data,deg_count,deg_data = x_data[data_mask],deg_count[data_mask],deg_data[data_mask[1:]]
    deg_count = np.log10(deg_count)
    
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0,0])      

    # PLOT DONUTS
    radius = 0.4
    width = radius*0.6
    colors = color_dict[len(inf_types)]['diverging']['Spectral']
    for i,coord in enumerate(zip(x_data[1:],deg_count[1:])):
        print(coord)
        print(deg_data[i])
        wedges, *_ = ax.pie(deg_data[i],colors = colors, center=coord, radius=radius,wedgeprops={'width':width,'linewidth': 0},frame=True)
    ax.legend(wedges, list(inf_dict.keys()),loc='lower left', prop={'size': 6})

    # plot global curve
    ax.scatter(x_data,deg_count,c = 'k')
    # X AXIS STUFF
    ax.set_xlabel('Number of relatives per participant')
    ax.set_xlim((min(x_data)-1,max(x_data)+1))
    ax.set_xticks(x_data)
    ax.set_xticklabels(xlabels[data_mask])

    # y AXIS STUFF
    ax.set_ylim((1,1+np.max(deg_count)))   
    ax.set_ylabel(r'Number of participants')
    yticks = list(range(0,int(max(deg_count)) +2))
    ax.set_yticks(yticks)
    ylabels = [r"$10^{{ {:2d} }}$".format(exponent) for exponent in yticks]
    ax.set_yticklabels(ylabels)

    for tick in ax.xaxis.get_major_ticks() + ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(6)
        
    fig.savefig(args.degree_fig)
    plt.close()
    print('done')

        
def plot_kinship(args):
    '''
    Plots the ditribution of kinship values.
    '''

    args.kinship_fig = os.path.join(args.out_path,args.prefix +'_kinship_distribution.pdf')
    if  os.path.isfile(args.kinship_fig) and not args.force:
        print(3,'kinship plot already done')
        return
    
    else:
       print('plotting kinship...')
   
    kin_dump = os.path.join(args.misc_path, args.prefix +'_kinship.npy')
    if not os.path.isfile(kin_dump) or  args.force: 
        print('uploading kinship data')
        kin_data = pd.read_csv(args.kin_file,sep = identify_separator(args.kin_file), usecols = ['Kinship']).values.flatten()
        print('dumping it')
        np.save(kin_dump,kin_data)
    else:

        kin_data = np.load(kin_dump, allow_pickle = True)
        
    entries = len(kin_data)

    xPos = [0.0442,0.0884,0.177,0.354,0.51]
       
    bin_data_file = os.path.join(args.kinship_path,'bin_average.npy')
    plot_data_file = os.path.join(args.kinship_path,'plot_data.npy')
    print('importing bin data..')
    if not os.path.isfile(bin_data_file) or not os.path.isfile(plot_data_file) or args.force:
        bins = binner.Bins(float,min(xPos),1,'log',1.05)
        print('bin data missing, generating...')
        countNotNormalized = bins.bin_count_divide(kin_data)
        count = np.array(binner.normalize(list(countNotNormalized)))
        binAvg = bins.bin_average(zip(kin_data,kin_data))
        binMask = ~np.ma.getmask(binAvg)
        plot_data = count[binMask]
        bin_data = binAvg[binMask]
        bin_data.dump(bin_data_file)
        plot_data.dump(plot_data_file)
    else:
        bin_data = np.load(bin_data_file,allow_pickle = True)
        plot_data = np.load(plot_data_file,allow_pickle = True)
        
    print('plotting...')

    
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(1,1)
   

    texts = ['3rd','2nd','FS/PO','MZ/Duplicates']
    ax = fig.add_subplot(gs[0,0])
    ax2 = ax.twinx()
    ax2.set_ylim([0,1])
    cum_data = np.cumsum(plot_data)
    ax.plot(bin_data,plot_data, '--')
    ax2.plot(bin_data,cum_data, color = 'red', linewidth = 1)

    ax.set_xticks(xPos[:-1])
    ax.set_xlabel('kinship')
    ax.set_ylabel(r'P(k)')

    ax2.set_yticks([0,1])
    #ax2.set_yticklabels([0,entries])
    ax2.set_ylabel(r'$P(k < x) $')

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation(30)
        tick.label.set_fontsize(7)     

    fig.savefig(args.kinship_fig)
    plt.close()
    print('done')

    
def plot_batch_data(args):
    """
    Plots some network properties of the kinship networks of different batches.
    """
    
    args.batch_fig = os.path.join(args.out_path,args.prefix +'_batches.pdf')
    if os.path.isfile(args.batch_fig):
        print(args.batch_fig + " already generated")
        return

    plot_data = os.path.join(args.misc_path, args.prefix +'_batches_degree.npy')
    if not os.path.isfile(plot_data):
        G = read_network(args)
        sample_batch_dict = read_batch_data(args)
        batches = sorted(set(sample_batch_dict.values()))
        print(batches)
        network_batch_data = np.zeros((len(batches),3))
        for i,batch in enumerate(batches):
            batch_nodes = [s for s in sample_batch_dict if sample_batch_dict[s] == batch]
            z = nx.subgraph(G,batch_nodes)
            print(batch,z.number_of_nodes(),z.number_of_edges())
            avg_deg = 2*z.number_of_edges()/z.number_of_nodes()
            clustering = nx.average_clustering(z)
            network_batch_data[i] = (avg_deg,clustering,z.number_of_nodes())
        np.save(plot_data,network_batch_data)
    else:
        network_batch_data = np.load(plot_data)    
    
    print(network_batch_data)
    batches = np.loadtxt(os.path.join(args.misc_path,'cohorts.txt'),dtype = str)
    n = int(np.average(network_batch_data[:,2]))
    max_avg_deg = network_batch_data[:,0].max()
    print(n,max_avg_deg)
    def markers_iterator():
        for marker in Line2D.filled_markers:yield marker

    mit = markers_iterator()  
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(1,1)
    ax = fig.add_subplot(gs[0,0])

    for i,elem in enumerate(network_batch_data):
        x,y,s = elem
        batch = batches[i]
        axiom = True if batch in [f"R{i}" for i in range(2,20)] else False
        if axiom:
            ax.scatter(x,y,s/1000,label = batch,color = 'black',marker = next(mit),alpha = 0.5)
        else:
            ax.scatter(x,y,s/1000,label = batch)


    # RANDOM GRAPH
    print(max_avg_deg)
    p_max = max_avg_deg/(n-1)
    print('pmax:',p_max)
    random_data = []
    for p in np.linspace(0,p_max,100):
        avg_deg = p*(n-1)
        random_data.append([avg_deg,p])
    data = np.array(random_data)
    ax.plot(data[:,0],data[:,1],'k--',label = 'random graph',linewidth = 0.5)   
  
    ax.set_xlabel('Average degree')
    ax.set_ylabel('Average clustering coefficient')
    ax.legend(loc='upper left', prop={'size': 6})

    
    fig.savefig(args.batch_fig)
    plt.close()
    print('done')
    
def read_batch_data(args):

    sample_info = pd.read_csv(args.sample_info,sep = identify_separator(args.sample_info))
    sample_batch_dict = pd.Series(sample_info.COHORT.values,index=sample_info.IID).to_dict()
    return sample_batch_dict


def read_network(args):

    edgelist = os.path.join(args.misc_path,'kinship.edgelist')
    if not os.path.isfile(edgelist):
        header = return_header(args.kin_file)
        idx = [header.index(elem) for elem in ['ID1','ID2','InfType']]
        deg_iterator = basic_iterator(args.kin_file,columns = idx,skiprows=1)
        G = nx.Graph()
        for id1,id2,inf in deg_iterator:
            G.add_edge(id1,id2,weight = inf_dict[inf])
        nx.write_weighted_edgelist(G,edgelist)

    G = nx.read_edgelist(edgelist,data=(('weight',int),))
    return G
