import numpy as np
from collections import defaultdict as dd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(palette='Set2')
import matplotlib.ticker as ticker
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import os,unidecode,natsort
from itertools import combinations
from pca_scripts.color_dict import color_dict
from utils import make_sure_path_exists,identify_separator,mapcount,tmp_bash
from collections import Counter
import conda
# conda_file_dir = conda.__file__
# conda_dir = conda_file_dir.split('lib')[0]
# proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
# os.environ["PROJ_LIB"] = proj_lib



#MAP STUFF
import cartopy,os
import cartopy.crs as crs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature


locations = "Aland,Mariehamn,60.0971,19.9348\nCentral Finland,Jyvaskyla,62.2393,25.745951\nCentral Ostrobothnia,Kokkola,63.836583,23.133625\nEtelä-Savo,Mikkeli,61.687796,27.272657\nKainuu,Kajaani,64.224087,27.733423\nKanta-Häme,Hameenlinna,60.981091,24.458264\nKymenlaakso,Kotka,60.467423,26.945084\nLapland,Rovaniemi,66.497621,25.71921\nNorth Karelia,Joensuu,62.600816,29.760536\nNorth Ostrobothnia,Oulu,65.011873,25.471681\nOstrobothnia,Vaasa,63.081821,21.479812\nPirkanmaa,Tampere,61.498021,23.760312\nPohjois-Savo,Kuopio,62.824142,27.594561\nPäijät-Häme,Lahti,60.983876,25.656181\nSatakunta,Pori,61.486724,21.791002\nSouth Karelia,Lappeenranta,61.058242,28.18753\nSouth Ostrobothnia,Seinajoki,62.79541,22.844202\nUusimaa,Helsinki,60.167409,24.942568\nVarsinais-Suomi,Turku,60.451753,22.267052"


conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib


########################
#---PLOTTING PCA MAP---#
########################

def pca_map(args):
    loc_df = get_loc_df()
    print(loc_df)
    
def get_loc_df():
    """
    Create region to lat/lon mapping
    """

    #data = [elem.strip() for elem in data[1:]]
    a = {}
    for entry in locations.split('\n'):
        region,_,*coord = entry.split(',')
        a[region] = list(map(float,coord))
    loc_df = pd.DataFrame.from_dict(a,orient='index',columns=['lat','lon'])
    loc_df.index.names = ['regionofbirth']

    return loc_df


########################
#--PLOTTING FINAL PCA--#
########################

def plot_final_pca(args):

    outlier_plot_data = os.path.join(args.plot_path,'plot_data')
    make_sure_path_exists(outlier_plot_data)

    pca_pairwise =os.path.join(args.plot_path,args.name+ '_final_pca_pairwise.pdf')
    pca_plot =  os.path.join(args.plot_path,args.name+'_final_pca.pdf')
    tag_plot =  os.path.join(args.plot_path,args.name+'_tags.pdf')

    if not np.all(list(map(os.path.isfile,[pca_pairwise,pca_plot,tag_plot]))):
        print('loading pc data...')
        tags,pc_data = return_tags_df(args)
        print('done')
    if not os.path.isfile(pca_plot):
        plot_3d(pc_data,pca_plot,tags)

    pcs = ["PC1",'PC2','PC3']
    if not os.path.isfile(tag_plot):
        print(tag_plot)
        plot_tag_averages(pc_data,tags,tag_plot,pcs)

    if not os.path.isfile(pca_pairwise):
        plot_2d(pc_data,pca_pairwise,tags)


        
def return_tags_df(args):
    """
    Returns a pandas df where the cohort info is in a column.
    N.B. the eigenvec dataframe and the batches dataframe are merged on IID, so if the samples is missing in the batches metadata, it will not be plotted.
    """

    out_file = os.path.join(args.plot_path,'plot_data',"pc_tags.csv")
    if not os.path.isfile(out_file):
        #read in pc_data
        pc_data = pd.read_csv(args.eigenvec,sep = '\t',usecols = ['IID',"PC1",'PC2','PC3'], dtype = {pc: np.float64 for pc in ["PC1",'PC2','PC3']})
        print(len(pc_data))
        sample_data = pd.read_csv(args.sample_info,sep='\t',index_col=0,header=0,names=["IID","TAG"])
        # mapping smallest batchets to other if too many
        while len(set(sample_data['TAG'].values))  > max(color_dict.keys()):
            count = Counter(sample_data['TAG'].values)
            smallest_cohort = min(count,key = count.get)
            sample_data.loc[sample_data['zAG'] == smallest_cohort,'TAG'] = 'Other'

        
        pc_data = pc_data.merge(sample_data,on = "IID")
        print(pc_data)
        pc_data.to_csv(out_file,index=False)

    else:
        pc_data = pd.read_csv(out_file)

    final_cohorts = natsort.natsorted(set(pc_data['TAG']))
    print(final_cohorts,len(pc_data))
    return final_cohorts,pc_data


def plot_tag_averages(pc_data,tags,out_file,pc_columns):

    print('plotting tags...')
    fig, axes = plt.subplots(nrows=len(pc_columns), sharex=True)

    for i,column in enumerate(pc_columns):
        ax = axes[i]
        ax.set_ylabel(column)
        violin_data = []
        for i,tag in enumerate(tags):
            tag_data = pc_data[column][pc_data["TAG"] == tag]
            violin_data.append(tag_data)

        sns.violinplot(data=violin_data,ax = ax,scale = 'count')
        ax.set(xticklabels=tags)
        ax.tick_params(axis='both', which='major', labelsize=6, rotation=45)

    fig.savefig(out_file)
    plt.close()


       
#######################
#--PLOTTING TEMPLATE--#
#######################
def plot_3d(pc_data,out_file,tags,pc_columns = ['PC1','PC2','PC3'],pc_tags = None,color_map= None,alpha_map = None,size_map = None,legend_fontsize = 4,label_fontsize = 5,random_samples = 5000,tag_column="TAG",azim=-135,elev = 25,max_size = 3000,max_map = None):

    '''
    Inputs:
    -- pc_file : name of file where to fetch data
    -- out_file : name of figure for saving
    -- tag_column : columns of the dataframe to plot
    -- tags: list of tags to plot in tag_column
    -- pc_columns : name of the columns in the file
    -- pc_tags : name of the pcs to plt
    -- colors_map : annotation to color_dict
    '''

    if not pc_tags: pc_tags = pc_columns

    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(1,1)
    ax = Axes3D(fig, azim=azim, elev=elev)

    # init tag sizes
    sizes = dd(lambda :.5)
    if size_map:
        for key in size_map:sizes[key] = size_map[key]
    alphas = dd(lambda:0.2)
    if alpha_map:
        for key in alpha_map:alphas[key] = alpha_map[key]
    max_sizes = dd(lambda:max_size)
    if max_map:
        for key in max_map: max_sizes[key] = max_map[key]

    if not color_map:
        color_maps = list(color_dict[len(tags)]['qualitative'].keys())
        cm = 'Set1' if 'Set1' in color_maps else color_maps[0]
        color_map = {tag:color_dict[len(tags)]['qualitative'][cm][i] for i,tag in enumerate(tags)}

    for i,tag in enumerate(tags):
        tag_data = pc_data[pc_data[tag_column] == tag]
        max_data = max_sizes[tag]
        tag_data = tag_data.head(n=max_data)
        color = color_map[tag]
        size = sizes[tag]
        alpha = alphas[tag]
        print(tag,len(tag_data),alpha,size,max_data)
        ax.scatter(tag_data[pc_columns[0]],tag_data[pc_columns[1]],tag_data[pc_columns[2]], s= size,alpha = alpha,color = color,label = tag)


    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.linspace(start,end,5))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3f'))
    for ticks in [ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks(),ax.zaxis.get_major_ticks()]:
        for tick in ticks:
            tick.label.set_fontsize(label_fontsize)

    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.linspace(start,end,5))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3f'))

    start, end = ax.get_zlim()
    ax.zaxis.set_ticks(np.linspace(start,end,5))
    ax.zaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3f'))


    ax.set_xlabel(pc_tags[0],fontsize=10)
    ax.set_ylabel(pc_tags[1],fontsize=10)
    ax.set_zlabel(pc_tags[2],fontsize=10)

    trim_axis(ax)
    leg = ax.legend(loc='upper left', numpoints=1, fancybox = True,prop={'size': legend_fontsize})
    for lh in leg.legendHandles:
        lh.set_alpha(1)
        lh._sizes = [50]
        
    fig.savefig(out_file)
    plt.close()




def plot_2d(pc_data,out_file,tags,pc_columns = ['PC1','PC2','PC3'],pc_tags = None,color_map= None,alpha_map = None,size_map = None,legend_fontsize = 6,label_fontsize = 5,tag_column = "TAG",max_size = 3000,max_map = None,axis_legend =2,legend_location = "lower left",rescale = 4):
    '''
     Inputs:
    -- pc_file : name of file where to fetch data
    -- out_file : name of figure for saving
    -- tag_column : columns of the dataframe to plot
    -- tags: list of tags to plot in tag_column
    -- pc_columns : name of the columns in the file
    -- pc_tags : name of the pcs to plt
    -- colors_map : annotation to color_dict
    '''

    if not pc_tags: pc_tags = pc_columns

    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(3,1)

    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0],sharex = ax1)
    ax3 = fig.add_subplot(gs[2,0],sharex = ax1)
    axes=[ax1,ax2,ax3]

    
    # init tag sizes
    sizes = dd(lambda :.5)
    if size_map:
        for key in size_map:sizes[key] = size_map[key]
    alphas = dd(lambda:0.2)
    if alpha_map:
        for key in alpha_map:alphas[key] = alpha_map[key]
    max_sizes = dd(lambda:max_size)
    if max_map:
        for key in max_map: max_sizes[key] = max_map[key]

    if not color_map:
        color_maps = list(color_dict[len(tags)]['qualitative'].keys())
        cm = 'Set1' if 'Set1' in color_maps else color_maps[0]
        color_map = {tag:color_dict[len(tags)]['qualitative'][cm][i] for i,tag in enumerate(tags)}

    for i,tag in enumerate(tags):
        tag_data = pc_data[pc_data[tag_column] == tag]
        max_data = max_sizes[tag]
        n = min(max_data,len(tag_data))
        tag_data = tag_data.sample(n=n)
        color = color_map[tag]
        size = sizes[tag]
        alpha = alphas[tag]
        print(tag,len(tag_data),alpha,size,n)
        for i,pcs in enumerate(list(combinations(pc_columns,2)) ):
            ax = axes[i]
            plot_x,plot_y = tag_data[pcs[0]],tag_data[pcs[1]]
            ax.scatter(plot_x,plot_y, s= size,alpha = alpha,color = color,label = tag)        


    for i,tags in enumerate(list(combinations(pc_tags,2))):
        tag1,tag2 = tags
        ax = axes[i]
        ax.set_xlabel(tag1)
        ax.set_ylabel(tag2)
        trim_axis(ax)
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3f'))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.3f'))
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(6)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(6)

    leg_ax = axes[axis_legend]
    leg = leg_ax.legend(loc=legend_location, numpoints=1, fancybox = True,prop={'size': legend_fontsize})
    for lh in leg.legendHandles:
        lh.set_alpha(1)
        if rescale:
            lh._sizes = [lh._sizes[0]*7]
        else:
            lh._sizes = [50]

    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    fig.savefig(out_file)
    plt.close()



def trim_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    try:
        ax.get_yaxis().tick_left()
    except:
        pass
