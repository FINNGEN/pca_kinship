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
import os,unidecode
from itertools import combinations
from pca_scripts.color_dict import color_dict
from utils import make_sure_path_exists,identify_separator,mapcount
from collections import Counter

import conda
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
#from mpl_toolkits.basemap import Basemap



#########################
#--PLOTTING REGION PCA--#
#########################

def get_loc_df(args):

    with open(os.path.join(args.data_path,'regionlocation.txt')) as i: data = [elem.strip() for elem in i.readlines()[1:]]
    #data = [elem.strip() for elem in data[1:]]
    a = {}
    for entry in data:
        region,_,*coord = entry.split(',')
        a[region] = list(map(float,coord))
    loc_df = pd.DataFrame.from_dict(a,orient='index',columns=['lat','lon'])
    loc_df.index.names = ['regionofbirth']
    return loc_df

def plot_map(args,pc_list = [1,2,3]):

    save_path = os.path.join(args.plot_path, args.name + f"_pc_map.pdf")
    if not args.meta or os.path.isfile(save_path):
        return

    region_plot_data = os.path.join(args.plot_path,'plot_data')

    pc_avg = os.path.join(args.plot_path,'plot_data','pc_averages.csv')
    if not os.path.isfile(pc_avg):

        # read region of births from samples and get counter
        birth_df = pd.read_csv(args.meta,sep = identify_separator(args.meta), index_col='FINNGENID')[['regionofbirthname']]

        # return valid regions
        loc_df = get_loc_df(args)
        print(loc_df)

        usecols = ['PC' + str(pc) for pc in pc_list]
        eigenvec = pd.read_csv(args.eigenvec, sep = '\t',index_col = 'IID')[usecols]
        eigenvec.index.names = ['FINNGENID']
        print(eigenvec.head())
        print(birth_df.head())

        count_df =  pd.concat([birth_df,eigenvec],join = 'inner',axis = 1)['regionofbirthname'].value_counts().to_frame('count')
        count_df.index.names = ['regionofbirth']

        tmp_avg = pd.concat([birth_df,eigenvec],join = 'inner',axis = 1).set_index('regionofbirthname').groupby('regionofbirthname').mean()
        print(tmp_avg.head())
        region_avg = pd.concat([loc_df,tmp_avg,count_df],join = 'inner',axis = 1)
        region_avg.index.names = ['regionofbirth']
        print(region_avg)

        region_avg.to_csv(pc_avg,index=True)
        print(region_avg.shape)
    df = pd.read_csv(pc_avg,index_col = 0)
    print(df)

    resize = 100/max(df['count'].values)
    m = Basemap(projection='lcc', resolution= 'l', lat_0=np.average(df['lat']), lon_0=np.average(df['lon']), width=1E6, height=1.2E6)
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(1,len(pc_list))
    for i,pc in enumerate(list(map(str,pc_list))):
        ax = fig.add_subplot(gs[0,i])
        m.scatter(df['lon'].values,df['lat'].values, latlon=True,s=df['count'].values*resize,c=df['PC' + str(pc)].values,cmap='seismic', alpha=0.8, edgecolors = 'k', linewidth = .5,zorder= 10)

        #make legend with dummy points
        line1=[]
        popList =[100, 1000, 10000]
        for a in popList:
            locLine= plt.scatter([], [], c='k', alpha=0.5, s= a*resize, label=str(a) + ' samples')
            line1.append(locLine)
        legend1 = plt.legend(line1,[str(a)+ ' samples' for a in popList],loc = 'lower right',
        fontsize = 6)

        m.shadedrelief()
        m.drawcoastlines(color='gray',linewidth = 0.3)
        m.drawcountries(color='gray')
    fig.savefig(save_path, bbox_inches = 'tight', pad_inches = 0)



########################
#--PLOTTING FINAL PCA--#
########################
def plot_final_pca(args):

    outlier_plot_data = os.path.join(args.plot_path,'plot_data')
    make_sure_path_exists(outlier_plot_data)

    pca_pairwise =os.path.join(args.plot_path,args.name+ '_final_pca_pairwise.pdf')
    pca_plot =  os.path.join(args.plot_path,args.name+'_final_pca.pdf')
    cohort_plot =  os.path.join(args.plot_path,args.name+'_cohorts.pdf')

    if not np.all(list(map(os.path.isfile,[pca_pairwise,pca_plot,cohort_plot]))):
        print('loading pc data...')
        tags,pc_data = return_cohorts_df(args)
        print('done')

    if not os.path.isfile(pca_plot):
        plot_3d(pc_data,pca_plot,tags)

    pcs = ["PC1",'PC2','PC3']
    if not os.path.isfile(cohort_plot):
        print(cohort_plot)
        plot_cohort_averages(pc_data,tags,cohort_plot,pcs)

    if not os.path.isfile(pca_pairwise):
        plot_2d(pc_data,pca_pairwise,tags)


def return_cohorts_df(args):
    """
    Returns a pandas df where the cohort info is in a column.
    N.B. the eigenvec dataframe and the batches dataframe are merged on IID, so if the samples is missing in the batches metadata, it will not be plotted.
    """

    out_file = os.path.join(args.plot_path,'plot_data',"pc_cohorts.csv")
    if not os.path.isfile(out_file):
        pc_data = pd.read_csv(args.eigenvec,sep = '\t',usecols = ['IID',"PC1",'PC2','PC3'], dtype = {pc: np.float64 for pc in ["PC1",'PC2','PC3']})
        print(len(pc_data))
        cohort_data = pd.read_csv(args.sample_info)

        # mapping smallest batchets to other if too many
        while len(set(cohort_data['COHORT'].values))  > max(color_dict.keys()):
            count = Counter(cohort_data['COHORT'].values)
            smallest_cohort = min(count,key = count.get)
            cohort_data.loc[cohort_data['COHORT'] == smallest_cohort,'COHORT'] = 'Other'

        print(len(cohort_data))
        pc_data = pc_data.merge(cohort_data,on = "IID")
        pc_data.to_csv(out_file,index=False)

    else:
        pc_data = pd.read_csv(out_file)

    final_cohorts = set(pc_data['COHORT'])
    print(final_cohorts)
    return sorted(final_cohorts),pc_data

def plot_cohort_averages(pc_data,cohorts,out_file,pc_columns):

    print('plotting cohorts...')
    fig, axes = plt.subplots(nrows=len(pc_columns), sharex=True)

    for i,column in enumerate(pc_columns):
        ax = axes[i]
        ax.set_ylabel(column)
        violin_data = []
        for i,cohort in enumerate(cohorts):
            print(cohort)
            cohort_data = pc_data[column][pc_data["COHORT"] == cohort]
            violin_data.append(cohort_data)

        sns.violinplot(data=violin_data,ax = ax,scale = 'count')
        ax.set(xticklabels=cohorts)
        ax.tick_params(axis='both', which='major', labelsize=6, rotation=45)

    fig.savefig(out_file)
    plt.close()

##########################
#--PLOTTING FIN EUR PCA--#
##########################

def return_fin_eur_df(args):
    """
    Returns pandas dataframe with eur/fin/inlier/outlier data.
    """
    out_file = os.path.join(args.plot_path,'plot_data',"pc_eur.csv")

    if not os.path.isfile(out_file):
        pc_avg = ["PC1_AVG",'PC2_AVG','PC3_AVG']
        df_list = []
        eur_outliers = np.loadtxt(args.finngen_eur_outliers,dtype = str)
        for tag in ['eur','fin','finngen']:
            score_file =  os.path.join(args.pca_outlier_path, "eur_pca/", tag + '.sscore' )
            df =  pd.read_csv(score_file,sep = '\t',usecols = ['IID'] + pc_avg, dtype = {pc: np.float64 for pc in pc_avg}).rename(columns = {pc: pc.replace("_AVG","") for pc in pc_avg})
            df["TAG"] = tag.upper()
            if tag =='finngen':
                df['TAG'] = np.where(df['IID'].isin(eur_outliers), "outliers", "inliers")
            df_list.append(df)

        pc_data = pd.concat(df_list)
        print(pc_data.head())
        pc_data.to_csv(out_file,index=False)

    else:
        pc_data = pd.read_csv(out_file)


    final_tags = set(pc_data.TAG)
    final_tags = ['inliers','FIN','outliers','EUR']
    print(final_tags)
    return final_tags,pc_data

def plot_fin_eur_outliers(args):
    '''
    Plots eur/fin outlier detection results
    '''

    if mapcount(args.finngen_eur_outliers) == 0:
        return
    outlier_plot_data = os.path.join(args.plot_path,'plot_data')
    make_sure_path_exists(outlier_plot_data)

    outliers_plot = os.path.join(args.plot_path,args.name+'_eur_outliers.pdf')
    outliers_2d = os.path.join(args.plot_path,args.name +'_eur_outliers_pairwise.pdf')
    if not os.path.isfile(outliers_plot) or not os.path.isfile(outliers_2d):
        tags,pc_data = return_fin_eur_df(args)

        alpha_map = {'inliers':0.1,'EUR':1,'FIN':1,'outliers':0.1}
        size_map = {'inliers':1,'EUR':3,'FIN':3,'outliers':1}
        colors= color_dict[len(tags)]['qualitative']['Set1']
        [red,blue,green,purple] = colors
        color_map = {'inliers':red,'EUR':blue,'FIN':purple,'outliers':green}
        print('ready to plot')

    if not os.path.isfile(outliers_plot):
        plot_3d(pc_data,outliers_plot,tags,alpha_map= alpha_map,size_map = size_map,color_map = color_map,tag_column="TAG")
    else:
        args.v_print(3,'eur outliers 3d plot already done.')

    if not os.path.isfile(outliers_2d):
        plot_2d(pc_data,outliers_2d,tags,alpha_map= alpha_map,size_map = size_map,color_map = color_map,tag_column="TAG", legend_location = "lower right",legend_fontsize=10)
    else:
        args.v_print(3,'eur outliers pairwise plot already done.')



#########################
#--PLOTTING ETHNIC PCA--#
#########################

def return_outliers_df(args):
    """
    Returns pandas df with ethnic outlier info
    """
    out_file = os.path.join(args.plot_path,'plot_data',"pc_ethnic.csv")
    if not os.path.isfile(out_file):
        tg_pca_file=  os.path.join(args.pca_outlier_path, '1k_pca/',args.name)
        eigenvec_path =  tg_pca_file  + '.eigenvec'
        #import metadata about samples
        outlier_info = tg_pca_file + '_outlier_samples.tsv'
        samples = np.loadtxt(args.sample_fam,usecols = 1,dtype = str)
        superpops = set(np.loadtxt(args.data_path + 'superpop.csv',dtype = str,delimiter =',',usecols = 1))
        # read pc data
        pc_data = pd.read_csv(eigenvec_path,sep = '\t',usecols = ['IID',"PC1",'PC2','PC3'], dtype = {pc: np.float64 for pc in ["PC1",'PC2','PC3']})
        # set finngen samples as "FINNGEN"
        outlier_data = pd.read_csv(outlier_info,dtype = str,sep = '\t',usecols = ['IID',"outlier",'SuperPops']).rename(columns ={"SuperPops":"TAG"})

        outlier_data.loc[((outlier_data['IID'].isin(samples)) & (outlier_data["outlier"] == "TRUE")),"TAG"] = "FINNGEN_OUTLIER"
        outlier_data.loc[((outlier_data['IID'].isin(samples)) & (outlier_data["outlier"] == "FALSE")),"TAG"] = "FINNGEN_INLIER"
        print(outlier_data.head())
        pc_data = pc_data.merge(outlier_data,on = "IID")
        print(pc_data.head())
        pc_data.to_csv(out_file,index=False)

    else:
        pc_data = pd.read_csv(out_file)

    print(pc_data.dtypes)
    tags = set(pc_data.TAG)
    tag_size = [len(pc_data[pc_data.TAG == tag]) for tag in tags]
    
    # plot them by size!
    final_tags = [tag for _,tag in sorted(zip(tag_size,tags),reverse = True)]
    color_maps = list(color_dict[len(tags)]['qualitative'].keys())
    cm = 'Set1' if 'Set1' in color_maps else color_maps[0] 
    colors= color_dict[len(tags)]['qualitative'][cm]
    [red,blue,green,purple,orange,yellow,brown,pink] = colors
    color_map = {'FIN':purple,'FINNGEN_INLIER':red,'FINNGEN_OUTLIER':green,'EUR':blue,'AFR':pink,'EAS':yellow,'SAS':brown,'AMR':orange}
    #color_map ={final_tags[i]:color for i,color in enumerate(colors)}

    print(color_map)
    return final_tags,pc_data,color_map


def plot_first_round_outliers(args):
    '''
    Plots outlier detection results.
    '''

    ethnic_plot = os.path.join(args.plot_path,args.name +'_ethnic_outliers.pdf')
    ethnic_2d = os.path.join(args.plot_path,args.name +'_ethnic_outliers_pairwise.pdf')

    # create data
    if not os.path.isfile(ethnic_plot) or not os.path.isfile(ethnic_2d):
        tags,pc_data,color_map = return_outliers_df(args)


    #plot 3d
    if not os.path.isfile(ethnic_plot):
        #build super pop dict for plotting
        plot_3d(pc_data,ethnic_plot,tags,tag_column="TAG",color_map = color_map)

    else:
        args.v_print(3,'ethnic outliers 3d plot already done.')

    if not os.path.isfile(ethnic_2d):
        plot_2d(pc_data,ethnic_2d,tags,tag_column="TAG",color_map=color_map)
    else:
        args.v_print(3,'ethnic outliers pairwise plot already done.')

        
#######################
#--PLOTTING TEMPLATE--#
#######################
def plot_3d(pc_data,out_file,tags,pc_columns = ['PC1','PC2','PC3'],pc_tags = None,color_map= None,alpha_map = None,size_map = None,legend_fontsize = 4,label_fontsize = 5,random_samples = 5000,tag_column="COHORT",azim=-135,elev = 25,max_size = 3000,max_map = None):

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
    #tag colors
    print(alphas)
    print(sizes)

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




def plot_2d(pc_data,out_file,tags,pc_columns = ['PC1','PC2','PC3'],pc_tags = None,color_map= None,alpha_map = None,size_map = None,legend_fontsize = 6,label_fontsize = 5,tag_column = "COHORT",max_size = 3000,max_map = None,axis_legend =2,legend_location = "lower left",rescale = 4):
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
    #tag colors
    print(alphas)
    print(sizes)

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
