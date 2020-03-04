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
from pca.color_dict import color_dict
from utils import make_sure_path_exists,identify_separator,return_header,basic_iterator
from collections import Counter

import conda
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
from mpl_toolkits.basemap import Basemap

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
    if not args.cov or os.path.isfile(save_path):
        return
    
  
    region_plot_data = os.path.join(args.plot_path,'plot_data')    

    pc_avg = os.path.join(args.plot_path,'plot_data','pc_averages.csv')
    if not os.path.isfile(pc_avg):

        # read region of births from samples and get counter
        birth_df = pd.read_csv(args.cov,sep = identify_separator(args.cov), index_col='FINNGENID')[['regionofbirthname']]
        count_df = birth_df['regionofbirthname'].value_counts().to_frame('count')
        count_df.index.names = ['regionofbirth']

        # return valid regions
        loc_df = get_loc_df(args)
        print(loc_df)
        
        usecols = ['PC' + str(pc) for pc in pc_list]
        eigenvec = pd.read_csv(args.eigenvec, sep = '\t',index_col = 'IID')[usecols]
        eigenvec.index.names = ['FINNGENID']
        print(eigenvec.head())
        print(birth_df.head())
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
        legend1 = plt.legend(line1,[str(a)+ ' samples' for a in popList],loc = 'lower right',fontsize = 4)

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
    cohort_plot =  os.path.join(args.plot_path,args.name+'_cohorts_plot.pdf')

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
        eur_outliers = np.loadtxt(os.path.join(args.pca_outlier_path, 'finngen_eur_outliers.txt'),dtype = str,usecols =1)

        for tag in ['eur','fin','finngen']:
            score_file = args.eur_outlier_path + tag + '.sscore' 
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
        
    final_tags = set(pc_data['TAG']) 
    print(final_tags)
    return final_tags,pc_data
    
def plot_fin_eur_outliers(args):
    '''
    Plots eur/fin outlier detection results
    '''
    
    outlier_plot_data = os.path.join(args.plot_path,'plot_data')
    make_sure_path_exists(outlier_plot_data)

    
    outliers_plot = os.path.join(args.plot_path,args.name+'_eur_outliers.pdf')
    outliers_2d = os.path.join(args.plot_path,args.name +'_eur_outliers_pairwise.pdf')
    if not os.path.isfile(outliers_plot) or not os.path.isfile(outliers_2d):
        tags,pc_data = return_fin_eur_df(args)
            
        alpha_map = {'inliers':0.1,'EUR':1,'FIN':1,'outliers':0.1}
        size_map = {'inliers':0.1,'EUR':3,'FIN':3,'outliers':1}
        color_map = {'inliers':"red",'EUR':"blue",'FIN':"purple",'outliers':"green"}
        print('ready to plot')

    if not os.path.isfile(outliers_plot):
        plot_3d(pc_data,outliers_plot,tags,alpha_map= alpha_map,size_map = size_map,color_map = color_map,tag_column="TAG")
    else:
        args.v_print(3,'eur outliers 3d plot already done.')
    
    if not os.path.isfile(outliers_2d):
        plot_2d(pc_data,outliers_2d,tags,alpha_map= alpha_map,size_map = size_map,color_map = color_map,tag_column="TAG")
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
        eigenvec_path = os.path.join(args.pca_outlier_path,'first_round','first_round_' + args.name + '.eigenvec')
        #import metadata about samples
        outlier_info = os.path.join(args.pca_outlier_path,'first_round','first_round_' + args.name + '_outlier_samples.tsv')
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
    final_tags = set(pc_data['TAG']) 
    print(final_tags)
    return sorted(final_tags),pc_data
           
    return
    
    

def plot_first_round_outliers(args):
    '''
    Plots outlier detection results.
    '''
           
    ethnic_plot = os.path.join(args.plot_path,args.name +'_ethnic_outliers.pdf')
    ethnic_2d = os.path.join(args.plot_path,args.name +'_ethnic_outliers_pairwise.pdf')

    # create data
    if not os.path.isfile(ethnic_plot) or not os.path.isfile(ethnic_2d):
        tags,pc_data = return_outliers_df(args)
       
    #plot 3d
    if not os.path.isfile(ethnic_plot):
        #build super pop dict for plotting
        plot_3d(pc_data,ethnic_plot,tags,tag_column="TAG")
        
    else:
        args.v_print(3,'ethnic outliers 3d plot already done.')
    
    if not os.path.isfile(ethnic_2d):
        plot_2d(pc_data,ethnic_2d,tags,tag_column="TAG")
    else:
        args.v_print(3,'ethnic outliers pairwise plot already done.')
        
#######################
#--PLOTTING TEMPLATE--#
#######################
def plot_3d(pc_data,out_file,tags,pc_columns = ['PC1','PC2','PC3'],pc_tags = None,color_map= None,alpha_map = None,size_map = None,legend_fontsize = 4,label_fontsize = 5,random_samples = 5000,tag_column="COHORT"):
    
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
    ax = Axes3D(fig, azim=-135, elev=25)

    # init tag sizes
    sizes = dd(lambda :1)
    if size_map:
        for key in size_map:sizes[key] = size_map[key]
    alphas = dd(lambda:0.2)
    if alpha_map:
        for key in alpha_map:alphas[key] = alpha_map[key]
    #tag colors
    
    if not color_map:
        color_maps = list(color_dict[len(tags)]['qualitative'].keys())
        if 'Set1' in color_maps:
            cm = 'Set1'
        else:
            cm = color_maps[0]
        color_map = {tag:color_dict[len(tags)]['qualitative'][cm][i] for i,tag in enumerate(tags)}

    for i,tag in enumerate(tags):
        tag_data = pc_data[pc_data[tag_column] == tag]
        tag_data = tag_data.head(n=3000)
        print(tag,len(tag_data))
        color = color_map[tag]
        ax.scatter(tag_data[pc_columns[0]],tag_data[pc_columns[1]],tag_data[pc_columns[2]], s= sizes[tag],alpha = alphas[tag],color = color,label = tag)

    
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




def plot_2d(pc_data,out_file,tags,pc_columns = ['PC1','PC2','PC3'],pc_tags = None,color_map= None,alpha_map = None,size_map = None,legend_fontsize = 4,label_fontsize = 5,tag_column = "COHORT"):
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
    sizes = dd(lambda :1)
    if size_map:
        for key in size_map:sizes[key] = size_map[key]
    alphas = dd(lambda:0.2)
    if alpha_map:
        for key in alpha_map:alphas[key] = alpha_map[key]
    #tag colors
    if not color_map:
        color_maps = list(color_dict[len(tags)]['qualitative'].keys())
        if 'Set1' in color_maps:
            cm = 'Set1'
        else:
            cm = color_maps[0]
        color_map = {tag:color_dict[len(tags)]['qualitative'][cm][i] for i,tag in enumerate(tags)}

    for i,tag in enumerate(tags):
        tag_data = pc_data[pc_data[tag_column] == tag]
        tag_data = tag_data.head(n=3000)
        color = color_map[tag]
        for i,pcs in enumerate(list(combinations(pc_columns,2)) ):
            ax = axes[i]
            plot_x,plot_y = tag_data[pcs[0]],tag_data[pcs[1]]
            ax.scatter(plot_x,plot_y, s= sizes[tag],alpha = alphas[tag],color = color,label = tag)

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
            
    leg = ax1.legend(loc='upper left', numpoints=1, fancybox = True,prop={'size': legend_fontsize})
    for lh in leg.legendHandles: 
        lh.set_alpha(1)
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
