import numpy as np
from collections import defaultdict as dd
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
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
import scipy.stats as st

##########################
#--- PLOTTING PCA MAP ---#
##########################


import cartopy
import cartopy.crs as crs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

locations = "Aland,Mariehamn,60.0971,19.9348\nCentral Finland,Jyvaskyla,62.2393,25.745951\nCentral Ostrobothnia,Kokkola,63.836583,23.133625\nEtelä-Savo,Mikkeli,61.687796,27.272657\nKainuu,Kajaani,64.224087,27.733423\nKanta-Häme,Hameenlinna,60.981091,24.458264\nKymenlaakso,Kotka,60.467423,26.945084\nLapland,Rovaniemi,66.497621,25.71921\nNorth Karelia,Joensuu,62.600816,29.760536\nNorth Ostrobothnia,Oulu,65.011873,25.471681\nOstrobothnia,Vaasa,63.081821,21.479812\nPirkanmaa,Tampere,61.498021,23.760312\nPohjois-Savo,Kuopio,62.824142,27.594561\nPäijät-Häme,Lahti,60.983876,25.656181\nSatakunta,Pori,61.486724,21.791002\nSouth Karelia,Lappeenranta,61.058242,28.18753\nSouth Ostrobothnia,Seinajoki,62.79541,22.844202\nUusimaa,Helsinki,60.167409,24.942568\nVarsinais-Suomi,Turku,60.451753,22.267052"

def get_loc_df():
    """
    Create region to lat/lon mapping
    """
    a = {}
    for entry in locations.split('\n'):
        region,_,*coord = entry.split(',')
        a[region] = list(map(float,coord))
    loc_df = pd.DataFrame.from_dict(a,orient='index',columns=['lat','lon'])
    loc_df.index.names = ['regionofbirth']
    return loc_df


def save_data(region_data,eigenvec_data,out_path):

    out_csv = os.path.join(out_path,"region_pc.csv")
    if not os.path.isfile(out_csv):
        # read in location to lat/lon df
        loc_df = get_loc_df()

        # read in birth data of samples
        birth_df = pd.read_csv(region_data,sep='\t',usecols = ["FINNGENID","regionofbirthname"],index_col='FINNGENID')
        usecols = ['PC' + str(pc) for pc in ["1","2","3"]]
        eigenvec = pd.read_csv(eigenvec_data, sep = '\t',index_col = 'IID')[usecols]
        eigenvec.index.names = ['FINNGENID']
        print(birth_df.head())
        print(eigenvec.head())

        #count people in each region
        count_df =  pd.concat([birth_df,eigenvec],join = 'inner',axis = 1)['regionofbirthname'].value_counts().to_frame('count')
        count_df.index.names = ['regionofbirth']
        print(count_df.head())

        # get averages of each region
        tmp_avg = pd.concat([birth_df,eigenvec],join = 'inner',axis = 1).set_index('regionofbirthname').groupby('regionofbirthname').mean()
        print(tmp_avg.head())
        #merge averages with location
        region_avg = pd.concat([loc_df,tmp_avg,count_df],join = 'inner',axis = 1)
        region_avg.index.names = ['regionofbirth']
        print(region_avg)
        region_avg.to_csv(out_csv,index=True)
        print(region_avg.shape)
        
    df = pd.read_csv(out_csv,index_col = 0)
    return df,out_path

def get_axis_limits(ax, scale=.9):
    return ax.get_xlim()[1]*(1-scale), ax.get_ylim()[1]*scale



def plot_pca_map(df,args):

    save_fig = os.path.join(args.plot_path, args.name + f"_pc_map.pdf")
    print(save_fig)

    if os.path.isfile(save_fig):
        return

    # set ip figure
    pc_list = [1,2,3]
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(1,len(pc_list))
    #set up min/max lon/lat
    lat1,lon1,lat2,lon2 = min(df.lat)-1,min(df.lon)-1,max(df.lat)+4,max(df.lon)+2
    print(lat1,lon1,lat2,lon2)

    # read in map
    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    fname = os.path.join(__location__, 'map','gadm36_FIN_2.shp')
    shape_feature = ShapelyFeature(Reader(fname).geometries(),crs.PlateCarree())
    stamen_terrain = cimgt.Stamen('terrain-background')

    # get extreme values of data for plotting
    pc_columns = ["PC" + str(elem) for elem in pc_list]
    values = df[pc_columns]
    limits = values.min().min(),values.max().max()
    minmax = abs(max(limits,key=abs))
    vmin,vmax = -minmax,minmax

    # resize dots based on population size
    resize = 200./max(df['count'])
    print(df['count'].values*resize)
    
    # turn count into
    axes = []
    for pc in pc_list:
        print(pc)
        ax = fig.add_subplot(gs[0,pc-1],projection=stamen_terrain.crs)
        axes.append(ax)
        ax.set_extent([lon1,lon2,lat1,lat2], crs=crs.PlateCarree())

        transform = crs.PlateCarree()._as_mpl_transform(ax)
        ax.annotate(str(pc), xy=(lon2- 0.05*(lon2-lon1),lat1+.08*(lat2-lat1)),xycoords=transform,ha='right', va='top')
        ax.add_feature(shape_feature,linewidth = 0.2, facecolor = (1, 1, 1, 0),edgecolor = "black")
        ax.add_feature(cartopy.feature.COASTLINE,linewidth=0.1,alpha =0.1)
        ax.add_feature(cartopy.feature.BORDERS, linestyle=':',linewidth=0.3,alpha =0.5)
        ax.add_feature(cartopy.feature.LAND)
        ax.add_feature(cartopy.feature.LAKES, alpha=0.3)
        ax.add_feature(cartopy.feature.RIVERS,alpha=0.3)

        im = plt.scatter(x=df.lon, y=df.lat,s=df['count'].values*resize, alpha=1,transform=crs.PlateCarree(),c = df['PC' + str(pc)].values,cmap='seismic',vmin=vmin,vmax=vmax)
        plt.scatter(x=df.lon, y=df.lat,s=df['count'].values*resize, alpha=1,facecolors='none', edgecolors='k',transform=crs.PlateCarree())

    magnitude = np.floor(np.log10(max(df['count'].values)))
    popList = [int(np.power(10,x)) for x in [magnitude -2, magnitude -1,magnitude]]
    print(popList)
    lines = [plt.scatter([], [], c='black', alpha=0.8, s= legend_pop*resize, label=str(legend_pop) + ' samples') for legend_pop in popList]
    
    legend1 = plt.legend(lines,[str(a)+ ' samples' for a in popList],loc = 'upper right',fontsize = 6)

    p0 = axes[0].get_position().get_points().flatten()
    p1 = axes[1].get_position().get_points().flatten()
    p2 = axes[2].get_position().get_points().flatten()
    ax_cbar = fig.add_axes([p0[0], .05, p2[2]-p0[0], 0.051])
    cbar=plt.colorbar(im, cax=ax_cbar, orientation='horizontal',aspect=20)
    cbar.set_ticks([vmin,0,vmax])
    for t in cbar.ax.get_xticklabels():t.set_fontsize(4)

    print('saving...')
    fig.savefig(save_fig)


def pca_map(args):
    loc_df = get_loc_df()
    print(loc_df)
    df,out_path = save_data(args.meta,args.eigenvec,os.path.join(args.plot_path,'plot_data'))
    print(df)
    plot_pca_map(df,args)


    

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
        eur_outliers = np.loadtxt(args.finngen_eur_outliers,dtype = str,usecols = 0)
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

    #if mapcount(args.finngen_eur_outliers) == 0:
        #return
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
        #args.v_print(3,'eur outliers 3d plot already done.')
        pass
    if not os.path.isfile(outliers_2d):
        plot_2d(pc_data,outliers_2d,tags,alpha_map= alpha_map,size_map = size_map,color_map = color_map,tag_column="TAG", legend_location = "lower right",legend_fontsize=10)
    else:
        ##args.v_print(3,'eur outliers pairwise plot already done.')
        pass


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
        # read pc data
        pc_data = pd.read_csv(eigenvec_path,sep = '\t',usecols = ['IID',"PC1",'PC2','PC3'], dtype = {pc: np.float64 for pc in ["PC1",'PC2','PC3']})
        # set finngen samples as "FINNGEN"
        outlier_data = pd.read_csv(outlier_info,dtype = str,sep = '\t',usecols = ['IID',"outlier",'SuperPops']).rename(columns ={"SuperPops":"TAG"})
        
        outlier_data.loc[((outlier_data['IID'].isin(samples)) & (outlier_data["outlier"] == "TRUE")),"TAG"] = "FINNGEN_OUTLIER"
        outlier_data.loc[((outlier_data['IID'].isin(samples)) & (outlier_data["outlier"] == "FALSE")),"TAG"] = "FINNGEN_INLIER"
        print(outlier_data.head())
        pc_data = pc_data.merge(outlier_data,on = "IID").set_index("IID")
        pc_data.to_csv(out_file)
        print(pc_data.head())


    pc_data = pd.read_csv(out_file,index_col =0)
    print(pc_data)
    #UPDATE TAGS SO I CAN REPLOT WITHOUT RECALCULATING OUTLIERS
    with open(args.annot_pop) as i: tag_dict = {sample:pop for sample,pop in (line.strip().split() for line in i)}    
    for index,sample in pc_data.iterrows():
        if sample.TAG not in ["FINNGEN_OUTLIER","FINNGEN_INLIER"]:
            pc_data.at[index,"TAG"] = tag_dict[index]

    tags = set(pc_data.TAG)
    print(tags)
    tag_size = [len(pc_data[pc_data.TAG == tag]) for tag in tags]
    # plot them by size!
    final_tags = [tag for _,tag in sorted(zip(tag_size,tags),reverse = True)]
    color_maps = list(color_dict[len(tags)]['qualitative'].keys())
    cm = 'Set1' if 'Set1' in color_maps else color_maps[0] 
    colors= color_dict[len(tags)]['qualitative'][cm]

    #[red,blue,green,purple,orange,yellow,brown,pink] = colors   
    #color_map = {'FIN':purple,'FINNGEN_INLIER':red,'FINNGEN_OUTLIER':green,'EUR':blue,'AFR':pink,'EAS':yellow,'SAS':brown,'AMR':orange}
    
    color_map ={final_tags[i]:color for i,color in enumerate(colors)}

    print(color_map)
    return final_tags,pc_data,color_map


def plot_first_round_outliers(args):
    '''
    Plots outlier detection results.
    '''

    ethnic_plot = os.path.join(args.plot_path,args.name +'_ethnic_outliers.pdf')
    ethnic_2d = os.path.join(args.plot_path,args.name +'_ethnic_outliers_pairwise.pdf')


    size_map = {"STRANGE_FIN":10}
    alpha_map = {"STRANGE_FIN":1}
    # createdata
    if not os.path.isfile(ethnic_plot) or not os.path.isfile(ethnic_2d):
        tags,pc_data,color_map = return_outliers_df(args)

    
    #plot 3d
    if not os.path.isfile(ethnic_plot):
        #build super pop dict for plotting
        plot_3d(pc_data,ethnic_plot,tags,tag_column="TAG",color_map = color_map)

    else:
        pass

    if not os.path.isfile(ethnic_2d):
        plot_2d(pc_data,ethnic_2d,tags,tag_column="TAG",color_map=color_map,size_map=size_map,alpha_map=alpha_map)
    else:
        pass



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
    fig.savefig(out_file.replace('.pdf','.png'),dpi=300)
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
    fig.savefig(out_file.replace('.pdf','.png'),dpi=300)
    plt.close()



def trim_axis(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.get_xaxis().tick_bottom()
    try:
        ax.get_yaxis().tick_left()
    except:
        pass



def plot_2d_density(pc_data,out_file,tags,pcs,color_map=None,tag_column="TAG",max_size = 3000,max_map=None,pc_tags=None,axis_legend =2,legend_location = "lower left",legend_fontsize=6,linewidths=None,levels = 3):


    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(3,1)

    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0],sharex = ax1)
    ax3 = fig.add_subplot(gs[2,0],sharex = ax1)
    axes=[ax1,ax2,ax3]

    if not pc_tags: pc_tags = pc_columns
    line_dict = dd(lambda:.6)
    if linewidths:
        for tag in linewidths :line_dict[tag] = linewidths[tag]

    if not color_map:
        color_maps = list(color_dict[len(tags)]['qualitative'].keys())
        cm = 'Set1' if 'Set1' in color_maps else color_maps[0]
        color_map = {tag:color_dict[len(tags)]['qualitative'][cm][i] for i,tag in enumerate(tags)}


    max_sizes = dd(lambda:max_size)
    if max_map:
        for key in max_map: max_sizes[key] = max_map[key]

    import matplotlib.patches as mpatches

    handles = []
    for i,tag in enumerate(tags):
        tag_data = pc_data[pc_data[tag_column] == tag]
        color = color_map[tag]
        max_data = max_sizes[tag]
        n = min(max_data,len(tag_data))
        tag_data = tag_data.sample(n=n)
        patch = mpatches.Patch(color=color, label=tag)
        handles.append(patch)
        print(tag,len(tag_data))
        for i,pcs in enumerate(list(combinations(pc_columns,2)) ):
            ax = axes[i]
            print(tag,pcs)

            x = tag_data[pcs[0]]
            y = tag_data[pcs[1]]
            xmin,xmax=min(x),max(x)
            ymin,ymax=min(y),max(y)
            # Peform the kernel density estimate
            xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]

            # DUMP DATA
            positions = np.vstack([xx.ravel(), yy.ravel()])
            values = np.vstack([x, y])
            kernel = st.gaussian_kde(values)
            f = np.reshape(kernel(positions).T, xx.shape)
            
            #data_path =os.path.splitext(out_file)[0]
            ##print(f"generating data --> {save_dump}")
            #np.savetxt(save_dump,f)
            #print(f"loading {save_dump}")
            #f = np.loadtxt(save_dump)

            ax.contour(xx, yy, f,colors=[color],linewidths=line_dict[tag],levels =levels,linestyles = 'dashed' )



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
    leg = leg_ax.legend(loc=legend_location,handles=handles,numpoints=1, fancybox = True,prop={'size': legend_fontsize})
    

    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    fig.savefig(out_file)
    fig.savefig(out_file.replace('.pdf','.png'),dpi=300)
    plt.close()
    

def pc_marginal(pc_data,bin_func,out_file,tags,pcs,color_map =None,alpha_map = None,size_map = None,max_map=None,max_size = 3000,tag_column = "TAG"):

    fig = plt.figure()

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

        
    gs = gridspec.GridSpec(2, 2, width_ratios=[1,3], height_ratios=[3,1])
    ax = plt.subplot(gs[0,1])
    axl = plt.subplot(gs[0,0], sharey=ax)
    axb = plt.subplot(gs[1,1], sharex=ax)
    for i,tag in enumerate(tags):
        tag_data = pc_data[pc_data[tag_column] == tag]
        max_data = max_sizes[tag]
        n = min(max_data,len(tag_data))
        tag_data = tag_data.sample(n=n)
        color = color_map[tag]
        size = sizes[tag]
        alpha = alphas[tag]
        print(tag,len(tag_data),alpha,size,n)
        plot_x,plot_y = tag_data[pcs[0]],tag_data[pcs[1]]
        ax.scatter(plot_x,plot_y, s= size,alpha = alpha,color = color,label = tag)
        bins = int(len(tag_data)/20)
        bin_data,plot_data = bin_func(tag_data[pcs[1]],bins)
        axl.plot(plot_data,bin_data,color=color)
        bin_data,plot_data = bin_func(tag_data[pcs[0]],bins)
        axb.plot(bin_data,plot_data,color=color)

    ax.tick_params(labelbottom=False)
    ax.tick_params(labelleft=False)

    ax.set_xlabel(pcs[1])
    ax.set_ylabel(pcs[0])
    
    print(out_file)
    fig.savefig(out_file)
    fig.savefig(out_file.replace('.pdf','.png'),dpi=300)
    
    return
    
def plot_2d_marginal(pc_data,bin_func,out_file,tags,pc_columns = ['PC1','PC2','PC3'],pc_tags = None,color_map= None,alpha_map = None,size_map = None,legend_fontsize = 6,label_fontsize = 5,tag_column = "TAG",max_size = 3000,max_map = None,axis_legend =2,legend_location = "lower left",rescale = 4):
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

    
    for i,pcs in enumerate(list(combinations(pc_columns,2)) ):
        pc_out = out_file.replace('.pdf',f'.{"_".join(pcs)}.pdf')
        pc_marginal(pc_data,bin_func,pc_out,tags,pcs)
        
