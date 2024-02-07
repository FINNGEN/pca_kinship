import argparse, os,subprocess,shlex,pickle,pathlib
import pandas as pd
import numpy as np
from collections import Counter
from collections import defaultdict
from scipy.stats import chi2
from scipy.spatial.distance import cdist
from utils import file_exists,make_sure_path_exists,pretty_print,basename,mapcount,progressBar
from pca_scripts.plot import plot_2d_density,plot_2d,trim_axis,plot_2d_marginal
from pca_scripts.color_dict import color_dict
from verkko.binner import binner

import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(palette='Set2')

def read_in_tags(fam_file,tags):

    # read in data from sample meta
    tag_dict = defaultdict(lambda : "Other")
    with open(tags) as i:
        for line in i:
            sample,tag= line.strip().split('\t')
            tag_dict[sample] = tag

    # assign to each sample in bed file a tag
    fam_tag_dict = {}
    with open(fam_file) as i:
        for line in i:
            sample = line.strip().split()[1]
            fam_tag_dict[sample] = tag_dict[sample]
    print(Counter(fam_tag_dict.values()))
    return fam_tag_dict



def merge_pca(pca_root,ref_bed,proj_bed,plink_cmd,extract=None,force=False):
    """
    Here we merge two datasets and run PCA the merged set.
    """
    # i defined the merged file only using input files names, so i don't regenerate them everytime.
    plink_root = os.path.join(os.path.dirname(os.path.dirname(pca_root)),'plink')
    make_sure_path_exists(plink_root)
    merged_plink = os.path.join(plink_root,pathlib.Path(ref_bed).stem + "_" + pathlib.Path(proj_bed).stem + "_merged.bed")
    if not os.path.isfile(merged_plink):
        print ('Merged dataset missing')
        cmd = f"plink --bfile {basename(ref_bed)} --bmerge {basename(proj_bed)} --make-bed --out {basename(merged_plink)}"
        print(cmd)
        subprocess.call(shlex.split(cmd))
    else :
        print(f"{merged_plink} bed file already generated")

    freq = basename(merged_plink) + ".afreq"
    if not os.path.isfile(freq):
        cmd = f"plink2 --bfile {basename(merged_plink)} --freq --out {basename(merged_plink)}"
        subprocess.call(shlex.split(cmd))
    # PCA
    eigenvec = pca_root + '.eigenvec'
    if not os.path.isfile(eigenvec) or force:
        approx = "approx" if mapcount(basename(merged_plink) +'.fam') > 5000 else ""
        extract = f" --extract {extract} " if extract else "" 
        cmd = f"{plink_cmd} --bfile {basename(merged_plink)} {extract} --write-snplist   --pca 3 {approx} biallelic-var-wts -out {pca_root} --read-freq {freq}"
        print(cmd)
        subprocess.call(shlex.split(cmd))
    else:
        print('PCA already calculated')

    # SPLIT AND CHANGE HEADER TO MATCH PROJ DATA
    ref_score = pca_root +  '_ref.sscore'
    proj_score =pca_root +  "_proj.sscore"
    proj_iids = set(np.loadtxt(basename(proj_bed) + '.fam',dtype = str,usecols = 1))
    ref_iids = set(np.loadtxt(basename(ref_bed) + '.fam',dtype = str,usecols = 1))

    with open(ref_score,'wt') as ref,open(proj_score,'wt') as proj,open(eigenvec) as i:
        new_header = '\t'.join([elem +"_AVG"  if elem.startswith("PC") else elem for elem in next(i).strip().split()]) + '\n'
        proj.write(new_header)
        ref.write(new_header)
        for line in i:
            line = line.strip().split()
            iid = line[1]
            if line [0] == "0": line[0] = iid
            line = '\t'.join(line) + '\n'            
            if iid in ref_iids:ref.write(line)
            elif iid in proj_iids:proj.write(line)
            else:"iid missing, we have a problem!"
                
    return ref_score,proj_score


def run_pca(pca_root,ref_bed,proj_bed,plink_cmd,flip,extract = None):
    """
    Here we calculate the PC on the ref set and project the other group onto it.
    """
    eigenvec = pca_root + '.eigenvec.var'
    print(eigenvec)
    if not os.path.isfile(eigenvec):
        extract = f" --extract {extract} " if extract else ""
        pca_bed = proj_bed if flip else ref_bed
        approx = "approx" if mapcount(pca_bed.replace('.bed','.fam')) > 5000 else ""
        freq = f"--read-freq {pca_bed.replace('.bed','.afreq')}" if {os.path.isfile(pca_bed.replace('.bed','.afreq'))} else " "
        cmd = f"{plink_cmd} --bfile {basename(pca_bed)} {extract}  --pca 3  biallelic-var-wts -out {pca_root} {freq}"
        
        subprocess.call(shlex.split(cmd))
    else:
        print('PCA already calculated')
  
    proj_cmd = f"plink2 --score {eigenvec}  2 3 header-read  no-mean-imputation   --score-col-nums 5-7  "
    #FREQ FILES
    for bed_file in [ref_bed,proj_bed]:
        bed_root=basename(bed_file)
        if not  os.path.isfile(bed_root + '.afreq'):
            subprocess.call(shlex.split(f"plink2 --bfile {bed_root} --freq --out {bed_root}"))
    #REF BED
    ref_score = pca_root +  '_ref.sscore'
    if not os.path.isfile(ref_score):
        freq_cmd = f"plink2 --bfile {basename(ref_bed)}"    
        ref_cmd = proj_cmd + f" --bfile {basename(ref_bed)} --out {basename(ref_score)} --read-freq {basename(ref_bed)}.afreq"
        subprocess.call(shlex.split(ref_cmd))
    else :
        print(f"{ref_score} already projected.")
              
    proj_score =pca_root +  "_proj.sscore"
    if not os.path.isfile(proj_score):
        test_cmd = proj_cmd + f" --bfile {basename(proj_bed)} --out {basename(proj_score)} --read-freq {basename(proj_bed)}.afreq"
        subprocess.call(shlex.split(test_cmd))
    else :
        print(f"{proj_score} already projected.")

    return ref_score,proj_score


def plot_projection(ref_scores,proj_scores,plot_root,tag_dict,top_regions):

    # read in data
    plot_data = plot_root + '_proj.csv'
    if not os.path.isfile(plot_data):
        pc_avg = ["PC1_AVG",'PC2_AVG','PC3_AVG']
        ref_df =  pd.read_csv(ref_scores,index_col = 0,sep = '\t',usecols = ['IID'] + pc_avg, dtype = {pc: np.float64 for pc in pc_avg}).rename(columns = {pc: pc.replace("_AVG","") for pc in pc_avg})
        ref_df['TAG'] = "core"
        test_df =  pd.read_csv(proj_scores,sep = '\t',index_col = 0,usecols = ['IID'] + pc_avg, dtype = {pc: np.float64 for pc in pc_avg}).rename(columns = {pc: pc.replace("_AVG","") for pc in pc_avg})
        test_df['TAG'] = "proj"
        df = pd.concat([ref_df,test_df])
        df.to_csv(plot_data)

    df = pd.read_csv(plot_data,index_col=0)
    tags = list(set(df.TAG))
    scatter_fig = plot_root + '_scatter_all.pdf'
    density_fig = plot_root + '_scatter_all_density.pdf'
    print('reading in data')
    color_map = {"proj":(1,0,0),'core':(0,0,1)}
    print(tags)
    plot_2d(df,scatter_fig,tags=tags,color_map=color_map,max_size = 20000,alpha_map={"core":.1,'proj':.3})
    #plot_2d_density(df,density_fig,tags=tags,color_map=color_map,max_size=20000)



    #### TAG/POP STUFF ####
    tag_df = pd.DataFrame(tag_dict.items(),columns=["IID","TAG"]).set_index("IID")
    df = pd.read_csv(plot_data,index_col=0)
    df.update(tag_df)
    tags = list(set(df.TAG))
    count = dict(df.TAG.value_counts())
    print(count)
    
    # use all tags if low enough else keep only relevant ones
    tags = [tag for tag in tags if count[tag] > 10]
    tags = tags if len(tags) <= 10 else ["proj"] + top_regions
    
    tag_scatter = plot_root + '_scatter_tags.pdf'
    tag_density = plot_root + '_scatter_tags-density.pdf'
    
    plot_2d(df,tag_scatter,tags=tags,max_size = 20000)
    lw= {elem:.3 for elem in tags}
    lw["proj"] = 1

    plot_2d_marginal(df,return_bin_data,scatter_fig,tags=tags,color_map=color_map,max_size = 20000,alpha_map={"core":.1,'proj':.3})
    plot_tags(df,plot_root,tags)
    #plot_2d_density(df,tag_density,tags=tags,max_size=np.inf,linewidths=lw,levels = 2)
    
    

def return_bin_data(data,n_bins =40):
    xmin,xmax = data.to_numpy().min(),data.to_numpy().max()
    bins = binner.Bins(float,xmin,xmax,'lin',n_bins)
    countNotNormalized = bins.bin_count_divide(data)
    count = np.array(binner.normalize(list(countNotNormalized)))
    binAvg = bins.bin_average(zip(data,data))
    binMask = ~np.ma.getmask(binAvg)
    plot_data = count[binMask]
    bin_data = binAvg[binMask]
    return bin_data,plot_data


def plot_tags(df,plot_root,tags):

    save_fig = plot_root + '_tags_pc_density.pdf'

    tags.insert(0, tags.pop(tags.index("proj")))

    print(len(tags))
    print(tags)
    scheme = "Set1" if len(tags) <10 else "Set3"
    colors = color_dict[len(tags)]['qualitative'][scheme]
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(3,1)
    pc_tags =["PC1",'PC2','PC3']
    for i,pc in enumerate(pc_tags):
        print(pc)
        ax = fig.add_subplot(gs[i,0])
        ax.set_ylabel(r'P('+pc+')')
        pc_data = df[pc]              
        x_min,x_max,y_max = [],[],[]
        for j,tag in enumerate(tags):
            tag_data = pc_data[df.TAG==tag]
            bin_data,plot_data = return_bin_data(tag_data)           
            lw = 1 if tag =='proj' else .7
            ls = '-' if tag =='proj' else '--'
            ax.plot(bin_data,plot_data,ls,color = colors[j],label=tag,linewidth=lw)
            y_max.append(max(plot_data))
            x_min.append(bin_data[0])
            x_max.append(bin_data[-1])
        bin_data,plot_data = return_bin_data(pc_data,100)
        ax.plot(bin_data,plot_data,1,color = 'k',label="core",linewidth=1)
        ax.set_xlim(min(x_min),max(x_max))
        ax.set_ylim(0,max(max(y_max),max(plot_data)))
        
        trim_axis(ax)
        for tick in ax.xaxis.get_major_ticks():tick.label1.set_fontsize(6)
        for tick in ax.yaxis.get_major_ticks():tick.label1.set_fontsize(6)
        
    handles,labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    leg = ax.legend(by_label.values(),by_label.keys(),loc="lower left", numpoints=1, fancybox = True,prop={'size':4})
    
    fig.savefig(save_fig)
    fig.savefig(save_fig.replace('.pdf','.png'))
    plt.close()

    return
                

def generate_tag_data(tag_dict,ref_scores):
    """
    Here i generate the probs for all tags.
    """

    summary = {}
    pc_data = pd.read_csv(ref_scores,usecols=['IID','PC1_AVG','PC2_AVG','PC3_AVG'],index_col = 0,sep='\t')

    tags = set(tag_dict.values())
    print(tags)
    count = Counter(tag_dict.values())
    print(count)
    tags = [tag for tag in count if count[tag] > 10]
    print(tags)

    for tag in tags:
        samples = [elem for elem in tag_dict if tag_dict[elem] == tag]
        tag_data = pc_data.loc[samples].to_numpy()
        tag_avg = np.reshape(np.average(tag_data,axis=0),(1,tag_data.shape[1]))
        tag_cov = np.linalg.inv(np.cov(tag_data.T))
        summary[tag]=[tag_avg,tag_cov]

    return summary


def generate_df(proj_scores,ref_scores,tag_dict,out_root):

    # read in pc data
    proj_data = pd.read_csv(proj_scores,usecols=['IID','PC1_AVG','PC2_AVG','PC3_AVG'],index_col = 0,sep='\t')
    print(proj_data)
    df_prob = pd.DataFrame(index=proj_data.index)
    # samples
    samples = proj_data.index.values
    # convert df to array
    proj_data = proj_data.to_numpy()
        
    print(f"Calculating tag avg/cov...")
    # calculate average and covariance for each core group
    tag_avg_cov = generate_tag_data(tag_dict,ref_scores)
    for tag in tag_avg_cov:
        print(tag)
        avg,cov = tag_avg_cov[tag]
        tag_dist = cdist(proj_data,avg,metric = 'mahalanobis',VI = cov).flatten()**2
        tag_prob = 1 - chi2.cdf(tag_dist,3)
        df_prob[tag] = tag_prob

    if "NA" in df_prob: df_prob = df_prob.drop(["NA"],axis=1)
    print(df_prob)

    df_prob.to_csv(out_root + '_probs.txt')
    return df_prob,out_root

def calculate_probs(df_prob,cutoffs,out_root):
    # top region
    hit_regions = df_prob.idxmax(axis=1) # get column of max value
    top_regions = [elem for elem  in Counter(hit_regions).most_common()]

    region_root = out_root + "_samples_most_likely_region.txt"
    with open(region_root,'wt') as o:
        for entry in hit_regions.to_dict().items():
            o.write('\t'.join(list(map(str,entry))) + '\n')

    # now i filter the df for all cutoff values
    for cutoff in sorted(cutoffs):
        tmp_df = df_prob.copy(deep=True)
        with open(region_root.replace(".txt",f"_{cutoff}.txt"),'wt') as o:
            tmp_df[df_prob<cutoff] = 0 #set values under threshold to 0
            # zip together index,max val and region of max val
            for entry in zip(df_prob.index,tmp_df.max(axis=1),tmp_df.idxmax(axis=1)): 
                sample,value,region = entry
                # if value is non 0 it means it's above threshold
                region = region if float(value) else "NA"
                o.write('\t'.join([sample,region]) + '\n')   


    if "FIN" in list(df_prob.columns.values):
        df = df_prob.drop(columns=["FIN"])
        calculate_probs(df,cutoffs,out_root + "_finless")
                
    return top_regions


    

def main(args):
    pretty_print("TAG DICT")
    tag_dict = read_in_tags(args.ref_bed.replace('.bed','.fam'),args.sample_info)

    pretty_print("PCA")
    pca_path = os.path.join(args.out_path,'pca')
    make_sure_path_exists(pca_path)
    pca_root = os.path.join(pca_path,args.name)

    plink_cmd = f"plink2 {args.pargs}"
    if args.merge:
        pretty_print("MERGE-PCA")
        ref_scores,proj_scores =  merge_pca(pca_root,args.ref_bed,args.proj_bed,plink_cmd,args.extract,args.force)
    else:
        ref_scores,proj_scores = run_pca(pca_root,args.ref_bed,args.proj_bed,plink_cmd,args.flip,args.extract)

    plot_path = os.path.join(args.out_path,'plot')
    plot_root = os.path.join(plot_path,args.name)

    
    df,out_root = generate_df(proj_scores,ref_scores,tag_dict,os.path.join(args.out_path,args.name))
    top_regions = calculate_probs(df,args.cutoffs,out_root)
    if args.plot:
        pretty_print("PLOT PROJECTION")
        make_sure_path_exists(plot_path)
        df = plot_projection(ref_scores,proj_scores,plot_root,tag_dict,list(list(zip(*top_regions))[0]))

    pretty_print("TAG PROJECTIONS")
    print(top_regions)
    return



if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description="Classify samples based on PC.")

    # BED FILES
    parser.add_argument("--ref-bed", type=file_exists, help = "Bed file for reference PC", required = True)
    parser.add_argument('--proj-bed',type = file_exists,help = 'Bed file for samples that need to be projected.')

    # GENERAL PARAMS
    parser.add_argument("--name", type=str,help ="prefix of output files")
    parser.add_argument('-o',"--out_path",type = str, help = "Folder in which to save the results", required = True)
    parser.add_argument("--sample-info", type=file_exists, help =  "Tsv file with sample data, used for grouping.", required = True)
    parser.add_argument('--plot',action = 'store_true',help = 'Plotting',default = False)
    parser.add_argument('--force',action = 'store_true',help = 'Recompute PCA',default = False)
    parser.add_argument('--flip',action = 'store_true',help = 'Inverts choice of ref/proj for calculating PCs',default = False)
    parser.add_argument( '--cutoffs', nargs='+', type = str, default=[0.5,0.8,0.9,0.99])

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--proj', action='store_true')
    group.add_argument('--merge', action='store_true')

    parser.add_argument("--pargs", type=str,help ="Plink pca args",default = " ")
    parser.add_argument("--extract", type=file_exists, help =  "Snps to use.", required = False)
    args = parser.parse_args()

    make_sure_path_exists(args.out_path)
    name_tag = "_merged" if args.merge else "_proj"
    flip_tag = "_flip" if args.flip else ""
    args.name += name_tag + flip_tag

    main(args)
