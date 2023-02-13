import argparse, os,subprocess,shlex
import pandas as pd
import numpy as np
from collections import defaultdict
from utils import file_exists,make_sure_path_exists,pretty_print,basename
from pca_scripts.plot import plot_2d_density,plot_2d


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
            
    return fam_tag_dict


def run_pca(bed,pca_path):

    pca_root = os.path.join(pca_path,'ref_pca')
    eigenvec = pca_root + '.eigenvec.var'
    if not os.path.isfile(eigenvec):
        cmd = f"plink2 --bfile {basename(bed)} --pca 10 approx biallelic-var-wts -out {pca_root}"
        subprocess.call(shlex.split(cmd))
    else:
        print('PCA already calculated')
    return eigenvec


def project(eigenvec,pca_path,ref_bed,test_bed):
    '''
    Project all samples onto same space
    '''
    proj_cmd = f"plink2 --score {eigenvec}  2 3 header-read no-mean-imputation variance-standardize  --score-col-nums 5-14  "

    #FREQ FILES
    for bed_file in [ref_bed,test_bed]:
        bed_root=basename(bed_file)
        if not  os.path.isfile(bed_root + '.afreq'):
            subprocess.call(shlex.split(f"plink2 --bfile {bed_root} --freq --out {bed_root}"))
    #REF BED
    ref_score = os.path.join(pca_path, 'ref_proj.sscore')
    if not os.path.isfile(ref_score):
        freq_cmd = f"plink2 --bfile {basename(ref_bed)}"    
        ref_cmd = proj_cmd + f" --bfile {basename(ref_bed)} --out {basename(ref_score)} --read-freq {basename(ref_bed)}.afreq"
        subprocess.call(shlex.split(ref_cmd))
    else :
        print(f"{ref_score} already projected.")
              
    test_score = os.path.join(pca_path, "test_proj.sscore")
    if not os.path.isfile(test_score):
        test_cmd = proj_cmd + f" --bfile {basename(test_bed)} --out {basename(test_score)} --read-freq {basename(test_bed)}.afreq"
        subprocess.call(shlex.split(test_cmd))
    else :
        print(f"{test_score} already projected.")

    # MERGE SCORES
    all_scores = os.path.join(pca_path,'all.sscore')
    with open(all_scores,'wt') as o, open(ref_score) as ref,open(test_score) as test:
        for line in ref:o.write(line)
        next(test)
        for line in test:o.write(line)
        
    return all_scores,ref_score,test_score



def plot_projection(ref_scores,test_scores,plot_path):

    # read in data
    plot_data = os.path.join(plot_path,'proj.csv')
    if not os.path.isfile(plot_data):
        pc_avg = ["PC1_AVG",'PC2_AVG','PC3_AVG']
        ref_df =  pd.read_csv(ref_scores,index_col = 0,sep = '\t',usecols = ['IID'] + pc_avg, dtype = {pc: np.float64 for pc in pc_avg}).rename(columns = {pc: pc.replace("_AVG","") for pc in pc_avg})
        ref_df['TAG'] = "core"
        test_df =  pd.read_csv(test_scores,sep = '\t',index_col = 0,usecols = ['IID'] + pc_avg, dtype = {pc: np.float64 for pc in pc_avg}).rename(columns = {pc: pc.replace("_AVG","") for pc in pc_avg})
        test_df['TAG'] = "proj"
        df = pd.concat([ref_df,test_df])
        df.to_csv(plot_data)
    else:
        print('reading in data')
        df = pd.read_csv(plot_data,index_col=0)
    print(df)
    
    color_map = {"proj":'red','core':'blue'}
    tags = list(set(df.TAG))
    print(tags)
    
    fig_path = os.path.join(plot_path,'projection.pdf')
    print(fig_path)
    if not os.path.isfile(fig_path):
        alpha_map = {"core":.1,'proj':.3}
        plot_2d(df,fig_path,tags=set(df.TAG),color_map=color_map,max_size = 10000,alpha_map=alpha_map)
    fig_path = os.path.join(plot_path,'projection_density.pdf')
    print(fig_path)
    if not os.path.isfile(fig_path):
        plot_2d_density(df,fig_path,tags=set(df.TAG),data_path=plot_path,color_map=color_map,max_size=np.inf)

def main(args):
    pretty_print("TAG DICT")
    tag_dict = read_in_tags(args.ref_bed.replace('.bed','.fam'),args.sample_info)

    pretty_print("PCA")
    pca_path = os.path.join(args.out_path,'pca')
    make_sure_path_exists(pca_path)
    eigenvec= run_pca(args.ref_bed,pca_path)

    pretty_print("PROJ")
    all_scores,ref_scores,test_scores = project(eigenvec,pca_path,args.ref_bed,args.test_bed)

    pretty_print("PLOT PROJECTION")
    plot_path = os.path.join(args.out_path,'plot')
    make_sure_path_exists(plot_path)
    plot_projection(ref_scores,test_scores,plot_path)
    return



if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description="Classify samples based on PC.")

    # BED FILES
    parser.add_argument("--ref-bed", type=file_exists, help = "Bed file for reference PC", required = True)
    parser.add_argument('--test-bed',type = file_exists,help = 'Bed file for samples that need to be projected.')

    # GENERAL PARAMS
    parser.add_argument('-o',"--out_path",type = str, help = "Folder in which to save the results", required = True)
    parser.add_argument("--sample-info", type=file_exists, help =  "Tsv file with sample data, used for grouping.", required = True)


    args = parser.parse_args()
    make_sure_path_exists(args.out_path)
    main(args)
