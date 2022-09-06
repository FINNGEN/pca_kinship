from utils import np,mapcount,basic_iterator,return_header,write_fam_samplelist,make_sure_path_exists
import pandas as pd
import csv
from scipy.spatial.distance import cdist
from scipy.stats import chi2
import os,shlex,subprocess



def detect_ethnic_outliers(args):
    '''
    Returns list of FinngGen samples that are too far from the core PCs.
    '''
    tg_pca_path=  os.path.join(args.pca_outlier_path, '1k_pca/')
    make_sure_path_exists(tg_pca_path)
    tg_pca_file = pca_round(args,tg_pca_path)

    return None
def pca_round(args,tg_pca_path,remove_list=None):
    '''
    Function that runs PCA on the merged dataset, with option to remove samples for iterative runs.
    '''

    tg_pca_file = os.path.join(tg_pca_path,args.name)

    if mapcount(tg_pca_file+ '.eigenvec') != mapcount(args.merged_plink_file +'.fam') +1 or args.force:
        args.force = True 
        #individuals that need to be removed
        remove = ''
        if remove_list is not None:
            remove_file = args.misc_path +  'remove.fam'
            write_fam_samplelist(remove_file,remove_list)
            remove =  f' --remove {remove_file}'             
            args.logging.debug(remove)
            
        cmd = f'plink2 --bfile {args.merged_plink_file} --read-freq  {args.merged_plink_file}.afreq {remove} --pca {args.pca_components} approx biallelic-var-wts --threads {args.cpus}  -out {tg_pca_file}'
        args.logging.debug(cmd)
        subprocess.call(shlex.split(cmd))
        
    return tg_pca_file

def build_superpop(args):
    '''
    Writes pop info for outlier detection.
    Returns filepath of file needed for outlier detection.
    '''
    
    # sample to super pop mapping file
    annot_pop =  os.path.join(args.pca_outlier_path, args.name + '_sample_annot.tsv')
    if not os.path.os.path.isfile(annot_pop) or args.force or mapcount(annot_pop) < 2:
        args.force = True 
        args.logging.info('generating super pop dict {annot_pop}')

        # pop to superpop mapping
        with open(args.data_path + 'superpop.csv') as i:
            superpop_dict=  {rows[0]:rows[1] for rows in csv.reader(i)}
        print(superpop_dict)

        # sample info for tg
        tg_pop = os.path.join(args.data_path,'20130606_sample_info.txt')
           # get index of pop column and build sample to population dictionary
        cols = [return_header(tg_pop).index(elem) for elem in ['Sample','Population']]
        pop_dict = {sample:pop for sample,pop in basic_iterator(tg_pop,columns = cols)}

        # now i build a sample to pop dictionary where i keep superpop unless it's a Finn
        with open(annot_pop,'wt') as o:
            o.write('IID' +'\t' + 'SuperPops\n')

            # loop through input tg fam file and update population data
            for sample in basic_iterator(args.tg_bed.replace('.bed','.fam'),columns = 1):
                pop = pop_dict[sample] # fetch country
                super_pop = pop if pop =='FIN'else superpop_dict[pop]
                o.write('\t'.join([sample,super_pop]) +'\n')

            #add batch info from FINNGEN samples
            cols = [return_header(args.sample_info).index(elem) for elem in ['IID','BATCH']]
            for sample,batch in basic_iterator(args.sample_info,columns =cols,separator =','):
                o.write('\t'.join([sample,batch]) +'\n')
                  

    return annot_pop
