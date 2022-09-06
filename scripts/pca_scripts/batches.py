import os,subprocess,shlex
import pandas as pd
import numpy as np
from collections import Counter
from utils import mapcount,write_fam_samplelist,identify_separator,tmp_bash,mem_mib

def batches(args):
    '''
    Creates batch metadata and, if needed, a new plink file with a subset of samples.
    '''

    #FILTER METADATA
    new_sample = os.path.join(args.misc_path,'sample_info.csv')
    if not os.path.isfile(new_sample) or mapcount(new_sample) < 1 or args.force:
        args.force = True
        #save metadata to custom format
        sample_info = pd.read_csv(args.sample_info,sep = identify_separator(args.sample_info),usecols = ['BATCH','RELEASE','COHORT','FINNGENID']).rename(columns={"FINNGENID": "IID" })
        # for axiom batches, replease release as cohort
        axiom_mask =(sample_info['COHORT'] == 'Other') & sample_info['BATCH'].str.contains('Axiom')
        sample_info.loc[axiom_mask,'COHORT'] = sample_info.loc[axiom_mask,'RELEASE']             

        sample_info.loc[:,('BATCH','COHORT','IID')].to_csv(new_sample,index = False)
        
    # Replace args.sample_info with the filtered version
    args.sample_info = new_sample
        
    # write batches summary
    args.cohorts = os.path.join(args.misc_path,'cohorts.txt')
    args.batches = os.path.join(args.misc_path,'batches.txt')
    if not os.path.isfile(args.cohorts) or args.force or not os.path.isfile(args.batches):
        args.force = True
        sample_info = pd.read_csv(args.sample_info,sep = identify_separator(args.sample_info))
        batches = np.array(list(set(sample_info[['BATCH']].values.flatten())),dtype = str)
        np.savetxt(args.batches,batches,fmt = '%s' )
        cohorts = np.array(list(set(sample_info[['COHORT']].values.flatten())),dtype = str)
        np.savetxt(args.cohorts,cohorts,fmt = '%s' )
        
    print(f"{mapcount(args.cohorts)} cohorts ")
    print(f"{mapcount(args.batches)} batches")
    
def merge_beds(args):
    args.sample_fam =os.path.join(args.misc_path, 'samples.fam')                                                                  
    if not os.path.isfile(args.sample_fam) or mapcount(args.sample_fam) < 1 or args.force:
        args.force = True
        if args.test and mapcount(args.bed.replace('.bed','.fam')) > args.test :
            # filter fam file to random n samples
            tmp_bash(f"cat {args.bed.replace('.bed','.fam')} | shuf | head -n {args.test} > {args.sample_fam}")           
            cmd = f"plink2 --bfile {args.bed.replace('.bed','')}  --freq --memory {int(mem_mib)} --make-bed --out {os.path.join(args.plink_path,args.name)} --keep {args.sample_fam}"
            subprocess.call(shlex.split(cmd))
            
        else:
            cmd = f"cat {args.bed.replace('.bed','.fam')}  > {args.sample_fam}"
            tmp_bash(cmd)
        
    # UPDATE BED PATH IF RE-RUN
    if args.test:
        args.bed = os.path.join(args.plink_path,args.name) + '.bed'

    args.n_samples = mapcount(args.sample_fam)
    print(f"{args.n_samples} samples used for PCA")

        
    
