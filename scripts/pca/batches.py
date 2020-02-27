import os,subprocess,shlex
import pandas as pd
import numpy as np
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
        axiom_mask =(sample_info['COHORT'] == 'Other') & sample_info['BATCH'].str.contains('Axiom')
        sample_info.loc[axiom_mask,'COHORT'] = sample_info.loc[axiom_mask,'RELEASE']
        sample_info.loc[:,('BATCH','COHORT','IID')].to_csv(new_sample,index = False)
        args.v_print(1,sample_info.head())      

    args.sample_info = new_sample
    sample_info = pd.read_csv(args.sample_info,sep = identify_separator(args.sample_info))
   
    # write batches summary
    args.cohorts = args.misc_path +'cohorts.txt'
    args.batches = args.misc_path +'batches.txt'
    if not os.path.isfile(args.cohorts) or args.force or not os.path.isfile(args.batches):
        args.force = True
        batches = np.array(list(set(sample_info[['BATCH']].values.flatten())),dtype = str)
        np.savetxt(args.batches,batches,fmt = '%s' )
        cohorts = np.array(list(set(sample_info[['COHORT']].values.flatten())),dtype = str)
        np.savetxt(args.cohorts,cohorts,fmt = '%s' )
    print(f"{mapcount(args.cohorts)} cohort ")
    print(f"{mapcount(args.batches)} batches")
    
    # CREATE NEW BED FILE IF NEEDED 
    args.sample_fam =os.path.join(args.misc_path, 'samples.fam')                                                                  
    if not os.path.isfile(args.sample_fam) or mapcount(args.sample_fam) < 1 or args.force:
        args.force = True
        if args.test and mapcount(args.bed.replace('.bed','.fam')) > args.test :
            # filter fam file to random n samples
            tmp_bash(f"cat {args.bed.replace('.bed','.fam')} | shuf | head -n {args.test} > {args.sample_fam}")
            new_plink_file = os.path.join(args.plink_path,args.name)
            cmd = f"plink --bfile {args.bed.replace('.bed','')}  --freq --memory {int(mem_mib)} --make-bed --out {new_plink_file} --keep {args.sample_fam}"
            subprocess.call(shlex.split(cmd))
            args.bed = new_plink_file + '.bed'
        else:
            cmd = f"cat {args.bed.replace('.bed','.fam')}  > {args.sample_fam}"
            tmp_bash(cmd)
            
    args.n_samples = mapcount(args.sample_fam)
    print(f"{args.n_samples} samples uesd for PCA")

        
    
