import os
import pandas as pd
import numpy as np
from collections import Counter
from utils import mapcount,identify_separator

def batches(args):
    '''
    Creates batch metadata and, if needed, a new plink file with a subset of samples.
    '''

    #FILTER METADATA
    sample_info = os.path.join(args.misc_path,'sample_info.csv')
    if not os.path.isfile(sample_info) or mapcount(sample_info) < 1 or args.force:
        args.force = True
        args.logging.info(f"Saving sample metadata to {sample_info}")
        #save metadata to custom format
        sample_data = pd.read_csv(args.sample_info,sep = identify_separator(args.sample_info),usecols = ['BATCH','RELEASE','COHORT','FINNGENID']).rename(columns={"FINNGENID": "IID" })
        # for axiom batches, replace release as cohort
        axiom_mask =(sample_data['COHORT'] == 'Other') & sample_data['BATCH'].str.contains('Axiom')
        sample_data.loc[axiom_mask,'COHORT'] = sample_data.loc[axiom_mask,'RELEASE']             
        sample_data.to_csv(sample_info,index = False,columns = ('BATCH','COHORT','IID'))
        

    # write batches summary
    all_cohorts = os.path.join(args.misc_path,'cohorts.txt')
    all_batches = os.path.join(args.misc_path,'batches.txt')
    if not os.path.isfile(all_cohorts) or args.force or not os.path.isfile(all_batches):
        args.force = True
        args.logging.info(f"Saving batch/cohort data to {all_cohorts} {all_batches}")

        sample_info = pd.read_csv(args.sample_info,sep = identify_separator(args.sample_info),dtype=str)
        batches = np.array(list(set(sample_info[['BATCH']].values.flatten())),dtype = str)
        np.savetxt(all_batches,batches,fmt = '%s' )
        cohorts = np.array(list(set(sample_info[['COHORT']].values.flatten())),dtype = str)
        np.savetxt(all_cohorts,cohorts,fmt = '%s' )
        
    print(f"{mapcount(all_cohorts)} cohorts ")
    print(f"{mapcount(all_batches)} batches")

        
    return sample_info,all_cohorts,all_batches
