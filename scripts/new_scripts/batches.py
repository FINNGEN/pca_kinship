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
    all_cohorts = os.path.join(args.misc_path,'cohorts.txt')
    all_batches = os.path.join(args.misc_path,'batches.txt')

    if not all(list(map(os.path.isfile,[sample_info,all_cohorts,all_batches]))) or args.force:
        args.logging.debug(f"Saving sample metadata to {sample_info}")
        #save metadata to custom format
        sample_data = pd.read_csv(args.sample_data,sep = identify_separator(args.sample_data),usecols = ['BATCH','RELEASE','COHORT','FINNGENID']).rename(columns={"FINNGENID": "IID" })
        # for axiom batches, replace release as cohort
        axiom_mask =(sample_data['COHORT'] == 'Other') & sample_data['BATCH'].str.contains('Axiom')
        sample_data.loc[axiom_mask,'COHORT'] = sample_data.loc[axiom_mask,'RELEASE']             
        sample_data.to_csv(sample_info,index = False,columns = ('BATCH','COHORT','IID'))
        
        # write batches summary
        args.logging.debug(f"Saving batch/cohort data to {all_cohorts} {all_batches}")

        batches = sample_data.BATCH.unique()
        np.savetxt(all_batches,batches,fmt="%s")
        cohorts = sample_data.COHORT.unique()
        np.savetxt(all_cohorts,cohorts,fmt = '%s' )

    print(f"{mapcount(all_cohorts)} cohorts ")
    print(f"{mapcount(all_batches)} batches")

        
    return sample_info,all_cohorts,all_batches
