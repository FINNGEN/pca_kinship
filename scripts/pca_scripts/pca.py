from utils import mapcount,basic_iterator,progressBar,merge_files,return_header,mem_mib,write_fam_samplelist
import numpy as np
import os.path,shlex,subprocess

def build_inliers(args):
    '''
    Returns .fam file for running pca by excluding samples that are marked as related or as false finns. Also, it returns the intersection of the false_finns and duplicates, so that they can be discarded.
    '''

    #LOAD ALL RELEVANT LIST OF SAMPLES
    related_samples = np.loadtxt(args.related_individuals,dtype = str)
    all_outliers = np.loadtxt(args.all_outliers,dtype = str,usecols =0)
    duplicates = np.loadtxt(args.duplicates,dtype = str)
    # INLIERS
    unrelated_file = os.path.join(args.pca_path,args.name +'_unrelated.txt')
    related_file = os.path.join(args.pca_path,args.name +'_related.txt')
    rejected_file = os.path.join(args.pca_path,args.name +'_rejected.txt')
    final_samples = os.path.join(args.pca_path,args.name +'_final_samples.txt')

    n_samples = mapcount(args.sample_fam)
    if not all(map(os.path.isfile,[unrelated_file,related_file,rejected_file,final_samples])) or args.force:
        args.force = True 
        core,projected,rejected = [],[],[]
        print('Building list of unrelated core samples.')
        # loop through finngen samples
        sample_iterator = basic_iterator(args.sample_fam,count= True,columns = 0)
        for i,sample in sample_iterator:
            progressBar(i,n_samples)
            # if it's flagged as a non finn, he's out regardless of all the other flags
            if sample in all_outliers:
                rejected.append(sample)
            # if it's a finn and it's not related, then we keep the sample!
            elif sample not in related_samples:
                core.append(sample)
            #if it's a finn and it's related, then it it's a outlier/rejected depending on whether it's a duplicate
            elif sample in duplicates:
                # this means that if the other duplicate has been included in the core sample, their duplicate will be removed.
                rejected.append(sample)
            else:
                projected.append(sample)
        print('done.')
        for entry in [(unrelated_file,core),(related_file,projected),(rejected_file,rejected),(final_samples,core+projected)]:
            write_fam_samplelist(*entry)           
            
    else:
        pass

    assert np.sum([mapcount(unrelated_file),mapcount(related_file),mapcount(rejected_file)]) == n_samples
    print(f'PCA inliers (unrelated finns) : {mapcount(unrelated_file)}')
    print(f'PCA outliers (related finns) : {mapcount(related_file)}')
    print(f'PCA rejected (non finns and duplicates) : {mapcount(rejected_file)}')
    print(f'PCA final samples : {mapcount(final_samples)}')

    return unrelated_file,related_file,rejected_file,final_samples

def fast_pca_inliers(args,core_samples):
    '''
    Runs a fast pca based on inliers only.
    '''
    #ouput path for PCA and outlier detection
    args.pca_output_file = args.pca_path + args.name

    eigenvecs =  args.pca_output_file+ '.eigenvec'   
    if not os.path.isfile(eigenvecs) or args.force:
        args.force = True 
        cmd = f'plink2  --bfile {args.merged_plink_file} --keep {core_samples}  --read-freq {args.merged_plink_file}.afreq --pca  {args.pca_components}  approx biallelic-var-wts -out {args.pca_output_file}  --threads {args.cpus} --memory {int(mem_mib)}'
        print(cmd)
        subprocess.call(shlex.split(cmd))

    else:
        args.logging.info('PCA for inliers already calculated')

    return eigenvecs

def project_all(args,eigenvecs):
    '''
    Projects inliers and outliers on space generated by the fast_pca
    '''

    # deinfe projection command
    cmd_start =f'plink2 --bfile {args.merged_plink_file} --threads {args.cpus} --memory {int(mem_mib)}'
    cmd_args = f' --read-freq {args.merged_plink_file}.afreq --score {eigenvecs}.var 2 3 header-read no-mean-imputation variance-standardize --score-col-nums 5-{3+args.pca_components+1} --out  {args.pca_path + args.name}'

    #define projection outputs
    projected_pca_file = os.path.join(args.pca_path,args.name+ '_projected_proj.sscore')
    core_pca_file = os.path.join(args.pca_path,args.name+ '_core_proj.sscore')
    if not os.path.isfile(core_pca_file) or not os.path.isfile(projected_pca_file ) or args.force:
        args.force = True
        cmd = f'{cmd_start} --keep {args.related_file} {cmd_args + "_projected_proj"}'
        print(cmd)
        subprocess.call(shlex.split(cmd))
                
        cmd = f'{cmd_start} --keep {args.unrelated_file}  {cmd_args + "_core_proj"}'
        print(cmd)
        subprocess.call(shlex.split(cmd))
    else:
        args.logging.info('projection already calculated')


    # merge the files to have the final PCs
    final_eigenvec = os.path.join(args.pca_path, args.name + '_eigenvec.txt')
    if not os.path.isfile(final_eigenvec) or args.force:
        args.force = True
        
        # read header of projection file to get the columns of the pcs
        proj_header = return_header(core_pca_file)
        pc_cols = [0,1] + [i for i,s in enumerate(proj_header) if s.startswith("PC")]
        with open(eigenvecs) as i: original_header = i.readline()
        # merge both projections
        with open(final_eigenvec,'wt') as o:
            o.write(original_header)
            rejected = np.loadtxt(args.rejected_file,usecols = 0,dtype = str)
            for f in [projected_pca_file,core_pca_file]:
                sample_iterator = basic_iterator(f,skiprows =1 ,columns = pc_cols)
                for entry in sample_iterator:
                    sample = entry[0]
                    assert sample not in rejected
                    o.write('\t'.join(entry) +'\n')
        subprocess.call(shlex.split(f"mv { args.pca_output_file+'.eigenvec.var'} {args.pca_output_file+'_eigenvec.var'}"))
        
    else:
        args.logging.info('final eigenvec file already generated')

    return final_eigenvec



