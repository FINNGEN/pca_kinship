from utils import mapcount,basic_iterator,progressBar,merge_files,return_header,mem_mib,write_fam_samplelist
import numpy as np
import os.path,shlex,subprocess


def build_inliers(args):
    '''
    Returns .fam file for running pca by excluding samples that are marked as related or as false finns. Also, it returns the intersection of the false_finns and duplicates, so that they can be discarded.
    '''

    #LOAD ALL RELEVANT LIST OF SAMPLES
    related_samples = np.loadtxt(args.related_individuals,dtype = str,usecols = 0)
    false_finns = np.loadtxt(args.false_finns,dtype = str,usecols = 0)
    duplicates = np.loadtxt(args.duplicates,dtype = str,usecols = 0)
    # INLIERS
    args.inlier_file = os.path.join(args.pca_path,args.name +'_inliers.txt')
    args.outlier_file = os.path.join(args.pca_path,args.name +'_outliers.txt')
    args.rejected_file = os.path.join(args.pca_path,args.name +'_rejected.txt')
    args.final_samples = os.path.join(args.pca_path,args.name +'_final_samples.txt')
    if not all(map(os.path.isfile,[args.inlier_file,args.outlier_file,args.rejected_file,args.final_samples])) or args.force:
        args.force = True 
        inliers,outliers,rejected = [],[],[]
        print('Building list of inliers')
        # loop through finngen samples
        sample_iterator = basic_iterator(args.sample_fam,count= True,columns = 0)
        for i,sample in sample_iterator:
            progressBar(i,args.n_samples)
            # if it's flagged as a false finn, he's out regardless of all the other flags
            if sample in false_finns:
                rejected.append(sample)
            # if it's not a false finn and it's not related, then we keep the sample!
            elif sample not in related_samples:
                inliers.append(sample)
            #if it's not a false finn and it's related, then it it's a outlier/rejected depending on whether it's a duplicate
            elif sample in duplicates:
                rejected.append(sample)
            else:
                outliers.append(sample)

        for entry in [(args.inlier_file,inliers),(args.outlier_file,outliers),(args.rejected_file,rejected),(args.final_samples,inliers+outliers)]:
            write_fam_samplelist(*entry)           

    else:
        pass

    assert np.sum([mapcount(args.inlier_file),mapcount(args.outlier_file),mapcount(args.rejected_file)]) == args.n_samples
    print(f'PCA inliers (unrelated finns) : {mapcount(args.inlier_file)}')
    print(f'PCA outliers (related finns) : {mapcount(args.outlier_file)}')
    print(f'PCA rejected (non finns and duplicates) : {mapcount(args.rejected_file)}')
    print(f'PCA final samples : {mapcount(args.final_samples)}')
        
def fast_pca_inliers(args):
    '''
    Runs a fast pca based on inliers only.
    '''
    #ouput path for PCA and outlier detection
    args.pca_output_file = args.pca_path + args.name
    
    if not os.path.isfile( args.pca_output_file+ '.eigenval') or args.force:
        args.force = True 
        keep =  ' --keep ' + args.inlier_file
        cmd = f'plink2  --bfile {args.merged_plink_file}  {keep}  --read-freq {args.merged_plink_file}.afreq --pca  {args.pca_components}  approx biallelic-var-wts -out {args.pca_output_file}  --threads {args.cpus} --memory {int(mem_mib)}'
        print(cmd)
        subprocess.call(shlex.split(cmd))

    else:
        args.v_print(3,'PCA for inliers already calculated')


def project_all(args):
    '''
    Projects inliers and outliers on space generated by the fast_pca
    '''

    # deinfe projection command
    cmd_start =f'plink2 --bfile {args.merged_plink_file} --threads {args.cpus} --memory {int(mem_mib)}'
    cmd_args = f' --read-freq {args.merged_plink_file}.afreq --score {args.pca_output_file+".eigenvec.var"} 2 3 header-read no-mean-imputation variance-standardize --score-col-nums 5-{3+args.pca_components+1} --out  {args.pca_path + args.name}'

    #define projection outputs
    outliers_pca_file = os.path.join(args.pca_path,args.name+ '_outliers_proj.sscore')
    inliers_pca_file = os.path.join(args.pca_path,args.name+ '_inliers_proj.sscore')
    if not os.path.isfile(inliers_pca_file) or not os.path.isfile(outliers_pca_file ) or args.force:
        args.force = True
        cmd = f'{cmd_start} --keep {args.outlier_file} {cmd_args + "_outliers_proj"}'
        print(cmd)
        subprocess.call(shlex.split(cmd))
                
        cmd = f'{cmd_start} --keep {args.inlier_file}  {cmd_args + "_inliers_proj"}'
        print(cmd)
        subprocess.call(shlex.split(cmd))
    else:
        args.v_print(3,'inliers/outliers projected')
    # merge the files
    args.eigenvec = os.path.join(args.pca_path, args.name + '_final.eigenvec')
    if not os.path.isfile(args.eigenvec) or args.force:
        args.force = True
        
        # read header of projection file to get the columns of the pcs
        proj_header = return_header(inliers_pca_file)
        pc_cols = [0,1] + [i for i,s in enumerate(proj_header) if s.startswith("PC")]
        # read header of eigenvec file
        with open(args.pca_output_file + '.eigenvec','rt') as i: original_header = i.readline()
        # merge inliers and outliers file
        print(args.eigenvec)
    
        with open(args.eigenvec,'wt') as o:
            o.write(original_header)
            rejected = np.loadtxt(args.rejected_file,usecols = 0,dtype = str)
            for f in [outliers_pca_file,inliers_pca_file]:
                sample_iterator = basic_iterator(f,skiprows =1 ,columns = pc_cols)
                for entry in sample_iterator:
                    sample = entry[0]
                    assert sample not in rejected
                    o.write('\t'.join(entry) +'\n')
               
    else:
        args.v_print(3,'final eigenvec file already generated')
        return


