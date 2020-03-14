from utils import pretty_print,file_exists,make_sure_path_exists,get_filepaths,merge_files,mapcount,Logger
import multiprocessing,glob,argparse,os.path,subprocess,shlex,shutil,sys
from pca_scripts import batches,tg,kinship,pca,plot,true_finns
def main(args):

    
    # Get batch info
    pretty_print("BATCHES & SAMPLES")
    args.plink_path = os.path.join(args.out_path,'plink_files/') 
    make_sure_path_exists(args.plink_path)
    batches.batches(args)

    #downloads 1k genome project and builds a plink file that can be merged
    pretty_print('1k GENOME')
    tg.subset_1k(args)
    tg.merge_1k(args)
        
    # PCA for outliers
    pretty_print('TRUE FINNS') 
    args.pca_outlier_path = os.path.join(args.out_path, 'outliers_pca/')
    make_sure_path_exists(args.pca_outlier_path)
    true_finns.outlier_pca(args)
    true_finns.finn_or_not(args)
    true_finns.final_merge(args)
    
    #KINSHIP DATA
    pretty_print('KINSHIP')
    args.kinPath = os.path.join(args.out_path,'kinship/')
    make_sure_path_exists(args.kinPath)
    kinship.kinship(args)

    #RUN PCA
    pretty_print('PCA')
    args.pca_path = os.path.join(args.out_path,'pca/')
    make_sure_path_exists(args.pca_path)
    pca.build_inliers(args)
    pca.fast_pca_inliers(args)
    pca.project_all(args)  

    #PLOT
    args.plot_path = os.path.join(args.out_path,'plots')
    make_sure_path_exists(args.plot_path)
    plot.plot_final_pca(args)
    plot.plot_first_round_outliers(args)
    plot.plot_fin_eur_outliers(args)
    plot.plot_map(args)

    return True
def release(args):

    import glob
    doc_path = os.path.join(args.out_path,'documentation')
    make_sure_path_exists(doc_path)
    for f in get_filepaths(doc_path): os.remove(f) # clean path else shutil.copy might fail
    for pdf in glob.glob(os.path.join(args.plot_path,'*pdf')):
        shutil.copy(pdf,os.path.join(doc_path,os.path.basename(pdf)))
    outlier_pdf = os.path.join(args.pca_outlier_path, '1k_pca/',args.name + '_outlier_pcas.pdf')
    shutil.copy(outlier_pdf,os.path.join(doc_path,os.path.basename(outlier_pdf)))
    shutil.copy(args.log_file,os.path.join(doc_path,os.path.basename(args.log_file)))
        
    data_path = os.path.join(args.out_path,'data')
    make_sure_path_exists(data_path)
    for f in get_filepaths(data_path): os.remove(f) # clean path else shutil.copy might fail
    for f in [args.inlier_file,args.outlier_file,args.rejected_file,args.final_samples,args.duplicates,args.false_finns,args.eigenvec,args.pca_output_file + '.eigenval',args.pca_output_file + '_eigenvec.var']:
        shutil.copy(f,os.path.join(data_path,os.path.basename(f)))
    with open(args.log_file) as i:
        summary = i.read()
        
    readme = os.path.join(args.data_path,'pca.README') 
    with open(os.path.join(args.out_path,args.name + '_pca_readme'),'wt') as o, open(readme,'rt') as i:
        word_map = {'[PREFIX]':args.name,'[SUMMARY]': summary,'[N_SNPS]':mapcount(args.bed.replace('.bed','.bim'))}
        for line in i:
            for kw in word_map:
                if kw in line:
                    line = line.replace(kw,str(word_map[kw]))
            o.write(line)

if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description="Returning final list of variants after info_score filter and ld pruning")

    parser.add_argument("--bed", type=file_exists, help = "Folder in which the merged plink file is stored", required = True)
    parser.add_argument('-o',"--out_path",type = str, help = "Folder in which to save the results", required = True)
    parser.add_argument('--name',type = str,default = 'test',help = 'Name to append to output files ')
    parser.add_argument("--sample-info", type=file_exists, help =  "Path to csv file with sample,batch", required = True)
    parser.add_argument("--meta", type=file_exists, help =  "Path to file with regionofbirth info", required = False)

    #KINSHIP
    parser.add_argument('--degree',type=int,help='Degree for Kinship',default = 2)
    parser.add_argument("--kin", type=file_exists, help = "File with king related individuals")

    #PCA
    parser.add_argument('--pca-components',type=int,help='Components needed for pca',default = 20)

    #1k
    parser.add_argument('--tg-bed',type = file_exists,help = 'Plink 1k file')

    #OUTLIER DETECTION
    parser.add_argument('--finn-prob-filter',type = float,help = 'Filter falue to decide whether a finngen sample is part of the EUR or FIN centroid',default = 0.95)
    parser.add_argument('--pc-filter',type = int,help = 'Number of pcs on which to perform the outlier detection method',default = 3)
    
    #optional tags
    parser.add_argument('--test',type = int,default =0,help = 'Flag for quick pca_outlier. For testing purposes. It only keeps 10k sample.')
    parser.add_argument("--cpus",type = int, help = "Number of cpus to use (default available cpus)", default =  multiprocessing.cpu_count())
    parser.add_argument('--force',action = 'store_true',help = 'Replaces files by force',default = False)
    parser.add_argument('--release',action = 'store_true',help = 'Flag for data release',default = False)
    parser.add_argument('-v', '--verbosity', action="count", 
                        help="Increase output verbosity (e.g., -vv is more than -v)")
 
    args=parser.parse_args()
    args.name = args.name
    make_sure_path_exists(args.out_path)
  
    args.rootPath = '/'.join(os.path.realpath(__file__).split('/')[:-2]) + '/'
    args.data_path = os.path.join(args.rootPath,'data/')
    args.misc_path =  os.path.join(args.out_path, 'misc/')
    make_sure_path_exists(args.misc_path)

    if args.verbosity:
        def _v_print(*verb_args):
            if verb_args[0] > (3 - args.verbosity):
                print(verb_args[1])
    else:
        _v_print = lambda *a: None  # do-nothing function

    args.v_print = _v_print
         
    args.success = False
    success = main(args)
    
    args.log_file =os.path.join(args.out_path ,args.name + '.log' )
    if success:
        with  Logger(args.log_file,'wt'):
            args.force = False
            success = main(args)        
        
        if args.release:
            release(args)
