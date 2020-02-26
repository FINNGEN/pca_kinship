from utils import pretty_print,file_exists,make_sure_path_exists
import multiprocessing,glob,argparse,os.path,subprocess,shlex
from pca import batches,tg,kinship,true_finns

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
    true_finns.build_superpop(args)
    true_finns.outlier_pca(args)
    true_finns.finn_or_not(args)

    
    pretty_print('KINSHIP')
    args.kinPath = os.path.join(args.out_path,'kinship/')
    make_sure_path_exists(args.kinPath)
    kinship.kinship(args)


if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description="Returning final list of variants after info_score filter and ld pruning")

    parser.add_argument('-b',"--bed", type=file_exists, help = "Folder in which the merged plink file is stored", required = True)
    parser.add_argument('-o',"--out_path",type = str, help = "folder in which to save the results", required = True)
    parser.add_argument('--name',type = str,default = 'test',help = 'Name to append to output files ')
    parser.add_argument('-s',"--sample-info", type=file_exists, help =  "Path to csv file with sample,batch", required = True)

    #KINSHIP
    parser.add_argument('--degree',type=int,help='Degree for Kinship',default = 2)
    parser.add_argument('-k',"--kin-file", type=file_exists, help = "File with degree 3 kinship")

    #PCA
    parser.add_argument('--pca-components',type=int,help='Components needed for pca',default = 20)

    #1k
    parser.add_argument('--tg-bed',type = file_exists,help = 'Plink 1k file')

    #OUTLIER DETECTION
    parser.add_argument('--outlier-pca-rounds',type = int,help = 'How many roudns of PCA to perform in outlier detection',default =1,choices = (1,2) )
    parser.add_argument('--finn-prob-filter',type = float,help = 'Filter falue to decide whether a finngen sample is part of the EUR or FIN centroid',default = 0.95)
    parser.add_argument('--pc-filter',type = int,help = 'Number of pcs on which to perform the outlier detection method',default = 3)
    
    #optional tags
    parser.add_argument('--exclude',type = file_exists,help ='list of snps to exclude',default = None)
    parser.add_argument('--test',type = int,default =0,help = 'Flag for quick pca_outlier. For testing purposes. It only keeps 10k samples.')
    parser.add_argument('--plot',action = 'store_true',help = 'Flag for plotting')
    parser.add_argument("--cpus",type = int, help = "number of cpus to use", default =  multiprocessing.cpu_count())
    parser.add_argument('--force',action = 'store_true',help = 'Replaces files by force',default = False)
    parser.add_argument('-v', '--verbosity', action="count", 
                        help="increase output verbosity (e.g., -vv is more than -v)")
 
    args=parser.parse_args()

    # creating outputpath
    make_sure_path_exists(args.out_path)

    # write log file
    args.results_path = os.path.join(args.out_path,'results/')
    make_sure_path_exists(args.results_path)
    log_file =os.path.join(args.results_path ,args.name + '.log' )
    with open(log_file,'wt') as o:
        d = vars(args)
        for key in d:
            o.write(f'{key}:{d[key]} \n')   
    
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
    main(args)
    if args.success:
        print('Pipeline ran successfully, logging...')
        args.force = False
        import contextlib
        with open(log_file,'a') as f:
            with contextlib.redirect_stdout(f):
                main(args)


