import argparse,os,shlex,subprocess,sys,logging,multiprocessing,time
from utils import pretty_print,file_exists,make_sure_path_exists,log_levels,Logger
from new_scripts import batches,tg,ethnic_outliers
from pathlib import Path
from pca_scripts import plot

def main(args):

    # Get batch info
    #pretty_print("BATCHES & SAMPLES")
    #args.sample_info,args.cohorts,args.batches = batches.batches(args)

    #Merge bed and 1k genome
    pretty_print('1k GENOME')
    args.plink_path = os.path.join(args.out_path,'plink_files/')
    make_sure_path_exists(args.plink_path)
    args.merged_plink_file = tg.merge_1k(args)


    # PCA for genetic outliers
    pretty_print('ETHNIC FINNS') 
    args.pca_outlier_path = os.path.join(args.out_path, 'outliers_pca/')
    make_sure_path_exists(args.pca_outlier_path)
    args.annot_pop,args.fg_tags = ethnic_outliers.build_superpop(args)
    aberrant_output,args.ethnic_outliers = ethnic_outliers.detect_ethnic_outliers(args)
    args.finngen_eur_outliers = ethnic_outliers.detect_eur_outliers(args,aberrant_output)

    #PLOT
    args.plot_path = os.path.join(args.out_path,'plots')
    outlier_plot_data = os.path.join(args.plot_path,'plot_data')
    #make_sure_path_exists(outlier_plot_data)
    #make_sure_path_exists(args.plot_path)
    #plot.plot_first_round_outliers(args)
    #plot.plot_fin_eur_outliers(args)

    
    return True
    
if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description="FinnGen PCA pipeline.")

    # BED FILES
    parser.add_argument("--bed", type=file_exists, help = "Folder in which the merged plink file is stored", required = True)
    parser.add_argument('--tg-bed',type = file_exists,help = 'Plink 1k file')

    #SAMPLE DATA
    parser.add_argument("--sample-info", type=file_exists, help =  "Path to csv file with sample,batch", required = True)


    # NUMERIC PARAMS
    parser.add_argument('--pca-components',type=int,help='Components needed for pca',default = 20)
    parser.add_argument('--pc-filter',type = int,help = 'Number of pcs on which to perform the outlier detection method',default = 3)
    parser.add_argument('--finn-prob-filter',type = float,help = 'Filter falue to decide whether a finngen sample is part of the EUR or FIN centroid',default = 0.95)


    # GENERAL PARAMS
    parser.add_argument('-o',"--out_path",type = str, help = "Folder in which to save the results", required = True)
    parser.add_argument('--name',type = str,default = 'test',help = 'Name to append to output files')

    # LOGGING ET AL.
    parser.add_argument( "-log",  "--log",  default="warning", choices = log_levels, help=(  "Provide logging level. " "Example --log debug', default='warning'"))
    parser.add_argument('--force',action = 'store_true',help = 'Replaces files by force',default = False)
    parser.add_argument('--test',action = "store_true",help = 'Flag for quick pca_outlier method without plots.')
    parser.add_argument("--cpus",type = int, help = "Number of cpus to use (default available cpus)", default =  multiprocessing.cpu_count())



    args = parser.parse_args()    

    # logging level
    level = log_levels[args.log]
    logging.basicConfig(level=level,format="%(levelname)s: %(message)s")
    args.logging = logging

    args.parent_path = Path(os.path.realpath(__file__)).parent.parent
    args.data_path = os.path.join(args.parent_path,'data/')
    args.misc_path =  os.path.join(args.out_path, 'misc/')
    make_sure_path_exists([args.out_path,args.misc_path,args.data_path])
    args.sample_fam = args.bed.replace('.bed','.fam')

    args.success = False
    success = main(args)
    
    args.log_file =os.path.join(args.out_path ,args.name + '.log' )
    if success:
        args.logging.getLogger().setLevel(logging.INFO)
        with  Logger(args.log_file,'wt'):
            args.force = False
            #success = main(args)        

