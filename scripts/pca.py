import argparse,os,logging,multiprocessing,shutil
from utils import pretty_print,file_exists,make_sure_path_exists,log_levels,Logger,print_msg_box,get_filepaths,mapcount
from pca_scripts import tg,ethnic_outliers,kinship,pca,plot
from pathlib import Path

def main(args,do_plot=True):


    #Merge bed and 1k genome
    pretty_print('1k GENOME')
    args.plink_path = os.path.join(args.out_path,'plink_files/')
    make_sure_path_exists(args.plink_path)
    args.merged_plink_file = tg.tg_bed(args)


    # PCA for genetic outliers
    pretty_print('ETHNIC FINNS') 
    args.pca_outlier_path = os.path.join(args.out_path, 'outliers_pca/')
    make_sure_path_exists(args.pca_outlier_path)
    args.annot_pop,args.fg_tags = ethnic_outliers.build_superpop(args)
    aberrant_output,args.ethnic_outliers = ethnic_outliers.detect_ethnic_outliers(args)
    args.finngen_eur_outliers = ethnic_outliers.detect_eur_outliers(args,aberrant_output)
    args.all_outliers = ethnic_outliers.all_outliers(args,[args.ethnic_outliers,args.finngen_eur_outliers])
    
    #KINSHIP DATA
    pretty_print('KINSHIP')
    args.kinPath = os.path.join(args.out_path,'kinship/')
    make_sure_path_exists(args.kinPath)
    kinship.kinship(args)

    #PCA
    pretty_print('PCA')
    args.pca_path = os.path.join(args.out_path,'pca/')
    make_sure_path_exists(args.pca_path)
    args.unrelated_file,args.related_file,args.rejected_file,args.final_samples = pca.build_inliers(args)
    core_eigenvec = pca.fast_pca_inliers(args,args.unrelated_file)
    args.eigenvec = pca.project_all(args,core_eigenvec)  

    if do_plot:
        #PLOT
        pretty_print("PLOT")
        args.plot_path = os.path.join(args.out_path,'plots')
        outlier_plot_data = os.path.join(args.plot_path,'plot_data')
        make_sure_path_exists(outlier_plot_data)
        make_sure_path_exists(args.plot_path)
        plot.plot_first_round_outliers(args)
        plot.plot_fin_eur_outliers(args)
        plot.plot_final_pca(args)
        plot.pca_map(args)
    
    return True


def release(args):

    import glob
    doc_path = os.path.join(args.out_path,'documentation')
    data_path = os.path.join(args.out_path,'data')
    for path in [data_path,doc_path]:
        make_sure_path_exists(path)
        for f in get_filepaths(path): os.remove(f) # clean path else shutil.copy might fail

    # DOC
    for pdf in glob.glob(os.path.join(args.plot_path,'*pdf')):
        shutil.copy(pdf,os.path.join(doc_path,os.path.basename(pdf)))
    outlier_pdf = os.path.join(args.pca_outlier_path, '1k_pca/',args.name + '_outlier_pcas.pdf')
    shutil.copy(outlier_pdf,os.path.join(doc_path,os.path.basename(outlier_pdf)))
    shutil.copy(args.log_file,os.path.join(doc_path,os.path.basename(args.log_file)))

    # DATA
    for f in [args.unrelated_file,args.related_file,args.rejected_file,args.final_samples,args.duplicates,args.all_outliers,args.eigenvec,args.pca_output_file + '.eigenval',args.pca_output_file + '_eigenvec.var']:
        shutil.copy(f,os.path.join(data_path,os.path.basename(f)))

    # README
    readme = os.path.join(args.data_path,'pca.README') 
    with open(os.path.join(args.out_path,args.name + '_pca_readme'),'wt') as o, open(readme,'rt') as i:
        with open(args.log_file) as tmp: summary = tmp.read()
        word_map = {'[PREFIX]':args.name,'[SUMMARY]': summary,'[N_SNPS]':mapcount(args.bed.replace('.bed','.bim'))}
        for line in i:
            for kw in word_map:
                if kw in line:
                    line = line.replace(kw,str(word_map[kw]))
            o.write(line)


    
if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description="FinnGen PCA pipeline.")

    # BED FILES
    parser.add_argument("--bed", type=file_exists, help = "Folder in which the merged plink file is stored", required = True)
    parser.add_argument('--tg-bed',type = file_exists,help = 'Plink 1k file')

    #SAMPLE DATA
    parser.add_argument("--sample-info", type=file_exists, help =  "Path to csv file with sample,batch", required = True)
    parser.add_argument("--meta", type=file_exists, help =  "Path to file with regionofbirth info", required = True)
    parser.add_argument("--tg-pop", type=file_exists, help =  "Path to file with regionofbirth info")

    #KINSHIP
    parser.add_argument('--degree',type=int,help='Degree for Kinship',default = 2)
    parser.add_argument("--kin", type=file_exists, help = "File with king related individuals")

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
    parser.add_argument('--release',action = 'store_true',help = 'Flag for data release',default = False)



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
    if not args.tg_pop:
        args.tg_pop = os.path.join(args.data_path,'1kg_meta.txt')
    args.success = False
    success = main(args)
    
    args.log_file =os.path.join(args.out_path ,args.name + '.log' )
    if success:
        args.logging.getLogger().setLevel(logging.WARNING)
        with  Logger(args.log_file,'wt'):
            args.force = False
            print_msg_box("\n~ SUMMARY ~\n",indent =30)
            success = main(args,False)        
            
        if args.release:
            release(args)

