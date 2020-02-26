from utils import basic_iterator,np,mapcount,subprocess,return_header,tmp_bash,merge_files,make_sure_path_exists,return_header,identify_separator,write_fam_samplelist
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import chi2
import os,shlex

etnoDownload = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.txt'

def finn_or_not(args):
    '''
    Determines if remaining samples belong to the FIN or EUR cluster.
    '''
    args.eur_outlier_path = os.path.join(args.pca_outlier_path, "finngen_pca/")
    make_sure_path_exists(args.eur_outlier_path)
    args.false_finns = os.path.join(args.pca_outlier_path,args.name + '_false_finns.txt')

    #######################################
    # RETURN LIST OF EUR/FIN SAMPLES LEFT #
    #######################################    
    # euro samples that survived two rounds
    eur = args.eur_outlier_path + 'eur.txt'
    # fin samples that survived two rounds
    fin = args.eur_outlier_path + 'fin.txt'
    # finngen samples that survived two rounds
    finngen = args.eur_outlier_path + 'finngen.txt'
    if not os.path.isfile(eur) or not os.path.isfile(fin) or not os.path.isfile(finngen) or args.force:
        args.force = True
        batches = np.loadtxt(args.batches,dtype =str)
        # load summary statistics from the last round of outlier detection
        tag = 'first_round' if args.outlier_pca_rounds == 1 else 'second_round'
        last_round_samples = os.path.join(args.pca_outlier_path, tag ,  tag + '_' + args.name + '_outlier_samples.tsv')
        cols = [return_header(last_round_samples).index(elem) for elem in ['IID','outlier','SuperPops']]
        # keep inliers of latest round 
        outlier_iterator = basic_iterator(last_round_samples,columns = cols)  
        euros,finns,finngens =[],[],[]
        for entry in outlier_iterator:
            iid,outlier,superpop = entry
            if superpop in batches:
                if outlier == 'FALSE':
                    finngens.append(iid)               
            elif superpop == "EUR":euros.append(iid)                      
            elif superpop =="FIN":finns.append(iid)
                
        for entry in [(eur,euros),(fin,finns),(finngen,finngens)]:
            write_fam_samplelist(*entry)
            
    if os.path.getsize(eur) == 0:
        #If there are no EUR left, we save the outliers of the first two rounds as false_finns and exit
        outliers = np.loadtxt(args.ethnic_pca_outliers,dtype = str,usecols = [0])
        write_fam_samplelist(args.false_finns,outliers)
        return
    
    #######################################################################
    # PCA with finngen final samples and project everyone onto that space #
    #######################################################################
    pca_output_file = args.eur_outlier_path + 'finngen_' + args.name
    if not os.path.isfile(pca_output_file+ '.eigenval') or args.force:
        print('Calculating finngen pca...')
        args.force = True
        cmd = f'plink2 --bfile {args.tg_merged_plink_file}  --keep {finngen}  --read-freq {args.tg_merged_plink_file + ".frq"} --pca {args.pca_components}  approx biallelic-var-wts --threads {args.cpus} -out {pca_output_file}'
        subprocess.call(shlex.split(cmd))       
          
    else:
        args.v_print(3,'finngen pca already calculated.')
  
    #project finngen,eur,fin on space generated by finngen samples
    columns = range(5,5+args.pc_filter)
    for tag,f in [('finngen',finngen),('eur',eur),('fin',fin)]:
        project(args,tag,f,columns,pca_output_file)

    ###################
    # EUR VS FIN PROB #
    ###################
    
    # return euro outliers based on their projection
    finngen_eur_outliers_path = os.path.join(args.pca_outlier_path, 'finngen_eur_outliers.txt')
    if not os.path.isfile(finngen_eur_outliers_path) or args.force:
        args.force = True
        eur_outliers = fin_eur_probs(args.eur_outlier_path,args.pc_filter,args.finn_prob_filter)
        write_fam_samplelist(finngen_eur_outliers_path,eur_outliers)
    else:
        pass


    print(f'eur outliers : {mapcount(finngen_eur_outliers_path)}')

    # merge list of bayes and eur centroid outliers
    merge_files(args.false_finns,[finngen_eur_outliers_path,args.ethnic_pca_outliers])
    print(f'total false finns : {mapcount(args.false_finns)}')
   
  
#################
#-- BASIC PCA --#
#################
def outlier_pca(args):
    '''
    Performs single or double round pca outlier detection.
    Returns the final list of outliers in .fam format
    '''

    outliers = pca_round(args,'first_round')
    if args.outlier_pca_rounds > 1:
        # if requested, do another round
        outliers += pca_round(args,'second_round',remove_list = outliers)
    
    args.ethnic_pca_outliers = os.path.join(args.pca_outlier_path, args.name + '_ethnic_outliers.txt')
    if not os.path.os.path.isfile(args.ethnic_pca_outliers) or args.force:
        args.force = True 
        # remove tg samples from list of outliers and save in _ethnic outliers
        tg_sample_list = np.loadtxt(args.new_tg +'.fam',dtype = str,usecols = 1)
        false_finngen = [sample for sample in outliers if sample not in tg_sample_list]
        write_fam_samplelist(args.ethnic_pca_outliers,false_finngen)
    
    else:
        pass
    print(f'total finngen ethnic outliers : {mapcount(args.ethnic_pca_outliers)}' )
                                    
        
def pca_round(args,tag,remove_list = None):
    '''
    Performs PCA + outlier_detection and returns the list of outliers
    '''
    local_path=  os.path.join(args.pca_outlier_path, tag+'/')
    make_sure_path_exists(local_path)
    remove = ''
    if remove_list is not None:
        remove_file = tmpPath + tag + '_remove.fam'
        write_fam_samplelist(remove_file,remove_list)
        remove =  f' --remove {remove_file}'

    #######
    # PCA #
    #######
    pca_output_file = os.path.join(local_path,tag +'_' + args.name)
    if not os.path.isfile( pca_output_file+ '.eigenval') or args.force:
        args.force = True 
        #individuals that need to be removed
        print(remove)
        cmd = f'plink2 --bfile {args.tg_merged_plink_file} --read-freq  {args.tg_merged_plink_file}.frq {remove} --pca {args.pca_components} approx biallelic-var-wts --threads {args.cpus}  -out {pca_output_file}'
        print(cmd)
        subprocess.call(shlex.split(cmd))
        
    #####################
    # OUTLIER DETECTION #
    #####################
    outlier_samples = pca_output_file +'_outlier_samples.tsv'
    iterations = ' 1000 -p' if args.test else ' 3000 '

    if not os.path.isfile(outlier_samples) or args.force:
        args.force = True 
        print('generating outliers at ' + pca_output_file)
        cmd = f' -f {pca_output_file+".eigenvec"} -e {pca_output_file+".eigenval"} -s {args.annot_pop} --n_iterations {iterations}  -o {pca_output_file}'
        tmp_cmd  =f"Rscript {os.path.join(args.rootPath,'scripts/pca_outlier_detection/scripts/classify_outliers.R')} {cmd} "
        subprocess.call(shlex.split(tmp_cmd))

           
    #return outliers
    outlier_file = pca_output_file + '_outliers.txt'
    outliers = np.genfromtxt(outlier_samples, dtype=str, usecols=(0, 1), delimiter='\t')
    out_mask = (outliers[:,1] =='TRUE')
    outliers = outliers[:,0][out_mask]
    print(f'{tag} total outliers : {len(outliers)}')
    
    return list(outliers)



def build_superpop(args):
    '''
    Writes pop info for outlier detection
    '''
    
    # sample to super pop mapping file
    args.annot_pop =  os.path.join(args.pca_outlier_path, args.name + '_sample_annot.tsv')
    if not os.path.os.path.isfile(args.annot_pop) or args.force or mapcount(args.annot_pop) < 2:
        args.force = True 
        print('generating super pop dict')

        # pop to superpop mapping
        pop_iterator = basic_iterator(args.data_path + 'superpop.csv',',') #pop to superpop dict
        superpop_dict=  {pop:superpop for pop,superpop in pop_iterator}
        print(superpop_dict)

        # sample info for tg
        tg_pop = os.path.join(args.data_path,'20130606_sample_info.txt')
        if not os.path.isfile(tg_pop):
            cmd = f'wget {etnoDownload} -P {args.data_path}'
            subprocess.call(shlex.split(cmd))

        # get pop info of each samples
        cols = [return_header(tg_pop).index(elem) for elem in ['Sample','Population']]
        pop_dict = {sample:pop for sample,pop in basic_iterator(tg_pop,columns = cols)}
        with open(args.annot_pop,'wt') as o:
            o.write('IID' +'\t' + 'SuperPops\n')
            #write sample to superpop/batch mapping for 1k samples
            tg_sample_iterator = basic_iterator(args.new_tg +'.fam',columns = 1,separator = " ")
            for sample in tg_sample_iterator:
                pop = pop_dict[sample] # fetch country
                super_pop = pop if pop =='FIN'else superpop_dict[pop]
                o.write('\t'.join([sample,super_pop]) +'\n')
            #add batch info from FINNGEN samples
            cols = [return_header(args.sample_info).index(elem) for elem in ['IID','BATCH']]
            for sample,batch in basic_iterator(args.sample_info,columns =cols,separator =','):
                o.write('\t'.join([sample,batch]) +'\n')
                  
    else:
        args.v_print(3,'super pop dict already generated.')
        
def fin_eur_probs(eur_outlier_path,pc_filter,finn_prob_filter):

    '''
    Returns probability of being a FIN based on mahalanobis distance to EUR and FIN centroids of finngen samples
    '''
    
    # EUR
    eur_data = np.loadtxt(eur_outlier_path+ 'eur.eigenvec',dtype = float, skiprows = 1,usecols = range(1,pc_filter+1))
    eur_avg = np.reshape(np.average(eur_data,axis =0),(1,eur_data.shape[1]))
    eur_cov = np.linalg.inv(np.cov(eur_data.T))
    # FIN
    fin_data = np.loadtxt(eur_outlier_path+ 'fin.eigenvec',dtype = float, skiprows = 1,usecols = range(1,pc_filter+1))      
    fin_avg = np.reshape(np.average(fin_data,axis = 0),(1,fin_data.shape[1]))
    fin_cov = np.linalg.inv(np.cov(fin_data.T))
    # FINNGEN
    finngen_data = np.loadtxt(eur_outlier_path+ 'finngen.eigenvec',dtype = float, skiprows = 1,usecols = range(1,pc_filter+1))
    # MAHALANOBIS DISTANCES
    eur_dist = cdist(finngen_data,eur_avg,metric = 'mahalanobis',VI = eur_cov).flatten()**2
    fin_dist = cdist(finngen_data,fin_avg,metric = 'mahalanobis',VI = fin_cov).flatten()**2
    # PROBS
    p_eur = 1 -  chi2.cdf(eur_dist,pc_filter)
    p_fin =  1 - chi2.cdf(fin_dist,pc_filter)
    f_prob = p_fin/(p_fin + p_eur)
    fin_mask = (f_prob < finn_prob_filter)
    print(f'outliers: {fin_mask.sum()}')
    eur_outliers = np.loadtxt(eur_outlier_path+ 'finngen.eigenvec',dtype = str, skiprows = 1,usecols = 0 )[fin_mask]

    with open(eur_outlier_path + 'fin_bins.txt','wt') as o:
        fin_probs = f_prob[fin_mask]
        bins = np.linspace(0, 1, 11)
        digitized = np.digitize(fin_probs, bins)
        for j,bin_value in enumerate([len(fin_probs[digitized == i]) for i in range(1, len(bins))]):
            o.write(f"{round(bins[j],2)}\t{bin_value}\n")
        
    return eur_outliers


      
def project(args,tag,samples_keep,columns,pca_output_file):
    '''
    Project eur/fin/finngen samples onto the space generated by finngen data.
    '''

    # define projection command
    columns = list(map(str,columns))
    cut_columns =  ','.join(columns)
        
    if not os.path.isfile(args.eur_outlier_path+ tag + '.sscore' ) or args.force:
        args.force = True 
        print(tag)
        keep = f' --keep {samples_keep}'
        cmd = f'plink2 --bfile {args.tg_merged_plink_file} {keep} --score {pca_output_file+".eigenvec.var"} 2 3 header-read no-mean-imputation  --score-col-nums {columns[0]}-{columns[-1]} --out {args.eur_outlier_path+tag}'
        subprocess.call(shlex.split(cmd))
        cmd = f'cat {args.eur_outlier_path+tag + ".sscore"}  | cut -f2,{cut_columns} >  {args.eur_outlier_path+ tag + ".eigenvec"}'
        tmp_bash(cmd)
        
