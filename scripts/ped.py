import os,pickle, subprocess,shlex,argparse,multiprocessing,shutil
from utils import basic_iterator,return_header,mapcount,get_path_info,file_exists,make_sure_path_exists,cpus,tmp_bash,pretty_print,identify_separator,NamedTemporaryFile,get_filepaths
from collections import defaultdict as dd
import numpy as np
import pandas as pd
cpus = multiprocessing.cpu_count()

######################
#---BUILD BED FILE---#
######################
def build_bed(args):
    """ 
    Builds bed with hq variants for kinship. This is the core data set used for the all kinship analysis.    
    """

 
    args.kinship_bed =os.path.join(args.plink_path,args.prefix + '_kinship.bed')
    if not os.path.isfile(args.kinship_bed) or not mapcount(args.kinship_bed) or args.force:
        args.force = True
        # filter for only 30k samples 
        keep = f"--keep {args.fam}" if args.fam and mapcount(args.fam) > 100 else ""
        extract = f"--extract {args.extract}" if args.extract else ""
        cmd = f"plink2 --bfile {args.bed.replace('.bed','')} {extract} --threads {cpus}  --make-bed --out {args.kinship_bed.replace('.bed','')} {keep}"
        print(cmd)
        subprocess.call(shlex.split(cmd))
    else:
        print('bed file already generated')
    freq_file = args.kinship_bed.replace('.bed','.afreq')
    if not os.path.isfile(freq_file) or not mapcount(freq_file) :
        cmd = f"plink2 --bfile {args.kinship_bed.replace('.bed','')} --freq  --threads {cpus}  --out {args.kinship_bed.replace('.bed','')}"
        print(cmd)
        subprocess.call(shlex.split(cmd))


######################
#------KINSHIP-------#
######################
def kinship(args):
    """
    Returns degree 3 kinship data.
    """

    args.kinship_log_file = os.path.join(args.out_path,args.prefix + '_kinship.log')  
    args.kin_file = os.path.join(args.kinship_path,f"{args.prefix}.kin0")
    # RETURN RELATED AND PLOT FAMILIES
    if not os.path.isfile(args.kin_file) or mapcount(args.kin_file) < 1 or args.force:
        args.force = True
        cmd = f'king --cpus {cpus} -b {args.kinship_bed} --related --duplicate --degree 2 --prefix {os.path.join(args.kinship_path,args.prefix)} --rplot '
        print(cmd)
        with open(args.kinship_log_file,'at') as f: subprocess.call(shlex.split(cmd),stdout = f)
        # produce R scripts, fix them and run them
    else:
        print("related file already generated")
         
        
    if args.force:
        scriptFile = NamedTemporaryFile(delete=True)
        for f in [f for f in get_filepaths(args.kinship_path) if f.endswith('.R')]:
            file_path,file_root,file_extension = get_path_info(f)
            cmd = f" cat {f} | grep -v dev.off > {scriptFile.name} && Rscript {scriptFile.name} >& /dev/null  && ps2pdf {f.replace('.R','.ps')} {os.path.join(args.out_path,file_root)}.pdf && rm {f.replace('.R','.ps')} && rm {f} && rm *Rout"
            print(cmd)
            tmp_bash(cmd)
            
def plot_kinship(args):
    '''
    Plots the ditribution of kinship values.
    '''

    args.kinship_fig = os.path.join(args.out_path,args.prefix +'_kinship_distribution.pdf')
    if  os.path.isfile(args.kinship_fig):
        print(3,'kinship plot already done')
        return
    
    else:
       print('plotting kinship...')

    
    from verkko.binner import binner
   
    kin_dump = os.path.join(args.kinship_path, args.prefix +'_kinship.npy')
    if not os.path.isfile(kin_dump) or args.force: 
        print('uploading kinship data')
        kin_data = pd.read_csv(args.kin_file,sep = identify_separator(args.kin_file), usecols = ['Kinship']).values.flatten()
        print('dumping it')
        np.save(kin_dump,kin_data)
    else:
        kin_data = np.load(kin_dump)

    entries = len(kin_data)
    
    bin_data_file = os.path.join(args.kinship_path,'bin_average.npy')
    plot_data_file = os.path.join(args.kinship_path,'plot_data.npy')
    print('importing bin data..')
    if not os.path.isfile(bin_data_file) or not os.path.isfile(plot_data_file) or args.force:
        bins = binner.Bins(float,min(kin_data),max(kin_data),'log',1.05)
        print('bin data missing, generating...')
        countNotNormalized = bins.bin_count_divide(kin_data)
        count = np.array(binner.normalize(list(countNotNormalized)))
        binAvg = bins.bin_average(zip(kin_data,kin_data))
        binMask = ~np.ma.getmask(binAvg)
        plot_data = count[binMask]
        bin_data = binAvg[binMask]
        bin_data.dump(bin_data_file)
        plot_data.dump(plot_data_file)
    else:
        bin_data = np.load(bin_data_file,allow_pickle = True)
        plot_data = np.load(plot_data_file,allow_pickle = True)
        
    print('plotting...')

    import matplotlib as mpl
    mpl.use('Agg')
    from matplotlib import pyplot as plt
    import seaborn as sns
    sns.set(palette='Set2')
    import matplotlib.ticker as ticker
    fig = plt.figure()
    gs = mpl.gridspec.GridSpec(1,1)
    
    xPos = [0.0442,0.0884,0.177,0.354,0.51]
    texts = ['3rd Deg','2nd Deg','1st Deg','Twins/Duplicates']
    ax = fig.add_subplot(gs[0,0])
    ax2 = ax.twinx()
    ax2.set_ylim([0,1])
    cum_data = np.cumsum(plot_data)
    ax.plot(bin_data,plot_data, '--')
    ax2.plot(bin_data,cum_data, color = 'red', linewidth = 1)

    ax.set_xticks(xPos[:-1])
    ax.set_xlabel('kinship')
    ax.set_ylabel(r'P(k)')

    ax2.set_yticks([0,1])
    #ax2.set_yticklabels([0,entries])
    ax2.set_ylabel(r'$P(k < x) $')

    #add degree lines
    total = 0
    for i,pos in enumerate(xPos[:-1]):
        upper_mask = (kin_data >= pos)
        lower_mask = (kin_data < xPos[i+1])
        final_mask = upper_mask & lower_mask
        degree_data = str(len(kin_data[final_mask]))
        ax.axvline(x = pos, color = 'k',linestyle = '--',linewidth = 0.2)
        text_xpos = pos + (xPos[i+1] - pos)/2
        text = texts[i]
        ax.text(text_xpos,max(plot_data)*1.1,degree_data + '\n' +  text,fontsize = 5,horizontalalignment = 'center')

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_rotation(30)
        tick.label.set_fontsize(6)     

                 
    fig.savefig(args.kinship_fig)
    plt.close()
    print('done')

#######################
#------PEDIGREE-------#
#######################
    
    
def fix_fam(args):
    '''
    Adds sex info into a new fam file
    Sex code ('1' = male, '2' = female, '0' = unknown)
    '''
    args.new_fam =  args.kinship_bed.replace("kinship.bed",'pedigree.fam')
    pickle_path = os.path.join(args.out_path,'sex_dict.p')
    if not os.path.isfile(pickle_path) or args.force:
        sex_dict = dd(str)
        idx = [return_header(args.pheno_file).index(elem) for elem in ['FINNGENID','SEX']] # column indexes
        for fid,sex in basic_iterator(args.pheno_file,skiprows =1 ,columns = idx):
            sex_dict[fid] = '2' if sex == 'female' else '1'
        pickle.dump(sex_dict,open(pickle_path,'wb'))
    else:
        sex_dict = pickle.load(open(pickle_path,'rb'))

    print('generating new fam file ...')
    with open(args.new_fam,'wt') as o:
        for line in basic_iterator(args.kinship_bed.replace(".bed",'.fam')):
            new_sex = sex_dict[line[1]]
            if new_sex: #maybe iid is missing
                line[4] = str(new_sex)
            o.write('\t'.join(line) + '\n')
    print('done.')



    
def king_pedigree(args):
    """
    Return relatedness and build from king.
    """


    args.pedigree_log_file = os.path.join(args.out_path,args.prefix + '_pedigree.log')
    pedigree_root = os.path.join(args.pedigree_path, args.prefix +'_pedigree')
    pedigree_parents_file = pedigree_root + 'updateparents.txt'
    pedigree_ids_file = pedigree_root + 'updateids.txt'
    fix_fam(args)
    if not os.path.isfile(pedigree_parents_file) or mapcount(pedigree_parents_file) < 1 or args.force:
        args.force = True
        cmd= f'king -b {args.kinship_bed} --cpus {cpus}  --build --degree 3 --prefix {pedigree_root} --fam {args.new_fam} '
        print(cmd)
        with open(args.pedigree_log_file,'wt') as f: subprocess.call(shlex.split(cmd),stdout = f)

        # update fam file
        cmd = f"plink2 --fam {args.new_fam} --update-ids {pedigree_ids_file}  --make-just-fam --out {args.new_fam.replace('.fam','')}"
        print(cmd)
        subprocess.call(shlex.split(cmd))
        cmd = f"plink2  --fam {args.new_fam} --update-parents {pedigree_parents_file} --make-just-fam --out {args.new_fam.replace('.fam','')}"
        print(cmd)
        subprocess.call(shlex.split(cmd))
    
    else:
        print(f'pedigree files already generated')

       
    scriptFile = NamedTemporaryFile(delete=True)
    tmp_file = scriptFile.name


    args.log_file = os.path.join(args.out_path,args.prefix + '.log')
    with open(args.log_file,'wt') as o:
        o.write('#----KING SUMMARY-----#\n')

    tmp_bash(f"cat {args.kinship_log_file} | grep 'Relationship summary' -A 3| grep MZ  > {tmp_file}")
    tmp_bash(f"cat {args.kinship_log_file} | grep 'Relationship summary' -A 3| grep Inference  >> {tmp_file}")
    tmp_bash(f"head -n2 {tmp_file} | cut -f 2- >> {args.log_file}")
    
    with open(args.log_file,'at') as o:
        o.write('\n#----MANUAL COUNT-----#\n')
        desc_list = []        
                
        desc ='TRIOS: NUMBER OF FINNGEN_MOTHER FINNGEN_FATHER COUPLES WHO HAVE AT LEAST ONE SON IN FINNGEN'
        basic_cmd = f"""cat {pedigree_parents_file} |  awk '{{print $3"_"$4}}' | sort  """
        out_cmd = f" | wc -l >{tmp_file}"
        trio_cmd =  f""" {basic_cmd} |  uniq -c |  grep -o '\\bF\w*_FG\w*'  {out_cmd}""" 
        tmp_bash(trio_cmd)
        trios = int(open(tmp_file).read())
        o.write(':'.join(['Trios',str(trios)]) + '\n')
        desc_list.append(desc)
        
        desc = "ALL TRIOS: TOTAL NUMBER OF TRIOS (I.E. COUNTING MULTIPLES)"
        all_trio_cmd =  f""" {basic_cmd} |  uniq -c |  grep  '\\bFG\w*_FG\w*' |  awk '{{count+=$1}} END {{print count}}' > {tmp_file}""" 
        tmp_bash(all_trio_cmd)
        all_trios = int(open(tmp_file).read())
        o.write(':'.join(['All Trios',str(all_trios)]) + '\n')
        desc_list.append(desc)

        desc = "DUOS: NUMBER OF PARENT COUPLES WHERE ONLY ONE IS IN FINNGEN WHO HAVE A FINNGEN CHILD."
        duos_cmd =  f" {basic_cmd} |  uniq -c | grep -o '\bFG\w*\|\w*_FG\w*' | grep -v '\\bFG\w*_FG\w*'  {out_cmd} " 
        tmp_bash(duos_cmd)
        duos = int(open(tmp_file).read()) 
        o.write(':'.join(['Duos',str(duos)]) + '\n')
        desc_list.append(desc)
        
        desc = "ALL DUOS: TOTAL NUMBER OF DUOS COUNTING MULTIPLES"
        all_duos_cmd =  f" {basic_cmd} |  uniq -c | grep  '\bFG\w*\|\w*_FG\w*' | grep -v '\\bFG\w*_FG\w*'| awk '{{count+=$1}} END {{print count}}' > {tmp_file} " 
        tmp_bash(all_duos_cmd)
        all_duos = int(open(tmp_file).read()) 
        o.write(':'.join(['All Duos',str(all_duos)]) + '\n')
        desc_list.append(desc)
        
        desc = "SIBLINGS: NUMBER OF FAMILIES IN WHICH THERE ARE AT LEAST TWO FINNGEN SAMPLES WITH SAME PARENTS"
        sib_cmd = f"{basic_cmd} |  uniq -d  |  grep -v '0_' | grep -v '_0' {out_cmd}"
        tmp_bash(sib_cmd)
        sibs = int(open(tmp_file).read())
        o.write(':'.join(['Siblings',str(sibs)]) + '\n')
        desc_list.append(desc)
        
        desc = "TOTAL SIBLINGS : ALL SIBLINGS INCLUNDING MULTIPLE FROM EACH FAMILY"
        all_sib_cmd = f"{basic_cmd} |  uniq -cd  |  grep -v '0_' | grep -v '_0' | awk '{{count+=$1}} END {{print count}}' > {tmp_file}"
        tmp_bash(all_sib_cmd)
        all_sibs = int(open(tmp_file).read())
        o.write(':'.join(['All Siblings',str(sibs)]) + '\n')
        desc_list.append(desc)

        o.write('\n#----DESCRIPTION-----#\n')
        o.write('\n'.join(desc_list) + '\n')
        

def main(args):

    args.plink_path = os.path.join(args.out_path,'plink')
    make_sure_path_exists(args.plink_path)
    
    pretty_print("BUILD BED")
    build_bed(args)
      
    pretty_print("KINSHIP")
    args.kinship_path = os.path.join(args.out_path,'kinship')
    make_sure_path_exists(args.kinship_path)
    kinship(args)
    plot_kinship(args)
    
    pretty_print("PEDIGREE")
    args.pedigree_path = os.path.join(args.out_path,'pedigree')
    make_sure_path_exists(args.pedigree_path)
    king_pedigree(args)
 
    if args.release:
        # move data to doc
        import glob
        doc_path = os.path.join(args.out_path,'documentation')
        make_sure_path_exists(doc_path)
        for f in get_filepaths(doc_path): os.remove(f) # clean path else shutil.copy might fail
        for pdf in glob.glob(os.path.join(args.out_path,'*pdf')):
            shutil.copy(pdf,os.path.join(doc_path,os.path.basename(pdf)))
        for log in [args.log_file,args.kinship_log_file,args.pedigree_log_file]:
            shutil.copy(log,os.path.join(doc_path,os.path.basename(log)))
            
        data_path = os.path.join(args.out_path,'data')
        make_sure_path_exists(data_path)
        for f in get_filepaths(data_path): os.remove(f) # clean path else shutil.copy might fail
        all_files = get_filepaths(args.out_path)
        for f in [f for f in get_filepaths(args.out_path) if f.endswith(('kin0','con'))]:
            shutil.copy(f,os.path.join(data_path,os.path.basename(f)))
        endings = ('.bed','.fam','.bim','afreq',)
        for f in [f for f in get_filepaths(args.out_path) if f.endswith(endings) and ('kinship' in f or 'pedigree' in f)]:
            shutil.copy(f,os.path.join(data_path,os.path.basename(f)))
        
                 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="kinship analysis & pedigree")

    
    parser.add_argument("--extract", type=file_exists, help =  "Path to list of variants to include",default = False)
    parser.add_argument("--fam", type=file_exists, help =  "fam file to filter down")
    parser.add_argument("-f", '--pheno-file', metavar='F', type=file_exists, help="Phenotype filepath",required = True)
    parser.add_argument("-b", '--bed', metavar='F', type=file_exists, help="BED filepath", required=True)
    parser.add_argument('-o', "--out-path", type=str, help="folder in which to save the results", required=True)
    parser.add_argument('--prefix',  type=str, help="Output prefix", required=True)
    parser.add_argument('--force',help='Flag on whether to force run',action = "store_true")
    parser.add_argument('--release',help='Flag on whether to force run',action = "store_true")
    
    args = parser.parse_args()
    make_sure_path_exists(args.out_path)

    main(args)
    



