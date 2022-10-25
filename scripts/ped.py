import os,pickle, subprocess,shlex,argparse,shutil,glob,logging
from pathlib import Path
from utils import basic_iterator,return_header,mapcount,get_path_info,file_exists,make_sure_path_exists,cpus,tmp_bash,pretty_print,NamedTemporaryFile,get_filepaths,read_int,log_levels,print_msg_box
from collections import defaultdict as dd
from collections import Counter
import numpy as np
from pca_scripts import kinship_plots as kp

degree_dict = {'Dup/MZ':0,'PO':1,'FS':1,'2nd':2,'3rd':3,'4th':4}

######################
#---BUILD BED FILE---#
######################

def build_bed(args,name='kinship',kwargs = ""):
    """ 
    Builds bed with hq variants for kinship. This is the core data set used for the all kinship analysis.    
    """
    args.kinship_bed =os.path.join(args.plink_path,args.prefix + f"_{name}.bed")
    args.maj_bim = os.path.join(args.plink_path,args.prefix + f"_{name}_maj.bim")
    if not os.path.isfile(args.kinship_bed) or args.force:
        args.force = True
        keep = f"--keep {args.fam}" if args.fam and mapcount(args.fam) > 100 else ""
        extract = f"--extract {args.extract}" if args.extract else ""
        cmd = f"plink2 --bfile {args.bed.replace('.bed','')} {extract} --threads {cpus}  --make-bed --out {args.kinship_bed.replace('.bed','')} {keep} {kwargs}"
        logging.debug(cmd)
        subprocess.call(shlex.split(cmd))    

    print('done.')



    
######################
#------KINSHIP-------#
######################

def kinship(args):
    """
    Returns degree 3 kinship data.
    """

    args.kinship_log_file = os.path.join(args.out_path,args.prefix + '_kinship.log')  
    args.kin_file = os.path.join(args.kinship_path,f"{args.prefix}.kin0")
    args.dup_file = os.path.join(args.kinship_path,f"{args.prefix}.con")
    args.all_segs = os.path.join(args.kinship_path,f"{args.prefix}allsegs.txt")
    # RETURN RELATED AND PLOT FAMILIES
    if not os.path.isfile(args.kin_file) or mapcount(args.kin_file) < 1 or args.force:
        args.force = True
        cmd = f'king --cpus {cpus} -b {args.kinship_bed} --related --duplicate --degree 3 --prefix {os.path.join(args.kinship_path,args.prefix)} --rplot '
        logging.debug(cmd)
        with open(args.kinship_log_file,'wt') as f: subprocess.call(shlex.split(cmd),stdout = f)
        # filter un/4th
        kin_filter = args.kin_file.replace('kin0','tmp')
        cmd = f"cat {args.kin_file} | grep -vw 'UN'  | grep -vw '4th' > {kin_filter} && mv {kin_filter} {args.kin_file} && rm {kin_filter}"
        logging.debug(cmd)
        tmp_bash(cmd)

    else:
        print("related file already generated")

    # R SCRIPTS
    if args.force:
        scriptFile = NamedTemporaryFile(delete=True)
        for f in [f for f in get_filepaths(args.kinship_path) if f.endswith('.R')]:
            file_path,file_root,file_extension = get_path_info(f)
            cmd = f" cat {f} | grep -v dev.off > {scriptFile.name} && Rscript {scriptFile.name} >& /dev/null  && ps2pdf {f.replace('.R','.ps')} {os.path.join(args.out_path,file_root)}.pdf && rm {f.replace('.R','.ps')} && rm {f} && rm *Rout"
            logging.debug(cmd)
            tmp_bash(cmd)

            
def degree_summary(args):
    """
    Creates a summary of lowest degrees
    """
    args.degree_table = os.path.join(args.kinship_path,f"{args.prefix}_degree_summary.txt")
    if not os.path.isfile(args.degree_table) or args.force:
        header = return_header(args.kin_file)
        idx = [header.index(elem) for elem in ['ID1','ID2','InfType']]
        print(idx)

        degree_d = dd(lambda : np.inf)
        for key,val in degree_dict.items():degree_d[key] = val
        
        sample_deg_dict = dd(lambda:np.inf)
        deg_iterator = basic_iterator(args.kin_file,columns = idx,skiprows=1)
        for id1,id2,inf in deg_iterator:
            # match inftype to numerical degree
            deg = degree_dict[inf] 
            if deg < sample_deg_dict[id1] : sample_deg_dict[id1] = deg
            if deg < sample_deg_dict[id2] : sample_deg_dict[id2] = deg       

        with open(args.degree_table,'wt') as o:
            o.write('|Degree | Sample Count' + '|\n')
            o.write('|--|--|' + '\n')
            count = Counter(sample_deg_dict.values())
            for deg in sorted(count.keys()):
                if deg < np.inf:
                    o.write('|' + '|'.join(map(str,(deg,count[deg]))) + '|\n')         

    else:
        print(f"degree summary already generated")

    print('done.')


#######################
#------PEDIGREE-------#
#######################
    


def fix_fam(args):
    '''
    Adds sex info into a new fam file
    Sex code ('1' = male, '2' = female, '0' = unknown)
    '''

    sex_dict = dd(str)
    idx = [return_header(args.pheno_file).index(elem) for elem in ['FINNGENID','SEX']] # column indexes
    for fid,sex in basic_iterator(args.pheno_file,skiprows =1 ,columns = idx):
        sex_dict[fid] = '2' if sex == 'female' else '1'
    args.new_fam =  os.path.join(args.out_path,args.prefix + '_pedigree.fam')
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
    args.segs = pedigree_root + ".segments.gz"
    pedigree_parents_file = pedigree_root + 'updateparents.txt'
    pedigree_ids_file = pedigree_root + 'updateids.txt'

    fix_fam(args)

    king_pedigree_log = args.pedigree_log_file.replace('.log','.king.log')
    if not os.path.isfile(pedigree_parents_file)  or not os.path.isfile(args.segs) or args.force:
        args.force = True
        cmd= f'king -b {args.kinship_bed} --cpus {cpus}  --build --degree 3 --prefix {pedigree_root} --fam {args.new_fam}  2>&1 | tee {king_pedigree_log}'
        tmp_bash(cmd,True)
    else:
        print(f'pedigree files already generated')

    #update fam file and create logs
    tmp_log = args.pedigree_log_file.replace('.log','.plink.log')
    with open(tmp_log,'wt') as f:
        cmd = f"plink2 --fam {args.new_fam} --update-ids {pedigree_ids_file}  --make-just-fam --out {args.new_fam.replace('.fam','')}"
        logging.debug(cmd)
        subprocess.call(shlex.split(cmd),stdout = f,stderr =f)    
        cmd = f"plink2  --fam {args.new_fam} --update-parents {pedigree_parents_file} --make-just-fam --out {args.new_fam.replace('.fam','')}"
        logging.debug(cmd)
        subprocess.call(shlex.split(cmd),stdout = f,stderr =f)    

        cmd = f"cat {king_pedigree_log} > {args.pedigree_log_file} &&  cat {tmp_log} >> {args.pedigree_log_file} && rm {tmp_log}"
        tmp_bash(cmd,False)
        
    args.newparents = mapcount(pedigree_parents_file)
    args.newfids = mapcount(pedigree_ids_file)
    args.pedigree_parents_file = pedigree_parents_file

    print('done.')
     

     
def release_log(args):
    """
    Logs a bunch of stuff for the README.
    """
    scriptFile = NamedTemporaryFile(delete=True)
    tmp_file = scriptFile.name
    args.log_file = os.path.join(args.out_path,args.prefix + '.log')


    with open(args.log_file,'wt') as o:
        o.write('\n### Manual Count \n')
        
        # NUMBER OF COUPLES PER KINSHIP TYPE
        idx = return_header(args.kin_file).index('InfType')
        data = np.loadtxt(args.kin_file,usecols=idx,dtype =str)
        count =Counter(data)

        o.write('\n|Kinship Type|Number of couples|\n')
        o.write('|--|--|\n')
        for key in degree_dict:
            c = count[key]if key in count else 0
            o.write('|' + '|'.join([key,str(c)]) + '|\n')           

        
        # DUOS/TRIOS ETC
        o.write('\n|Family Structure Type|Count|Description|\n')
        o.write('|--|--|--|\n')

        desc_list = []
        desc ='Number of Finngen_mother Finngen_father couples who have at least one child in Finngen'
        basic_cmd = f"""cat {args.new_fam} |  awk '{{print $3"_"$4}}' | sort  """
        out_cmd = f" | wc -l >{tmp_file}"
        trio_cmd =  f""" {basic_cmd} |  uniq -c |  grep -o '\\bF\w*_FG\w*'  {out_cmd}""" 
        tmp_bash(trio_cmd)
        trios =  read_int(tmp_file)
        o.write('|' + '|'.join(['Trios',str(trios),desc]) + '|\n')
        
        desc = "Total number of trios (i.e. counting multiples)"
        all_trio_cmd =  f""" {basic_cmd} |  uniq -c |  grep  '\\bFG\w*_FG\w*' |  awk '{{count+=$1}} END {{print count}}' > {tmp_file}""" 
        tmp_bash(all_trio_cmd)
        all_trios = read_int(tmp_file)
        o.write('|' + '|'.join(['All Trios',str(all_trios),desc]) + '|\n')

        desc = "Parent - child duos where the other parent is not in Finngen"
        duos_cmd =  f" {basic_cmd} |  uniq -c | grep -o '\\bFG\w*\|\w*_FG\w*' | grep -v '\\bFG\w*_FG\w*'  {out_cmd} "
        tmp_bash(duos_cmd)
        duos =  read_int(tmp_file) 
        o.write('|' + '|'.join(['Duos',str(duos),desc]) + '|\n')
        
        desc = "Total number of duos counting multiples"
        all_duos_cmd =  f" {basic_cmd} |  uniq -c | grep  '\\bFG\w*\|\w*_FG\w*' | grep -v '\\bFG\w*_FG\w*'| awk '{{count+=$1}} END {{print count}}' > {tmp_file} "
        tmp_bash(all_duos_cmd)
        all_duos = read_int(tmp_file)  
        o.write('|' + '|'.join(['All Duos',str(all_duos),desc]) + '|\n')
        
        desc = "Number of families in which there are at least two finngen samples with same parents"
        sib_cmd = f"{basic_cmd} |  uniq -d  |  grep -v '0_' | grep -v '_0' {out_cmd}"
        tmp_bash(sib_cmd)
        sibs = read_int(tmp_file)  
        o.write('|' + '|'.join(['Siblings',str(sibs),desc]) + '|\n')

        desc = "Total number of siblings including multiple from each family"
        all_sib_cmd = f"{basic_cmd} |  uniq -cd  |  grep -v '0_' | grep -v '_0' | awk '{{count+=$1}} END {{print count}}' > {tmp_file}"
        tmp_bash(all_sib_cmd)
        all_sibs = read_int(tmp_file)
        o.write('|' + '|'.join(['All Siblings',str(all_sibs),desc]) + '|\n')    


        
        
def release(args):

    doc_path = os.path.join(args.out_path,'documentation')
    data_path = os.path.join(args.out_path,'data')
    for path in [doc_path,data_path]:
        make_sure_path_exists(path)
        for f in get_filepaths(path): os.remove(f) # clean path else shutil.copy might fail

    # DOC
    for pdf in glob.glob(os.path.join(args.out_path,'*pdf')):# copy figures
        shutil.copy(pdf,os.path.join(doc_path,os.path.basename(pdf)))
        
    for log in [args.log_file,args.kinship_log_file,args.pedigree_log_file,args.degree_table,args.all_segs]:#copy log files
        shutil.copy(log,os.path.join(doc_path,os.path.basename(log)))

    # DATA
    for f in [args.kin_file,args.dup_file,args.segs]: # copy KING outputs
        shutil.copy(f,os.path.join(data_path,os.path.basename(f)))       

    for ending in ('.bed','.fam','.bim','.afreq') : #copy plink files
        f = args.kinship_bed.replace('.bed',ending)
        shutil.copy(f,os.path.join(data_path,os.path.basename(f)))
    shutil.copy(args.new_fam,os.path.join(data_path,os.path.basename(args.new_fam)))

    # README FILE
    parent_path = Path(os.path.realpath(__file__)).parent.parent
    readme = os.path.join(parent_path,'data','kinship.README')    
    with open( os.path.join(args.out_path,args.prefix + '_kinship_readme'),'wt') as o, open(readme,'rt') as i:
        with open(args.log_file) as tmp:summary = tmp.read()
        with open(args.degree_table) as tmp:degree_summary = tmp.read()  
        word_map = {
            '[PREFIX]':args.prefix,
            '[NEW_PARENTS]':args.newparents,
            '[NEW_FID]':args.newfids,
            '[RELATED_COUPLES]':mapcount(args.kin_file) -1,
            '[DUPLICATES]':mapcount(args.kin_file.replace('kin0','con')) -1,
            '[SUMMARY]': summary,
            '[DEGREE_SUMMARY]':degree_summary,
            '[N_SNPS]': mapcount(args.kinship_bed.replace('.bed','.bim')),
            '[N_SAMPLES]': mapcount(args.kinship_bed.replace('.bed','.fam'))
        }
        for line in i:
            for kw in word_map:
                if kw in line:
                    line = line.replace(kw,str(word_map[kw]))
            o.write(line)
            
def main(args):    
  
    args.plink_path = os.path.join(args.out_path,'plink')
    make_sure_path_exists(args.plink_path)

    args.sample_info = args.meta
    args.test = True
    
    pretty_print("BUILD BED")
    build_bed(args,kwargs = '--maj-ref', name = 'maj_ref')

    pretty_print("KINSHIP")
    args.kinship_path = os.path.join(args.out_path,'kinship')
    make_sure_path_exists(args.kinship_path)
    kinship(args)
    degree_summary(args)

    pretty_print("PEDIGREE")
    args.pedigree_path = os.path.join(args.out_path,'pedigree')
    make_sure_path_exists(args.pedigree_path)
    king_pedigree(args)

    pretty_print("PLOTS")
    if args.plot:
        kp.plot_kinship(args)
        kp.plot_degree_dist(args)
        kp.plot_batch_data(args)
    
    if args.release:
        pretty_print("BUILD BED PLINK2")
        build_bed(args,kwargs = '--freq', name = 'kinship')
        release_log(args)
        release(args)

    return True
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="kinship analysis & pedigree")
    
    parser.add_argument("--extract", type=file_exists, help =  "Path to list of variants to include",default = False)
    parser.add_argument("--fam", type=file_exists, help =  "Optional .fam file to subset individuals")
    parser.add_argument('--pheno-file', type=file_exists, help="Phenotype filepath. Needs to contain SEX column",required = True)
    parser.add_argument('--meta', type=file_exists, help="File with batch info fo samples")

    parser.add_argument('--bed', type=file_exists, help="BED filepath", required=True)
    parser.add_argument("-o","--out-path", type=str, help="folder in which to save the results", required=True)
    parser.add_argument('--prefix',  type=str, help="Output prefix", required=True)
    parser.add_argument('--force',help='Flag on whether to force run',action = "store_true")
    parser.add_argument('--release',help='Flag to structure output for release',action = "store_true")
    parser.add_argument('--plot',help='Flag to plot',action = "store_true")
    # LOGGING ET AL.
    parser.add_argument( "-log",  "--log",  default="warning", choices = log_levels, help=(  "Provide logging level. " "Example --log debug', default='warning'"))

    args = parser.parse_args()

    # logging level
    level = log_levels[args.log]
    logging.basicConfig(level=level,format="%(levelname)s: %(message)s")
    args.logging = logging

    if args.release: args.plot = True
    if args.plot:
        if not args.meta:
            raise ValueError('batch data required')
        
    make_sure_path_exists(args.out_path)
    args.misc_path = os.path.join(args.out_path,'misc')
    make_sure_path_exists(args.misc_path)

    success =main(args)
    if success:
        args.release= False
        args.logging.getLogger().setLevel(logging.WARNING)
        print_msg_box("\n~ SUMMARY ~\n",indent =30)
        main(args)
    
    
