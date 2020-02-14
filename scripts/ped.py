import os,pickle, subprocess,shlex,argparse,multiprocessing
from utils import basic_iterator,return_header,mapcount,get_path_info,file_exists,make_sure_path_exists,cpus,tmp_bash,pretty_print,identify_separator,NamedTemporaryFile,get_filepaths
from collections import defaultdict as dd
import numpy as np
cpus = multiprocessing.cpu_count()

######################
#---BUILD BED FILE---#
######################
def build_bed(args):
    """ 
    Builds bed with hq variants for kinship. This is the core data set used for the all kinship analysis.    
    """

 
    args.kinship_bed =os.path.join(args.plink_path,args.prefix + '.kinship.bed')
    if not os.path.isfile(args.kinship_bed) or not mapcount(args.kinship_bed) or args.force:
        args.force = True
        # filter for only 30k samples 
        keep = ""
        if args.test:
            fam_file = args.bed.replace('.bed','.fam')
            samples = np.loadtxt(fam_file,dtype = str, delimiter = identify_separator(fam_file),usecols = [0,1])
            samples = samples[np.random.choice(samples.shape[0], 30000, replace=False), :]
            scriptFile = NamedTemporaryFile(delete=True)
            np.savetxt(scriptFile.name,samples,delimiter = '\t',fmt = '%s')
            keep = f"--keep {scriptFile.name}"

        extract = f"--extract {args.extract}" if hasattr(args,"extract") else ""
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

    log_file = os.path.join(args.out_path,args.prefix + '.kinship.log')
    with open(log_file,'wt') as o: o.write('KINSHIP\n')
    
    related_file = os.path.join(args.kinship_path,f"{args.prefix}.kin0")
    # RETURN RELATED AND PLOT FAMILIES
    if not os.path.isfile(related_file) or mapcount(related_file) < 1 or args.force:
        args.force = True
        cmd = f'king --cpus {cpus} -b {args.kinship_bed} --related --duplicate --degree 2 --prefix {os.path.join(args.kinship_path,args.prefix)} --rplot '
        print(cmd)
        with open(log_file,'at') as f: subprocess.call(shlex.split(cmd),stdout = f)
        # produce R scripts, fix them and run them
    else:
        print("related file already generated")
         
    unrelated_file =    os.path.join(args.kinship_path,args.prefix ) + 'unrelated.txt'
    # RETURN UNRELATED AND DUPLICATES
    if not os.path.isfile(unrelated_file) or args.force:
        args.force = True
        cmd =f"king -b {args.kinship_bed}  --unrelated  --degree 2 --cpus {cpus} --prefix {os.path.join(args.kinship_path,args.prefix)} "
 #       print(cmd)
#        with open(log_file,'at') as f: subprocess.call(shlex.split(cmd),stdout = f)
    else:
        print("unrelated already generated")

        
    if args.force:
        scriptFile = NamedTemporaryFile(delete=True)
        for f in [f for f in get_filepaths(args.kinship_path) if f.endswith('.R')]:
            file_path,file_root,file_extension = get_path_info(f)
            cmd = f" cat {f} | grep -v dev.off > {scriptFile.name} && Rscript {scriptFile.name} >& /dev/null  && ps2pdf {f.replace('.R','.ps')} {os.path.join(args.out_path,file_root)}.pdf && rm {f.replace('.R','.ps')} && rm {f}"
            print(cmd)
            tmp_bash(cmd)

    args.log_file = os.path.join(args.out_path,args.prefix + '.log')
    tmp_bash(f"cat {log_file} | grep 'Relationship summary' -A 3 > {args.log_file}")


#######################
#------PEDIGREE-------#
#######################
    
    
def fix_fam(args):
    '''
    Adds sex info into a new fam file
    Sex code ('1' = male, '2' = female, '0' = unknown)
    '''

    args.new_fam =  args.kinship_bed.replace(".bed",'.updated.fam')
    print('generating new fam file ...')
    pickle_path = os.path.join(args.out_path,'sex_dict.p')
    if not os.path.isfile(pickle_path) or args.force:
        sex_dict = dd(str)
        idx = [return_header(args.pheno_file).index(elem) for elem in ['FINNGENID','SEX']] # column indexes
        for fid,sex in basic_iterator(args.pheno_file,skiprows =1 ,columns = idx):
            sex_dict[fid] = '2' if sex == 'female' else '1'
        pickle.dump(sex_dict,open(pickle_path,'wb'))
    else:
        sex_dict = pickle.load(open(pickle_path,'rb'))

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


    log_file = os.path.join(args.out_path,args.prefix + '.pedigree.log')

    pedigree_root = os.path.join(args.pedigree_path, args.prefix +'.pedigree')
    pedigree_parents_file = pedigree_root + 'updateparents.txt'
    pedigree_ids_file = pedigree_root + 'updateids.txt'
    fix_fam(args)
    if not os.path.isfile(pedigree_parents_file) or mapcount(pedigree_parents_file) < 1 or args.force:
        args.force = True
        cmd= f'king -b {args.kinship_bed} --cpus {cpus}  --build --degree 3 --prefix {pedigree_root} --fam {args.new_fam} '
        print(cmd)
        with open(log_file,'wt') as f: subprocess.call(shlex.split(cmd),stdout = f)
      
    else:
        print(f'pedigree files already generated')


    
    # update fam file
    cmd = f"plink2 --fam {args.new_fam} --update-ids {pedigree_ids_file}  --make-just-fam --out {args.new_fam.replace('.fam','')}"
    print(cmd)
    subprocess.call(shlex.split(cmd))
    cmd = f"plink2  --fam {args.new_fam} --update-parents {pedigree_parents_file} --make-just-fam --out {args.new_fam.replace('.fam','')}"
    print(cmd)
    subprocess.call(shlex.split(cmd))

    
    scriptFile = NamedTemporaryFile(delete=True)
    tmp_file = scriptFile.name

    # NUMBER OF FINNGEN_MOTHER FINNGEN_FATHER COUPLES WHO HAVE AT LEAST ONE SON IN FINNGEN
    basic_cmd = f"""cat {pedigree_parents_file} |  awk '{{print $3"_"$4}}' | sort  """
    out_cmd = f" | wc -l >{tmp_file}"
    trio_cmd =  f""" {basic_cmd} |  uniq -c |  grep -o '\\bF\w*_FG\w*'  {out_cmd}""" 
    tmp_bash(trio_cmd)
    print(trio_cmd)
    trios = int(open(tmp_file).read())

    # TOTAL NUMBER OF TRIOS (I.E. COUNTING MULTIPLES)
    all_trio_cmd =  f""" {basic_cmd} |  uniq -c |  grep  '\\bFG\w*_FG\w*' |  awk '{{count+=$1}} END {{print count}}' > {tmp_file}""" 
    tmp_bash(all_trio_cmd)
    print(all_trio_cmd)
    all_trios = int(open(tmp_file).read())

    # NUMBER OF PARENT COUPLES WHERE ONLY ONE IS IN FINNGEN WHO HAVE A FINNGEN CHILD.
    duos_cmd =  f" {basic_cmd} |  uniq -c | grep -o '\bFG\w*\|\w*_FG\w*' | grep -v '\\bFG\w*_FG\w*'  {out_cmd} " 
    tmp_bash(duos_cmd)
    print(duos_cmd)
    duos = int(open(tmp_file).read()) 

    # TOTAL NUMBER OF DUOS COUNTING MULTIPLES
    all_duos_cmd =  f" {basic_cmd} |  uniq -c | grep  '\bFG\w*\|\w*_FG\w*' | grep -v '\\bFG\w*_FG\w*'| awk '{{count+=$1}} END {{print count}}' > {tmp_file} " 
    tmp_bash(all_duos_cmd)
    print(all_duos_cmd)
    all_duos = int(open(tmp_file).read()) 

    
    # SIBLINGS: NUMBER OF FAMILIES IN WHICH THERE ARE AT LEAST TWO FINNGEN SAMPLES WITH SAME PARENTS
    sib_cmd = f"{basic_cmd} |  uniq -d  |  grep -v '0_' | grep -v '_0' {out_cmd}"
    tmp_bash(sib_cmd)
    sibs = int(open(tmp_file).read())

    # TOTAL SIBLINGS
    all_sib_cmd = f"{basic_cmd} |  uniq -cd  |  grep -v '0_' | grep -v '_0' | awk '{{count+=$1}} END {{print count}}' > {tmp_file}"
    tmp_bash(all_sib_cmd)
    print(all_sib_cmd)
    all_sibs = int(open(tmp_file).read())
    
    a = zip(['Siblings','All Siblings','Duos','All Duos','Trios','All Trios'],[sibs,all_sibs,duos,all_duos,trios,all_trios])
    with open(args.log_file,'at') as o:
        o.write('\n\n#----MANUAL COUNT-----#\n')
        for elem in a:
            o.write(':'.join(map(str,elem)) + '\n')
            
    
def main(args):

    args.plink_path = os.path.join(args.out_path,'plink')
    make_sure_path_exists(args.plink_path)
    
    pretty_print("BUILD BED")
    build_bed(args)
      
    pretty_print("KINSHIP")
    args.kinship_path = os.path.join(args.out_path,'kinship')
    make_sure_path_exists(args.kinship_path)
    kinship(args)

    pretty_print("PEDIGREE")
    args.pedigree_path = os.path.join(args.out_path,'pedigree')
    make_sure_path_exists(args.pedigree_path)
    king_pedigree(args)
 
 
    
                 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="kinship analysis & pedigree")

    
    parser.add_argument("--extract", type=file_exists, help =  "Path to list of variants to include")
    parser.add_argument("-f", '--pheno-file', metavar='F', type=file_exists, help="Phenotype filepath",required = True)
    parser.add_argument("-b", '--bed', metavar='F', type=file_exists, help="BED filepath", required=True)
    parser.add_argument('-o', "--out-path", type=str, help="folder in which to save the results", required=True)
    parser.add_argument('--prefix',  type=str, help="Output prefix", required=True)
    parser.add_argument('--test',action = 'store_true',help = 'Flag for quick run with less samples')
    parser.add_argument('--force',help='Flag on whether to force run',action = "store_true")

    args = parser.parse_args()
    make_sure_path_exists(args.out_path)

    main(args)
    



