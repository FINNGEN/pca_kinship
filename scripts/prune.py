import os,subprocess,shlex,argparse,multiprocessing,shutil
from utils import get_path_info,file_exists,make_sure_path_exists,mapcount,tmp_bash,is_gz_file,get_filepaths
cpus = multiprocessing.cpu_count()
import numpy as np
from pathlib import Path


def filter_variants(args):
    """
    Returns hq variants based on min(INFO_SCORE) > VALUE
    """
    
    info_file,info_filter = args.info
    filtered_snps= os.path.join(args.out_path,info_filter + '_hq_snps.txt')
    if not os.path.isfile(filtered_snps) or mapcount(filtered_snps) <1 or args.force:
        args.force = True
        cat_cmd = "zcat " if is_gz_file(info_file) else "cat "
        cmd = cat_cmd + info_file +  " | awk 'BEGIN{min=9}{for(i=1;i<=NF;i++){min=(min<$i)?min:$i}print $1,min;min=9}'   | sed -E 1d | awk '{ if($2 >= " + info_filter +") { print $1}}' > " + filtered_snps
        print(cmd)
        tmp_bash(cmd)
    else:
        print(f"{filtered_snps} already generated")

    print(f"{mapcount(filtered_snps)} hq snps")
    return filtered_snps

def ld_pruning(args):
    """
    Iteratively ld prunes until target snpscount is reached.
    """     
    pruned_variants = f"{args.out_root}.{args.count}.prune.in"
    if not os.path.isfile(pruned_variants) or args.force:
        cmd = f"plink2  --out {args.out_root}.{args.count} --read-freq  {args.bed.replace('.bed','.afreq')} --threads {cpus} --bfile {args.bed.replace('.bed','')} --extract {args.snplist} --indep-pairwise {' '.join(map(str,args.ld))} {args.pargs} "
        print(cmd)
        subprocess.call(shlex.split(cmd))
    print(f"Variants left after pruning step {args.count}: {mapcount(pruned_variants)}")
    
    args.snplist = pruned_variants
    args.ld[-1] = round(args.ld[-1] - args.step,trailing_zeros(args.step))
    # keep pruning until we either go below the threshold or we run ot ouf r2 steps
    while args.target < mapcount(args.snplist) and float(args.ld[-1]) > 0 :
        args.count +=1
        ld_pruning(args)
             
def trailing_zeros(x):
    """
    Returns order of magnnitude of float number.
    """
    s = str(x).split('.')[-1]
    return len(s) - len(s.strip('0')) + 1

def main(args):

    if args.info:
        args.snplist = filter_variants(args)
    elif args.extract:
        args.snplist = args.extract
    else:
        args.snplist = args.bed.replace('.bed','.bim')

    args.initial_variants = mapcount(args.snplist)
    print(f"{mapcount(args.snplist)} initial variants.")   
    
    args.count= args.success = 0
    ld_pruning(args)

    args.pruned_variants = f"{args.out_root}.prune.in"
  
    # get closest value to target_snps and return random target_snps variants
    counts = {f"{args.out_root}.{i}.prune.in" : mapcount(f"{args.out_root}.{i}.prune.in") for i in range(args.count+1)}
    best_prune = [key for key in counts if counts[key] == min(counts.values(), key=lambda x:abs(x-args.target))][0]

    # copy plink log as final log
    args.log_file =f"{args.out_root}.prune.log" 
    cmd = f"cp {best_prune.replace('.prune.in','.log')} {args.log_file} "
    print(cmd)
    subprocess.call(shlex.split(cmd))
    
    print(f'closest pruning was {best_prune} with {counts[best_prune]} snps')
    cmd = f"cat {best_prune} | shuf | head -n {args.target} | sort > {args.pruned_variants} "
    tmp_bash(cmd)             

    args.final_variants = mapcount(args.pruned_variants)
    print(f'Final variants: {args.final_variants}')

    # return final value of ld params
    with open(args.log_file) as f:
        for line in f :
            if "--indep-pairwise" in line:
                args.final_ld = line.strip()
                break

    release(args)
    
def release(args):

    import glob
    doc_path = os.path.join(args.out_path,'documentation')
    data_path = os.path.join(args.out_path,'data')
    for path in [data_path,doc_path]:
        make_sure_path_exists(path)
        for f in get_filepaths(path): os.remove(f) # clean path else shutil.copy might fail

    for f in [args.pruned_variants]:
        shutil.copy(f,os.path.join(data_path,os.path.basename(f)))

    for f in [args.log_file]:
        shutil.copy(f,os.path.join(doc_path,os.path.basename(f)))


    # README
    if args.extract: variant_filter = f"{args.initial_variants} starting variants from {args.extract}"
    if args.info: variant_filter = f"{args.initial_variants} starting variants: only variants with a minimum info score of {args.info[1]} in all batches are kept."
    else: variant_filter = f"{args.initial_variants} starting variants."
    readme = os.path.join(args.data_path,'prune.README') 
    with open(os.path.join(args.out_path,args.prefix + '_prune_readme'),'wt') as o, open(readme,'rt') as i:
        with open(args.log_file) as tmp: summary = tmp.read()
        word_map = {
            '[PREFIX]':args.prefix,
            '[TARGET]':args.target,
            '[INITIAL_LD]':args.initial_ld,
            '[SNPS]':args.final_variants,
            '[FINAL_LD]':args.final_ld,
            '[PARGS]':args.pargs,
            '[STEP]':args.step,
            '[FILTER]':variant_filter
    }
        for line in i:
            for kw in word_map:
                if kw in line:
                    line = line.replace(kw,str(word_map[kw]))
            o.write(line)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LD pruning of a bed file")

     # mutually exclusive group. Either calculate the info score variants or pass directly pass the variants.
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--info',nargs =2, metavar = ('FILE','VALUE'), help = 'Info score file and filter value')
    group.add_argument("--extract", type=file_exists, help =  "Path to list of variants to include")

    # required args
    parser.add_argument('--bed',  type=file_exists, help="BED filepath", required=True)
    parser.add_argument("-o","--out-path", type=str, help="folder in which to save the results", required=True)
    parser.add_argument('--prefix',  type=str, help="Output prefix", required=True)

    # optional args
    parser.add_argument('--pargs',type = str,help='extra plink args',default = ' --snps-only --chr 1-22 --max-alleles 2 --maf 0.01  ')
    parser.add_argument('--force',help='Flag on whether to force run',action = "store_true")
    parser.add_argument('--target',type = int,help = "Target number of snps after pruning",default = 80000)
    parser.add_argument('--ld',nargs=4,type = float,metavar = ('SIZE','STEP','THRESHOLD','STEP2'),help ='size,step,threshold,threshold_step',default = [50,5,0.9,0.05])
    parser.add_argument('--release',action = 'store_true',help = 'Flag for data release',default = False)

    args = parser.parse_args()
    # PREPROCESSING
    *args.ld,args.step = args.ld
    args.initial_ld = list(args.ld)
    
    make_sure_path_exists(args.out_path)
    args.out_root = os.path.join(args.out_path,args.prefix)

    args.parent_path = Path(os.path.realpath(__file__)).parent.parent
    args.data_path = os.path.join(args.parent_path,'data/')
    
    main(args)

