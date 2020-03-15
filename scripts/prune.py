import os,subprocess,shlex,argparse,multiprocessing
from utils import get_path_info,file_exists,make_sure_path_exists,mapcount,tmp_bash,is_gz_file
cpus = multiprocessing.cpu_count()
import numpy as np

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

    args.count= args.success = 0
    ld_pruning(args)

    pruned_variants = f"{args.out_root}.prune.in"
  
    # get closest value to target_snps and return random target_snps variants
    counts = {f"{args.out_root}.{i}.prune.in" : mapcount(f"{args.out_root}.{i}.prune.in") for i in range(args.count+1)}
    best_prune = [key for key in counts if counts[key] == min(counts.values(), key=lambda x:abs(x-args.target))][0]
    
    print(f'closest pruning was {best_prune} with {counts[best_prune]} snps')
    cmd = f"cat {best_prune} | shuf | head -n {args.target} | sort > {pruned_variants} "
    tmp_bash(cmd)             
            
    print(f'Final variants: {mapcount(pruned_variants)}')
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LD pruning of a bed file")

     # mutually exclusive group. Either calculate the info score variants or pass directly pass the variants.
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--info',nargs =2, metavar = ('FILE','VALUE'), help = 'Info score file and filter value')
    group.add_argument("--extract", type=file_exists, help =  "Path to list of variants to include")

    # required args
    parser.add_argument('--bed',  type=file_exists, help="BED filepath", required=True)
    parser.add_argument("--out-path", type=str, help="folder in which to save the results", required=True)
    parser.add_argument('--prefix',  type=str, help="Output prefix", required=True)

    # optional args
    parser.add_argument('--pargs',type = str,help='extra plink args',default = ' --maf ')
    parser.add_argument('--force',help='Flag on whether to force run',action = "store_true")
    parser.add_argument('--target',type = int,help = "Target number of snps after pruning",default = 80000)
    parser.add_argument('--ld',nargs=4,type = float,metavar = ('SIZE','STEP','THRESHOLD','STEP2'),help ='size,step,threshold,threshold_step',default = [50,5,0.9,0.05])
    
    args = parser.parse_args()

    # PREPROCESSING
    *args.ld,args.step = args.ld
    make_sure_path_exists(args.out_path)
    args.out_root = os.path.join(args.out_path,args.prefix)
    main(args)
