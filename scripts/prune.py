import os,subprocess,shlex,argparse,multiprocessing
from utils import get_path_info,file_exists,make_sure_path_exists,mapcount,tmp_bash,is_gz_file
cpus = multiprocessing.cpu_count()

def filter_variants(args):
    """
    Returns hq variants based on min(INFO_SCORE) > VALUE
    """
    
    info_file,info_filter = args.info
    filtered_snps= os.path.join(args.variants_path,'hq_snps.txt')
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

    
    pruned_variants = args.out_root +'.prune.in'
    if not os.path.isfile(pruned_variants) or mapcount(pruned_variants) < 1 or args.force:
        cmd = f"plink2  --out {args.out_root} --read-freq  {args.bed.replace('.bed','.afreq')} --threads {cpus} --bfile {args.bed.replace('.bed','')} --extract {args.snplist} --indep-pairwise {' '.join(map(str,args.ld))} {args.pargs}"
        subprocess.call(shlex.split(cmd))
    print(f"Variants left after pruning: {mapcount(pruned_variants)}")

    if  mapcount(pruned_variants) < args.target_snps[0]:
        print('Starting ld params too stringent')
        return
    elif args.target_snps[0] <= mapcount(pruned_variants) <= args.target_snps[1]:
        print('done.')
        return
    else:
        # if number of SNPS is too high, iteratively ldprune lowering the threshold and using the previous pruned variants as the starting set
        args.ld[-1] = round(float(args.ld[-1]) - args.step,trailing_zeros(args.step))
        while not (args.target_snps[0] < mapcount(pruned_variants) < args.target_snps[1]) and float(args.ld[-1]) >= args.step:           
            args.force = True
            args.snplist = pruned_variants
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

    ld_pruning(args)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LD pruning of a bed file")

     # mutually exclusive group. Either calculate the info score variants or pass directly pass the variants.
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--info',nargs =2, metavar = ('FILE','VALUE'), help = 'Info score file and filter value')
    group.add_argument("--extract", type=file_exists, help =  "Path to list of variants to include")

    # required args
    parser.add_argument("-b", '--bed', metavar='F', type=file_exists, help="BED filepath", required=True)
    parser.add_argument('-o', "--out-path", type=str, help="folder in which to save the results", required=True)
    parser.add_argument('--prefix',  type=str, help="Output prefix", required=True)

    # optional args
    parser.add_argument('--pargs',type = str,help='extra plink args',default = '')
    parser.add_argument('--force',help='Flag on whether to force run',action = "store_true")
    parser.add_argument('--step',help="step by which to decrease the r2",type = float,default = 0.1)
    parser.add_argument('--target-snps',type = int,nargs = 2,help = "Target number of snps after ld",default = [100000,150000])
    parser.add_argument('--ld',nargs=3,metavar = ('SIZE','STEP','THRESHOLD'),help ='size,step,threshold',default = [1000,200,0.7])
    args = parser.parse_args()
    make_sure_path_exists(args.out_path)

    
    args.variants_path = os.path.join(args.out_path,'variants')
    make_sure_path_exists(args.variants_path)
    args.out_root = os.path.join(args.variants_path,args.prefix )

    main(args)
