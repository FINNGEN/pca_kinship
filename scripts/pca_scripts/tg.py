import os, shlex,subprocess
from utils import mem_mib,mapcount


#thousand genome project downloads
tg_download = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/'
tg_root,tg_ending  = 'ALL.chr','_GRCh38.genotypes.20170504'

def merge_1k(args):
    '''
    Merge finngen plink with tg plink so I can run kinship.
    '''
    args.merged_plink_file = args.plink_path + '1k_' + args.name
    if not os.path.isfile(args.merged_plink_file +'.bed') or args.force:
        args.force = True
        cmd = f'plink --bfile {args.bed.replace(".bed","")} --bmerge {args.new_tg}  --memory {int(mem_mib)} --make-bed --out {args.merged_plink_file} '
        args.v_print(3,cmd)
        subprocess.call(shlex.split(cmd))
        cmd = f'plink2 --bfile {args.merged_plink_file}   --memory {int(mem_mib)} --freq --out {args.merged_plink_file} '
        args.v_print(3,cmd)
        subprocess.call(shlex.split(cmd))

    else:
        args.v_print(3,'thousand genomes/finngen shared plink file already generated.')

    print(f"Total SNPs : {mapcount(args.merged_plink_file + '.bim')}")

def subset_1k(args):

    '''
     Subset 1kg data to smaller data set to speed up merging
    '''
    args.new_tg = args.plink_path + '1k_new'
    if not os.path.isfile(args.new_tg +'.bed') or args.force:
        print(args.new_tg +'.bed')
        args.force = True
        cmd = f'plink2  --bfile {args.tg_bed.replace(".bed","")}  --extract {args.bed.replace(".bed",".bim")} --memory {int(mem_mib)} --make-bed --out {args.new_tg} '
        print(cmd)
        subprocess.call(shlex.split(cmd))
    else:
        args.v_print(3,'thousand genomes subset plink file already generated.')


    
