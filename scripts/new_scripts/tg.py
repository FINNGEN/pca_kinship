import os, shlex,subprocess
from utils import mem_mib,mapcount


def merge_1k(args):
    '''
    Merge finngen plink with tg plink so I can run kinship.
    '''
    merged_plink_file = args.plink_path + '1k_' + args.name
    if not os.path.isfile(merged_plink_file +'.bed') or args.force:
        args.force = True
        args.logging.info(f"Building new plink data {merged_plink_file}")
        cmd = f'plink --bfile {args.bed.replace(".bed","")} --bmerge {args.tg_bed.replace(".bed","")}  --memory {int(mem_mib)} --make-bed --out {merged_plink_file} --extract {args.bed.replace(".bed",".bim")} '
        args.logging.debug(cmd)
        subprocess.call(shlex.split(cmd))

    else:
        args.logging.info(f'{merged_plink_file} already generated.')

    if not os.path.isfile(merged_plink_file +'.afreq') or args.force:
        args.force = True
        cmd = f'plink2 --bfile {merged_plink_file}   --memory {int(mem_mib)} --freq --out {merged_plink_file} '
        args.logging.debug(cmd)
        subprocess.call(shlex.split(cmd))

        
    print(f"Total SNPs : {mapcount(merged_plink_file + '.bim')}")
    return merged_plink_file
