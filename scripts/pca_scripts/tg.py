import os, shlex,subprocess
from utils import mem_mib,mapcount,tmp_bash


def tg_bed(args):
    '''
    Merge finngen plink with tg plink so I can run .
    '''

    merged_plink_file = args.plink_path + '1k_' + args.name
    
    # build sample lists
    fg_samples = merged_plink_file + ".fg.samples"
    if not os.path.isfile(fg_samples) or args.force :
        args.force = True
        cmd = f'plink2 --bfile {args.bed.replace(".bed","")} --memory {int(mem_mib)} --write-samples --no-id-header --out {merged_plink_file} && cut -f1 {merged_plink_file}.id > {fg_samples}'
        args.logging.debug(cmd)
        tmp_bash(cmd)
       
    if not os.path.isfile(merged_plink_file +'.bed') or args.force:
        args.force = True
        
        args.logging.info(f"Building new plink data {merged_plink_file}")
        extract = f" {args.extract}" if args.extract else f'{args.bed.replace(".bed",".bim")}'
        cmd = f'plink --bfile {args.bed.replace(".bed","")} --bmerge {args.tg_bed.replace(".bed","")}  --memory {int(mem_mib)} --make-bed --out {merged_plink_file} --extract  {extract}' 
        args.logging.debug(cmd)
        subprocess.call(shlex.split(cmd))

    else:
        args.logging.info(f'{merged_plink_file} already generated.')

    if not os.path.isfile(merged_plink_file +'.afreq') or args.force:
        #args.force = True
        cmd = f'plink2 --bfile {merged_plink_file}   --memory {int(mem_mib)} --write-samples --no-id-header --freq --out {merged_plink_file}.tmp '
        args.logging.debug(cmd)
        subprocess.call(shlex.split(cmd))
        os.rename(f"{merged_plink_file}.tmp.afreq",f"{merged_plink_file}.afreq")
        cmd = f"cut -f1 {merged_plink_file}.tmp.id > {merged_plink_file}.id && rm {merged_plink_file}.tmp.id"
        tmp_bash(cmd)
        
    print(f"Total SNPs : {mapcount(merged_plink_file + '.bim')}")

    return merged_plink_file
