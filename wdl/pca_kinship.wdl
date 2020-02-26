workflow pca_kinship {

    String prefix
    String docker
    Array[String] chrom_list 

    call prune_panel {
        input:
        docker = docker,
        prefix = prefix              
    }

    # create plink file base on pruned snplist
    scatter (chrom in chrom_list){
       call chrom_convert {
	   input :
           docker = docker,
	   chrom = chrom,
           variants = prune_panel.snplist           
       }     
    }
     
    # gathers the results from the file_size scatter
    call sum {input: values = chrom_convert.file_size,docker = docker}
    # merge vcf files
    call merge_plink {
        input:
	bed_files = chrom_convert.bed,
	bim_files = chrom_convert.bim,
	fam_files = chrom_convert.fam,
        docker = docker,
        name = prefix,
        total_size  = ceil(sum.out)
        }
        
    call kinship{
        input:
        bed_file = merge_plink.out_plink_merged_bed,
        fam_file = merge_plink.out_plink_merged_fam,
        bim_file = merge_plink.out_plink_merged_bim,
        freq_file = merge_plink.out_plink_merged_freq,
        prefix = prefix,
        docker = docker,
        }
}

task kinship{

    File bed_file
    File bim_file
    File fam_file
    File freq_file
    File pheno_file
    File? sub_fam

    String docker
    String? kinship_docker
    String? final_docker  = if defined(kinship_docker) then kinship_docker else docker
    String prefix

    Int disk_size = ceil(size(bed_file,'GB'))*6 + 10
    Int mem = ceil(size(bed_file,'GB')) + 10
    Int cpus = if defined(sub_fam) then 16 else 64

    String out_path = '/cromwell_root/kinship'
    command {
        python3  /scripts/ped.py \
        --bed ${bed_file} \
        ${'--fam '  + sub_fam} \
        -o ${out_path} \
        --prefix ${prefix} \
        --pheno-file ${pheno_file} 
    }
    
    runtime {
        docker: "${final_docker}"
        cpu: "${cpus}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory: "${mem} GB"
	preemptible: 0
    }

    output {
        Array[File] log_files = glob("${out_path}/.log")
        Array[File] pdf_files = glob("${out_path}/*.pdf")
        File kinship_bed =  "${out_path}/plink/${prefix}.kinship.bed"
        File kinship_fam =  "${out_path}/plink/${prefix}.kinship.fam"
        File pedigree_fam =  "${out_path}/plink/${prefix}.kinship.pedigree.fam"
        File kinship_bim =  "${out_path}/plink/${prefix}.kinship.bim"
        File kinship_afreq =  "${out_path}/plink/${prefix}.kinship.afreq"
        File kinship = "${out_path}/${prefix}.kin0"
        Array[File] pedigree_files = glob("${out_path}/*pedigree*")
        }
}

task prune_panel {

    String docker
    String? prune_docker
    String? final_docker =  if defined(prune_docker) then prune_docker else docker
    String prefix

    File info_score
    Float info_filter
    String plink_path
    Int cpus

    String pargs
    String ld_params
    String target_snps
    String step
    
    File bed_file = plink_path + ".bed"
    File bim_file = plink_path + ".bim"
    File fam_file = plink_path + ".fam"
    File freq_file = plink_path + ".afreq"

    Int disk_size = ceil(size(bed_file,'GB'))*2
    
    command {
        python3 /scripts/prune.py \
        --bed ${bed_file} \
        --info ${info_score} ${info_filter} \
        --prefix ${prefix} \
        --ld ${ld_params} \
        --step ${step} \
        --target-snps ${target_snps} \
        --pargs ${pargs} \
        -o "/cromwell_root/" \ 
        }


    runtime {
        docker: "${final_docker}"
        cpu: "${cpus}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory: "16 GB"
	preemptible: 1
    }

    output {

        File snplist = "/cromwell_root/variants/${prefix}.prune.in"
        }
}



task merge_plink {

    Array[File] bed_files 
    Array[File] bim_files
    Array[File] fam_files 
    
    String name
    String pargs
    Int total_size
    Int disk_factor
    Int mem
    
    String? merge_docker
    String docker
    String? final_docker = if defined(merge_docker) then merge_docker else docker
    
    Int disk_size = total_size*disk_factor + 20
    Int plink_mem = mem*1000 - 2000

    command <<<
    cat ${write_lines(bed_files)} | sed -e 's/.bed//g' > merge_list.txt
    plink --merge-list merge_list.txt ${pargs} --memory ${plink_mem} --make-bed --out ${name}
    plink2 --bfile ${name} --keep-allele-order --freq --out ${name}
    >>>

    runtime {
        docker: "${final_docker}"
	cpu: 16
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory:  "${mem}" + " GB"
        preemptible: 0

    }

    output {    
       File out_plink_merged_bed  = "${name}.bed"
       File out_plink_merged_bim  = "${name}.bim"
       File out_plink_merged_fam  = "${name}.fam"
       File out_plink_merged_freq = "${name}.afreq"
    }
}

task chrom_convert {

    String chrom
    String pargs
    Int disk_factor
    Int mem
    Int cpu

    String docker
    String? convert_docker
    String? final_docker = if defined(convert_docker) then convert_docker else docker

    File variants
    # get path to vcf
    String chromPath
    File cFile = sub(chromPath,"CHROM",chrom)
    Int chrom_size = ceil(size(cFile,"GB"))   
    
    Int disk_size = disk_factor * chrom_size + 20
    Int plink_mem = mem*1000 - 2000

    
    command <<<
    plink2 --vcf ${cFile} \
    ${pargs} --vcf-half-call h \
    --memory ${plink_mem}  --extract ${variants} \
    --make-bed --out ${chrom} 
    >>>

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
	memory:  "${mem}" + " GB"
        preemptible: 0

    }

    output {	
	File bed         = "${chrom}.bed"
	File bim         = "${chrom}.bim"
	File fam         = "${chrom}.fam"
        Int file_size = chrom_size
    }
}
task sum {

  Array[Float] values
  String docker

  command <<<
    python -c "print(${sep="+" values})"
  >>>

  output {
    Float out = read_float(stdout())
  }

  runtime {
    docker: "${docker}"
    memory: "1 GB"
    cpu: 1
    maxRetries: 1
  }
}
