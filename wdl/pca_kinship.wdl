workflow pca_kinship {

    String prefix
    String docker
    
    Array[String] chrom_list =  ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
    
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


    # filters 1k data to subset of snps. Speeds up massively the transfer to pca
    call filter_tg{
        input: docker = docker,prefix = prefix,snplist = prune_panel.snplist
        }
    
    call pca {
        input :
        docker = docker,
        bed_file = kinship.bed,
        fam_file = kinship.fam,
        bim_file = kinship.bim,
        freq_file = kinship.freq,
        prefix = prefix,
        kin_file = kinship.kin,
        tg_bed = filter_tg.bed,
        tg_fam = filter_tg.fam,
        tg_bim = filter_tg.bim,
        }
}


task pca {

    #finngen plink data
    File bed_file
    File bim_file
    File fam_file
    File freq_file
    # filtered tg plink data
    File tg_bed 
    File tg_fam 
    File tg_bim 
    # sample metadata
    File sample_file
    File kin_file
    String prefix
    # runtime data
    String docker
    String? pca_docker
    String? final_docker  = if defined(pca_docker) then pca_docker else docker
    Int disk_size = ceil(size(bed_file,'GB'))*6 + ceil(size(tg_bed,'GB')) + 10
    Int mem = ceil(size(bed_file,'GB')) + 10
    
    command {
        python3 /scripts/pca_main.py --bed ${bed_file} --tg-bed ${tg_bed}  -k ${kin_file}  -s ${sample_file}  --name ${prefix}  -o .
    }
    
    runtime {
        docker: "${final_docker}"
        cpu: 32
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory: "${mem} GB"
	preemptible: 0
    }
    output {
        Array[File] out = glob('/cromwell_root/plots/*pdf')
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

    Int disk_size = ceil(size(bed_file,'GB'))*4 + 20
    Int mem = ceil(size(bed_file,'GB')) + 10
    Int cpus = if defined(sub_fam) then 16 else 64

    String out_path = '/cromwell_root/'
    command {
        python3  /scripts/ped.py \
        --bed ${bed_file} \
        ${'--fam '  + sub_fam} \
        -o ${out_path} \
        --prefix ${prefix} \
        --pheno-file ${pheno_file} \
        --release
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
        # DATA
        File bed = "${out_path}/data/${prefix}_kinship.bed"
        File fam = "${out_path}/data/${prefix}_kinship.fam"
        File bim = "${out_path}/data/${prefix}_kinship.bim"
        File freq = "${out_path}/data/${prefix}_kinship.afreq"
        File new_fam = "${out_path}/data/${prefix}_pedigree.fam" 
        File kin = "${out_path}/data/${prefix}.kin0"
        File duplicates = "${out_path}/data/${prefix}.con"
        #DOCUMENTATION
        File log = "${out_path}/documentation/${prefix}.log"
        File kinship_log = "${out_path}/documentation/${prefix}_kinship.log"
        File pedigree_log = "${out_path}/documentation/${prefix}_pedigree.log"
        File duplicateplot = "${out_path}/documentation/${prefix}_duplicateplot.pdf"
        File relplot = "${out_path}/documentation/${prefix}_relplot.pdf"
        File famplot = "${out_path}/documentation/${prefix}_uniqfam0plot.pdf"
        File kinplot = "${out_path}/documentation/${prefix}_kinship_distribution.pdf"
        }
}


task filter_tg {

    String tg_root
    File tg_bed = tg_root + '.bed'
    File tg_fam = tg_root + '.fam'
    File tg_bim = tg_root + '.bim'
    File snplist
    
    String docker
    String prefix
    Int disk_size =  ceil(size(tg_bed,'GB'))*4 + 10
    String out_root = prefix + "_1kg"
    
    command {
        plink2 --bfile ${sub(tg_bed,'.bed','')} --extract ${snplist}  \
        --make-bed --out ${out_root}
        }

    runtime {
        docker: "${docker}"
        cpu: 32
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory: "16 GB"
	preemptible: 0
    }

    output {
        File bed = out_root + '.bed'
        File fam = out_root + '.fam'
        File bim = out_root + '.bim'
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
    Int mem
    
    String? merge_docker
    String docker
    String? final_docker = if defined(merge_docker) then merge_docker else docker
    
    Int disk_size = total_size*4 + 20
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
    Int disk_size = chrom_size*4  + 20
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
