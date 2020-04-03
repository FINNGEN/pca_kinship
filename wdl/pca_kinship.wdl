workflow pca_kinship {

    String prefix
    String docker
    File min_pheno
    File sample_data 
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
        total_size = sum.out
        }
        
    call kinship{
        input:
        bed_file = merge_plink.out_plink_merged_bed,
        fam_file = merge_plink.out_plink_merged_fam,
        bim_file = merge_plink.out_plink_merged_bim,
        prefix = prefix,
        pheno_file = min_pheno,
        docker = docker,
	metadata = sample_data
        }


    # filters 1k data to subset of snps. Speeds up massively the transfer to pca
    call filter_tg{
        input: docker = docker,prefix = prefix,snplist = prune_panel.snplist
        }
    
    call pca {
        input :
        docker = docker,
        metadata = min_pheno,
        prefix = prefix,
        # kinship data
        bed_file = kinship.bed,
        fam_file = kinship.fam,
        bim_file = kinship.bim,
        freq_file = kinship.freq,
        kin_file = kinship.kin,
        # tg data
        tg_bed = filter_tg.bed,
        tg_fam = filter_tg.fam,
        tg_bim = filter_tg.bim,
	sample_file = sample_data,
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
    File metadata
    String prefix
    
    # runtime data
    String docker
    String pca_docker
    String final_docker  = if pca_docker != "" then pca_docker else docker
    Int disk_size = ceil(size(bed_file,'GB'))*6 + ceil(size(tg_bed,'GB')) + 10
    Int mem = ceil(size(bed_file,'GB')) + 10
    Int cpu
    
    String out_path = '/cromwell_root/'
    String out_file = prefix + '_output.log'
    command {
        python3 /scripts/pca.py \
	--bed ${bed_file} \
	--tg-bed ${tg_bed} \
	--kin ${kin_file} \
	--sample-info ${sample_file} \
	--name ${prefix} \
	--meta ${metadata} \
	-o ${out_path} \
	--release  |& tee ${out_file}

	        mv ${out_file} /cromwell_root/documentation/
        
    }
    
    runtime {
        docker: "${final_docker}"
        cpu: "${cpu}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory: "${mem} GB"
	preemptible: 0
    }
    output {
        File readme = "${out_path}/${prefix}_pca_readme"
        #DATA
	Array[File] data  = glob("/cromwell_root/data/${prefix}*")	
        #DOCUMENTATION
	Array[File] doc = glob("/cromwell_root/documentation/${prefix}*")
    }
}

task kinship{

    File bed_file
    File bim_file
    File fam_file
    File pheno_file
    File metadata

    String docker
    String kinship_docker
    String final_docker  = if kinship_docker != "" then kinship_docker else docker
    String prefix

    Int disk_size = ceil(size(bed_file,'GB'))*4 + 20
    Int mem = ceil(size(bed_file,'GB')) + 10
    Int cpu 
    String out_path = '/cromwell_root/'
    
    command {
        python3  /scripts/ped.py \
        --bed ${bed_file} \
        --out-path ${out_path} \
        --prefix ${prefix} \
        --pheno-file ${pheno_file} \
	--meta ${metadata} \
        --release
    }
    
    runtime {
        docker: "${final_docker}"
        cpu: "${cpu}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory: "${mem} GB"
	preemptible: 0
    }

    output {
        File readme = "${out_path}/${prefix}_kinship_readme"
        # DATA
	Array[File] data  = glob("/cromwell_root/data/${prefix}*")	
        File bed = "${out_path}/data/${prefix}_kinship.bed"
        File fam = "${out_path}/data/${prefix}_kinship.fam"
        File bim = "${out_path}/data/${prefix}_kinship.bim"
        File freq = "${out_path}/data/${prefix}_kinship.afreq"
        File new_fam = "${out_path}/data/${prefix}_pedigree.fam" 
        File kin = "${out_path}/data/${prefix}.kin0"
        File con = "${out_path}/data/${prefix}.con"
        #DOCUMENTATION
	Array[File] doc = glob("/cromwell_root/documentation/${prefix}*")
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
    Int cpu
    Int mem
    
    Int disk_size =  ceil(size(tg_bed,'GB'))*4 + 10
    String out_root = prefix + "_1kg"

    command {
        plink2 --bfile ${sub(tg_bed,'.bed','')} --extract ${snplist}  \
        --make-bed --out ${out_root}
        }

    runtime {
        docker: "${docker}"
        cpu: "${cpu}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory: "${mem} GB"
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
    String prune_docker
    String final_docker =  if prune_docker != ""  then prune_docker else docker
    String prefix

    File info_score
    Float info_filter
    String plink_path
    Int cpu
    Int mem

    String pargs
    String ld_params
    String target
    
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
        --target ${target} \
        --pargs ${pargs} \
        --out-path "/cromwell_root/" \ 
        }

    runtime {
        docker: "${final_docker}"
        cpu: "${cpu}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory: "${mem} GB"
	preemptible: 1
    }

    output {
        File snplist = "/cromwell_root/${prefix}.prune.in"
        }
}


task merge_plink {

    Array[File] bed_files 
    Array[File] bim_files
    Array[File] fam_files 
    
    String name
    String pargs

    Int mem
    Int cpu
    
    String merge_docker
    String docker
    String final_docker = if merge_docker != ""  then merge_docker else docker

    Int total_size 
    Int disk_size = total_size*4 + 20
    Int plink_mem = mem*1000 - 2000

    command <<<
    cat ${write_lines(bed_files)} | sed -e 's/.bed//g' > merge_list.txt
    plink --merge-list merge_list.txt ${pargs} --keep-allele-order --memory ${plink_mem} --make-bed --out ${name} 
    >>>

    runtime {
        docker: "${final_docker}"
	cpu: "${cpu}"
        disks: "local-disk " + "${disk_size}" + " HDD"
        bootDiskSizeGb: 20
        memory:  "${mem}" + " GB"
        preemptible: 0
    }

    output {    
       File out_plink_merged_bed  = "${name}.bed"
       File out_plink_merged_bim  = "${name}.bim"
       File out_plink_merged_fam  = "${name}.fam"
    }
}

task chrom_convert {

    String chrom
    String pargs
    Int mem
    Int cpu

    String docker
    String convert_docker
    String final_docker = if convert_docker != "" then convert_docker else docker

    File variants
    # get path to vcf
    String chromPath
    File cFile = sub(chromPath,"CHROM",chrom)
    File tabix = cFile + '.tbi'
    Int disk_size = ceil(size(cFile,"GB")) *4  + 20
    Int plink_mem = mem*1000 - 2000
    
    command <<<
    cat ${variants} | grep chr${chrom}_ |  awk -F "_" '{print $1"\t"$2"\t"$2}' > tmp && split -n l/${cpu} -d tmp regions
    ls regions* | parallel -j ${cpu} "bcftools view ${cFile} -R {} -Oz -o chunk{}.vcf.gz && echo {}"
    bcftools concat -n -f <(ls chunk*vcf.gz) -Oz -o tmp.vcf.gz && rm chunk*
    
    plink2 --vcf tmp.vcf.gz \
    ${pargs} \
    --memory ${plink_mem}  \
    --extract ${variants} \
    --make-bed \
    --out ${chrom} 
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
        Int file_size = ceil(size(cFile,'GB'))
    }
}



task sum {

  Array[Int] values
  String docker

  command <<<
    python -c "print(${sep="+" values})"
  >>>

  output {
    Int out = read_int(stdout())
  }

  runtime {
    docker: "${docker}"
    memory: "1 GB"
    cpu: 1
    maxRetries: 1
  }
}
