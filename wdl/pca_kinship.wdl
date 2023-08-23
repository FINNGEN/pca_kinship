version 1.0

workflow pca_kinship {
  input {
    String prefix
    Boolean test
    String docker
    File min_pheno
    File sample_data 
    Array[String] chrom_list =  ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]
  }
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

  call merge_plink {
    input:
    bed_files = chrom_convert.bed,
    bim_files = chrom_convert.bim,
    fam_files = chrom_convert.fam,
    docker = docker,
    name = prefix,
    test= test,
  }

  # filters 1k data to subset of snps. Speeds up massively the transfer to pca
  call filter_tg{
    input: docker = docker,prefix = prefix,snplist = prune_panel.snplist
  }

  call kinship{
    input:
    bed_file = merge_plink.out_plink_merged_bed,
    fam_file = merge_plink.out_plink_merged_fam,
    bim_file = merge_plink.out_plink_merged_bim,
    prefix = prefix,
    pheno_file = min_pheno,
    docker = docker,
    metadata = sample_data,
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
  
  input {
    File bed_file
    File bim_file
    File fam_file
    File  freq_file
    
    File tg_bed
    File tg_fam
    File tg_bim
    
    File sample_file
    File kin_file
    File metadata
    
    String prefix
    String docker
    String? pca_docker
    Int cpu
  }
    
  Int disk_size =   ceil(size(bed_file,"GB"))*6 + ceil(size(tg_bed,"GB")) + 100
  Int mem = ceil(size(bed_file,"GB"))*3 + 10
  
  String out_path = "/cromwell_root/"
  String out_file = prefix + "_output.log"

  String? final_docker = if defined(pca_docker) then pca_docker else docker
  command <<<
  df -h
    python3.7 /scripts/pca.py \
    --bed  ~{bed_file} \
    --tg-bed ~{tg_bed} \
    --kin ~{kin_file} \
    --sample-info ~{sample_file} \
    --name ~{prefix} \
    --meta ~{metadata} \
    -o ~{out_path} \
    --release |& tee ~{out_file}
    
    mv ~{out_file} /cromwell_root/documentation/
    
  >>>
  runtime {
    docker: "~{final_docker}"
    cpu: "~{cpu}"
    disks:   "local-disk ~{disk_size} HDD"
    bootDiskSizeGb: 20
    memory: "~{mem} GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 0
  }
  
  output {
    File readme = "~{out_path}/${prefix}_pca_readme"
    #DATA
    Array[File] data =    glob("/cromwell_root/data/${prefix}*")
    #DOCUMENTATION
    Array[File] doc =    glob("/cromwell_root/documentation/${prefix}*")
  }
}

task kinship{
  input {
    File bed_file
    File bim_file
    File fam_file
    File pheno_file
    File metadata
    String docker
    String? kinship_docker
    String prefix
    Int cpu 
  }
  
  Int disk_size = ceil(size(bed_file,"GB"))*4 + 20
  Int mem = ceil(size(bed_file,"GB")) + 20

  String? final_docker = if defined(kinship_docker) then kinship_docker else docker
  command {
    python3.7  /scripts/ped.py \
    --bed ~{bed_file} \
    --out-path . \
    --prefix ~{prefix} \
    --meta ~{metadata} \
    --pheno-file ~{pheno_file} \
    --release
    
    ls ./data/
    
  }
  
  runtime {
    docker: "~{final_docker}"
    cpu: "~{cpu}"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    disks: "local-disk ~{disk_size} HDD"
    bootDiskSizeGb: 20
    memory: "~{mem} GB"
    preemptible: 0
  }
  
  output {
    File readme = "/cromwell_root/~{prefix}_kinship_readme"
    
    # DATA
    Array[File] data  = glob("/cromwell_root/data/~{prefix}*")	
    #DOCUMENTATION
    Array[File] doc = glob("/cromwell_root/documentation/~{prefix}*")

    File bed = "./data/${prefix}_kinship.bed"
    File fam = "./data/${prefix}_kinship.fam"
    File bim = "./data/${prefix}_kinship.bim"
    File freq = "./data/${prefix}_kinship.afreq"
    File kin = "./data/${prefix}.kin0"  
    
  }
}

task merge_plink {
  input {
    Array[File] bed_files 
    Array[File] bim_files
    Array[File] fam_files 
    File? exclusion_list
    String name
    String pargs
    
    String docker
    String? merge_docker
    Int mem
    Int cpu
    Boolean test
  }
  
  Int disk_size = ceil(size(bed_files,"GB"))*4+20
  Int plink_mem = mem*1000 - 2000

  String? final_docker = if defined(merge_docker) then merge_docker else docker
  command <<<
    cat ~{write_lines(bed_files)} | sed -e 's/.bed//g' > merge_list.txt
    cat ~{fam_files[0]} | shuf | head -n 50000  > test_fam.txt 
    plink --merge-list merge_list.txt ~{pargs} --keep-allele-order --memory ~{plink_mem} --make-bed --out ~{name} ~{if defined(exclusion_list) then  " --remove " + exclusion_list else "" } ~{if test then " --keep test_fam.txt " else ""}
  >>>
  
  runtime {
    docker: "~{final_docker}"
    cpu: "~{cpu}"
    disks: "local-disk ~{disk_size} HDD"
    bootDiskSizeGb: 20
    memory:  "~{mem}  GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 0
  }
  
  output {
    File out_plink_merged_bed  = "~{name}.bed"
    File out_plink_merged_bim  = "~{name}.bim"
    File out_plink_merged_fam  = "~{name}.fam"
  }
}


task chrom_convert {

  input{
    String chrom
    String pargs
    File variants
    String chromPath
    
    Int mem
    Int cpu   
    String docker
    String? convert_docker
  }
 
  File cFile = sub(chromPath,"CHROM",chrom)
  File tabix = cFile + '.tbi'
  Int disk_size = ceil(size(cFile,"GB")) *4  + 20
  Int plink_mem = mem*1000 - 2000

  String? final_docker = if defined(convert_docker) then convert_docker else docker
  command <<<
    cat ~{variants} | grep chr~{chrom}_ |  awk -F "_" '{print $1"\t"$2"\t"$2}' > tmp && split -n l/~{cpu} -d tmp regions
    ls regions* | parallel -j ~{cpu} "bcftools view ~{cFile} -R {} -Oz -o chunk{}.vcf.gz && echo {}"
    bcftools concat -n -f <(ls chunk*vcf.gz) -Oz -o tmp.vcf.gz && rm chunk*
        
    plink2 --vcf tmp.vcf.gz \
    ~{pargs} \
    --memory ~{plink_mem}  \
    --extract ~{variants} \
    --make-bed \
    --out ~{chrom} 
    
  >>>

  runtime {
    docker: "~{final_docker}"
    cpu: "~{cpu}"
    disks: "local-disk ~{disk_size}  HDD"
    bootDiskSizeGb: 20
    memory:  "~{mem} GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 0
    
  }
  
   output {
     File bed         = "~{chrom}.bed"
     File bim         = "~{chrom}.bim"
     File fam         = "~{chrom}.fam"
   }
}
 
task prune_panel {

  input {
    String docker
    String? prune_docker
    
    String prefix
  
    File info_score
    Float info_filter
    String plink_path
    Int cpu
    Int mem
    
    String pargs
    String ld_params
    String target
    
  }

  File bed_file = plink_path + ".bed"
  File bim_file = plink_path + ".bim"
  File fam_file = plink_path + ".fam"
  File freq_file = plink_path + ".afreq"
  
  Int disk_size = ceil(size(bed_file,'GB'))*2
  String? final_docker = if defined(prune_docker) then prune_docker else docker

  command {
    python3 /scripts/prune.py \
    --bed ~{bed_file} \
    --info ~{info_score} ~{info_filter} \
    --prefix ~{prefix} \
    --ld ~{ld_params} \
    --target ~{target} \
    --pargs ~{pargs} \
    --release \
    --out-path "/cromwell_root/" \ 
  }
  
  runtime {
    docker: "~{final_docker}"
    cpu: "~{cpu}"
    disks: "local-disk ~{disk_size} HDD"
    bootDiskSizeGb: 20
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory: "~{mem} GB"
    preemptible: 1
  }
  
  output {
    File readme  = "/cromwell_root/~{prefix}_prune_readme"
    File  snplist = "/cromwell_root/data/~{prefix}.prune.in"
    File log = "/cromwell_root/documentation/~{prefix}.prune.log"

  }
}


task filter_tg {

  input {
    
    String tg_root
    File snplist
    String docker
    String? tg_docker
    String prefix
    Int cpu
    Int mem
    
  }
  
  File tg_bed = tg_root + '.bed'
  File tg_fam = tg_root + '.fam'
  File tg_bim = tg_root + '.bim'
  
  Int disk_size =  ceil(size(tg_bed,'GB'))*4 + 11
  String out_root = prefix + "_1kg"

  String? final_docker = if defined(tg_docker) then tg_docker else docker
  command {
    plink2 --bfile ~{sub(tg_bed,'.bed','')} \
    --extract ~{snplist}  \
    --make-bed \
    --out ~{out_root}
  }
  
  runtime {
    docker: "~{final_docker}"
    cpu: "~{cpu}"
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory: "~{mem} GB"
    preemptible: 0
  }
  
  output {
    File bed = out_root + '.bed'
    File fam = out_root + '.fam'
    File bim = out_root + '.bim'
  }
}
