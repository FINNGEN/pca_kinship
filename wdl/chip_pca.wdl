version development

workflow pca_chip {

  input {
    String prefix
    Boolean test
    String docker
    File min_pheno
    File sample_data 
    String chip_plink_root
  }

  # no futher QC, we start pruning directly
  File chip_bim = chip_plink_root + ".bim"
  call prune_panel {
    input:
    docker = docker,
    filter_bim = chip_bim,
    prefix = prefix              
  }
  call filter_plink {
    input :
    test = test,
    snplist = prune_panel.snplist,
    chip_plink_root = chip_plink_root,
    docker = docker
  }

   call filter_tg{input: docker = docker,prefix = prefix,snplist = prune_panel.snplist}

  call kinship{
    input:
    test = test,
    bed_file = filter_plink.pruned_bed,
    fam_file = filter_plink.pruned_fam,
    bim_file = filter_plink.pruned_bim,
    prefix = prefix,
    pheno_file = min_pheno,
    metadata = sample_data,
    docker = docker,
  }
  call pca {
    input :
    docker = docker,
    metadata = min_pheno,
    sample_file = sample_data,
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
  }
}

task pca {
  
  input {
    File bed_file
    File bim_file
    File fam_file
    File freq_file
    
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
    Boolean test
    File bed_file
    File bim_file
    File fam_file
    File pheno_file
    File metadata
    String docker
    String? kinship_docker
    String prefix
  }

  Int cpu = if test then 32 else 95
  Int disk_size = ceil(size(bed_file,"GB"))*4 + 20
  Int mem = ceil(size(bed_file,"GB")) + 20

  String? final_docker = if defined(kinship_docker) then kinship_docker else docker
  command {
    python3.7  /scripts/ped.py \
    --bed ~{bed_file} \
    --out-path . \
    --prefix ~{prefix} \
    --pheno-file ~{pheno_file} \
    --meta ~{metadata} \
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
    File freq= "./data/${prefix}_kinship.afreq"
    File kin = "./data/${prefix}.kin0"  
    
  }
}

task filter_plink {
  
  input {
    Boolean test
    String docker
    String chip_plink_root
    File snplist
  }
  
  File chip_bed = chip_plink_root + ".bed"
  File chip_bim = chip_plink_root + ".bim"
  File chip_fam = chip_plink_root + ".FID.fam"

  Int mem = ceil(size(chip_bed,"GB"))
  Int disk_size = mem*3 + 10
  Int plink_mem = mem*1000 - 2000
  
  String out_root = basename(chip_bed,'.bed') + "_pruned"
  command <<<

   cat ~{chip_fam} | shuf | head -n 50000  > test_fam.txt 
  plink2 --bed ~{chip_bed} --bim ~{chip_bim} --fam ~{chip_fam}  --memory ~{plink_mem} --extract ~{snplist} --make-bed --out ~{out_root} ~{if test then " --keep test_fam.txt " else ""}
  >>>
  
  runtime {
    docker: "~{docker}"
    cpu: "32"
    disks: "local-disk ~{disk_size} HDD"
    bootDiskSizeGb: 20
    memory:  "~{mem}  GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    preemptible: 0
  }
  output {
    File pruned_bed  = "~{out_root}.bed"
    File pruned_bim  = "~{out_root}.bim"
    File pruned_fam  = "~{out_root}.fam"
  }
  
}
task prune_panel {

  input {
    String docker
    String? prune_docker
    
    String prefix
    File filter_bim
    String panel_root   
    String pargs
    String ld_params
    Int target
    
  }

  File bed_file = panel_root + ".bed"
  File bim_file = panel_root + ".bim"
  File fam_file = panel_root + ".fam"
  File freq_file = panel_root + ".afreq"

  Int mem = ceil(size(bed_file,'GB'))
  Int disk_size = mem*2
  String? final_docker = if defined(prune_docker) then prune_docker else docker

  command {
    python3 /scripts/prune.py \
    --bed ~{bed_file} \
    --extract ~{filter_bim} \
    --prefix ~{prefix} \
    --ld ~{ld_params} \
    --target ~{target} \
    --pargs ~{pargs} \
    --out-path "/cromwell_root/" \ 
  }
  
  runtime {
    docker: "~{final_docker}"
    cpu: "32"
    disks: "local-disk ~{disk_size} HDD"
    bootDiskSizeGb: 20
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory: "~{mem} GB"
    preemptible: 1
  }
  
  output {
    File  snplist = "/cromwell_root/data/~{prefix}.prune.in"
  }
  
}


task filter_tg {

  input {
    
    String tg_root
    File snplist
    String docker
    String? tg_docker
    String prefix
    
  }
  
  File tg_bed = tg_root + '.bed'
  File tg_fam = tg_root + '.fam'
  File tg_bim = tg_root + '.bim'

  Int mem = ceil(size(tg_bed,'GB'))
  Int disk_size =  mem*4 + 10
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
    cpu: 16
    disks: "local-disk ~{disk_size} HDD"
    bootDiskSizeGb: 20
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
