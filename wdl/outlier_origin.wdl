version 1.0

workflow outlier_origin {

  input {
    String docker
    String chip_outlier_root
    String name
  }

  # PRUNING OF CHIP OUTLIER DATA
  call prune_chip {
    input :
    docker=docker,
    plink_path = chip_outlier_root,
    prefix=name
  }
  # ASSIGN POP/SUBPOP TO CHIP OUTLIERS
  call first_round_chip {
    input:
    docker=docker,
    snps = prune_chip.snplist,
    chip_outlier_root=chip_outlier_root
  }

  # FROM MERGED PLINK OUTLIER DATA, REMOVE CHIP VARIANTS AND SPLIT CHIP/LEGACY SAMPLES
  call filter_imputed {
    input:
    docker=docker,
    chip_snps=chip_outlier_root + ".bim",
    chip_outliers_fam = chip_outlier_root + ".fam"
  }

  # ASSIGN ORIGIN TO LEGACY SAMPLES BASED ON PREVIOUS CHIP SAMPLES' ASSIGNMENT
  call second_round_imputed {
    input:
    docker=docker,
    chip_outliers_imputed=filter_imputed.chip_outliers_imputed,
    legacy_outliers_imputed=filter_imputed.legacy_outliers_imputed,
    chip_outliers_regions=first_round_chip.chip_origin
  }
}



task first_round_chip {
  input {
    String chip_outlier_root
    String panel_root
    File panel_info
    File snps
    String docker
  }
  #LOCALIZE PANEL
  File panel_bed = panel_root + '.bed'
  File panel_bim = panel_root + '.bim'
  File panel_fam = panel_root + '.fam'
  # LOCALIZE CHIP
  File chip_bed = chip_outlier_root + '.bed'
  File chip_bim = chip_outlier_root + '.bim'
  File chip_fam = chip_outlier_root + '.fam'
  # RUNTIME
  Int disk_size = (ceil(size(panel_bed,'GB')) + ceil(size(chip_bed,'GB')))
  String name = "first_round"

  command <<<
  python3.7 /scripts/project_ethnic.py --ref-bed ~{panel_bed} --proj-bed ~{chip_bed} --sample-info ~{panel_info} --extract ~{snps} -o . --name ~{name} --plot --merge

  >>>
  runtime {
    docker: "~{docker}"
    cpu: 32
    disks:   "local-disk ~{disk_size*2 +10} HDD"
    memory: "~{disk_size} GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
  }
  output {
    File chip_origin = "~{name}_merged_samples_most_likely_region.txt"
    Array[File] chip_plots = glob("./plot/*.p*")
  }
}


task second_round_imputed {
  input {
    String docker
    Array[File] chip_outliers_imputed
    Array[File] legacy_outliers_imputed
    File chip_outliers_regions
  }
  File chip_bed = chip_outliers_imputed[0]
  File legacy_bed= legacy_outliers_imputed[0]
  Int disk_size = (ceil(size(chip_bed,'GB')) + ceil(size(legacy_bed,'GB')))
  String name = "second_round"

  command <<<
  python3.7 /scripts/project_ethnic.py --ref-bed ~{chip_bed} --proj-bed ~{legacy_bed} --sample-info ~{chip_outliers_regions} -o . --name ~{name} --plot --merge

  >>>
  runtime {
    docker: "~{docker}"
    cpu: 32
    disks:   "local-disk ~{disk_size*2 +10} HDD"
    memory: "~{disk_size} GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
  }
  output {
    File chip_origin = "~{name}_merged_samples_most_likely_region.txt"
    Array[File] chip_plots = glob("./plot/*.p*")
  }

}

# SECOND ROUND STEP WHERE I NEED TO WORK WITH THE KINSHIP DATA,SPLITING CHIP AND NON CHIP AND REMOVING CHIP VARIANTS
task filter_imputed {

  input {
    String docker
    String plink_outliers_path
    File chip_outliers_fam
    File chip_snps
  }

  File bed = plink_outliers_path + ".bed"
  File bim = plink_outliers_path + ".bim"
  File fam = plink_outliers_path + ".fam"

  Int disk_size = ceil(size(bed,'GB'))
  command <<<
  plink2 --bfile ~{sub(bed,'.bed','')} --exclude ~{chip_snps} --keep   ~{chip_outliers_fam} --make-bed --out chip_outliers_imputed
  plink2 --bfile ~{sub(bed,'.bed','')} --exclude ~{chip_snps} --remove ~{chip_outliers_fam} --make-bed --out legacy_outliers_imputed
  >>>
  runtime {
    docker: "~{docker}"
    cpu: 8
    disks:   "local-disk ~{disk_size*2+10} HDD"
    memory: "~{disk_size} GB"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
  }
  output {
    Array[File] chip_outliers_imputed =   ["chip_outliers_imputed.bed",  "chip_outliers_imputed.fam",  "chip_outliers_imputed.bim"]
    Array[File] legacy_outliers_imputed = ["legacy_outliers_imputed.bed","legacy_outliers_imputed.fam","legacy_outliers_imputed.bim"]
    
  }
}

task prune_chip {
  input {
    String docker
    String prefix
    String plink_path
    Int cpu
    Int mem
    Int target
    String pargs
    String ld_params
  }
  File bed_file = plink_path + ".bed"
  File bim_file = plink_path + ".bim"
  File fam_file = plink_path + ".fam"
  File freq_file = plink_path + ".afreq"
  Int disk_size = ceil(size(bed_file,'GB'))*2

   command {
    python3 /scripts/prune.py \
    --bed ~{bed_file} \
    --prefix ~{prefix} \
    --ld ~{ld_params} \
    --target ~{target} \
    --pargs ~{pargs} \
    --out-path "/cromwell_root/" \ 
  }
  
  runtime {
    docker: "~{docker}"
    cpu: "~{cpu}"
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

  

