version 1.0

workflow outlier_origin {

  input {
    String docker
    String chip_outlier_root
    String panel_root
    Boolean plot
  }
  # PRUNING OF CHIP OUTLIER DATA
  call prune_chip {
    input :
    docker=docker,
    plink_path = chip_outlier_root,
  }

  # ASSIGN POP/SUBPOP TO CHIP OUTLIERS
  Array[File] chip_plink = [chip_outlier_root + ".bed",chip_outlier_root + ".fam",chip_outlier_root + ".bim"]
  Array[File] panel_plink = [panel_root + ".bed",panel_root + ".fam",panel_root + ".bim"]

  call assign_origins as first_round_chip {
    input:
    docker=docker,
    reference_plink = panel_plink,
    candidate_plink = chip_plink,
    name= "first_round_chip",
    plot=plot,
    snps = prune_chip.snplist
  }

  # FROM MERGED PLINK OUTLIER DATA, REMOVE CHIP VARIANTS AND SPLIT CHIP/LEGACY SAMPLES
  call filter_imputed {
    input:
    docker=docker,
    chip_snps=chip_outlier_root + ".bim",
    chip_outliers_fam = chip_outlier_root + ".fam"
  }

  # ASSIGN ORIGIN TO LEGACY SAMPLES BASED ON PREVIOUS CHIP SAMPLES' ASSIGNMENT
  call assign_origins as second_round_imputed {
    input:
    docker=docker,
    reference_plink = filter_imputed.chip_outliers_imputed,
    candidate_plink = filter_imputed.legacy_outliers_imputed,
    reference_regions=first_round_chip.chip_origin,
    plot = plot,
    name= "second_round_imputed"
  }
}

task assign_origins {
  input {
    String docker
    Array[File] reference_plink
    Array[File] candidate_plink
    File reference_regions
    File? snps
    String name
    Boolean plot
  }
  Int disk_size = (ceil(size(reference_plink,'GB')) + ceil(size(candidate_plink,'GB')))
  
  command <<<
  python3.7 /scripts/project_ethnic.py --ref-bed ~{reference_plink[0]} --proj-bed ~{candidate_plink[0]} --sample-info ~{reference_regions} -o . --name ~{name}  --merge ~{if defined(snps) then "--extract " + snps else ""} ~{if (plot) then "--plot" else ""}

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
    Array[File] chip_plots = glob("./plot/*")
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
    String plink_path
    Int target
    String pargs
    String ld_params
  }
  File bed_file = plink_path + ".bed"
  File bim_file = plink_path + ".bim"
  File fam_file = plink_path + ".fam"
  File freq_file = plink_path + ".afreq"
  Int disk_size = ceil(size(bed_file,'GB'))*2
  String prefix = "chip"
  command {
  python3 /scripts/prune.py --bed ~{bed_file}  --prefix ~{prefix}  --ld ~{ld_params}  --target ~{target}   --pargs ~{pargs}  --out-path "/cromwell_root/" 
  }
  
  runtime {
    docker: "~{docker}"
    cpu: 32
    disks: "local-disk ~{disk_size} HDD"
    zones: "europe-west1-b europe-west1-c europe-west1-d"
    memory: "8 GB"
    preemptible: 1
  }
  
  output {
    File  snplist = "/cromwell_root/data/~{prefix}.prune.in"
  }
}

  

