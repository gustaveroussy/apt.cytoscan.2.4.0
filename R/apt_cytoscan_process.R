apt.cytoscan.process <- function(CEL = NULL, samplename = NULL, dual.norm = FALSE, normal.diploid = FALSE, out.dir = getwd(), temp.files.keep = FALSE, force.OS = NULL, apt.build = "na33.r4") {
  
  # setwd("/home/job/svn/genomics/CGH/R/00_PIPELINE/TEST_ZONE/CSHD")
  # CEL <- "M2271_K03.CEL"
  # samplename <- "M2271_K03"
  # dual.norm = FALSE
  # normal.diploid = FALSE
  # out.dir = getwd()
  # temp.files.keep = FALSE
  # force.OS = NULL
  # apt.build = "na33.r4"
  
  if (is.null(CEL)) stop("A CEL file is required !")
  if (is.null(samplename)) stop("A samplename is required !")
  if (!file.exists(CEL)) stop(paste0("Could not find CEL file ", CEL, " !"))
  if (!dir.exists(out.dir)) stop(paste0("Output directory [", out.dir, "] does not exist !"))
  
  out.dir <- tools::file_path_as_absolute(out.dir)
  CEL <- tools::file_path_as_absolute(CEL)
  
  ## Checking build compatibility
  knownbuilds <- c("na33.r2", "na33.r4", "na36.r1")
  if (!(tolower(apt.build) %in% knownbuilds)) warning(paste0(" WARNING : The requested build ", apt.build, " is not in the validated list. Program may crash / fail !"))
  
  ## Checking apt-copynumber-cyto-ssa package loc
  tool.version <- "2.4.0"
  self.pkg.name <- paste0("apt.cytoscan.", tool.version)
  bin.dir <- system.file("apt/bin/", package = self.pkg.name)
  
  ## Checking apt-copynumber-cyto annotation package loc
  res.pkg.name <- paste0("CytoScanHD.Array.", tolower(apt.build))
  if (!(res.pkg.name %in% utils::installed.packages())) stop(paste0("Package ", res.pkg.name, " not found !"))
  res.dir <- system.file("apt/res/", package = res.pkg.name)
  suppressPackageStartupMessages(require(res.pkg.name, character.only = TRUE))
  apt.files <- annotation.set.describe()
  
  ## Checking annotation files availability
  for (f in names(apt.files)) { if (!file.exists(paste0(res.dir, "/", apt.files[[f]]))) stop(paste0("File ", apt.files[[f]], " is not available for ", apt.build, " !")) }

  ## Checking the OS
  # message("Identying OS ...")
  os.list <- c("linux", "windows", "osx")
  my.os <- get.os()
  tmsg(paste0("OS is reported as ", my.os))
  if (!is.null(force.OS)) {
    if (!(force.OS %in% os.list)) stop("Specified forced OS is not supported !")
    my.os <- force.OS
    tmsg(paste0("WARNING : Forcing OS to : ", my.os))
  } else if (!(my.os %in% os.list)) stop(paste0("Current OS [", my.os, "] not supported ! If you are sure of your OS support, use force.OS option with any of 'linux', 'windows', 'osx'"))
  
  if (my.os == "windows") my.os <- paste0(my.os, ".exe")
  
  oridir <- getwd()
  out.dir.p <- paste0(out.dir, "/", samplename)
  out.dir.w <- paste0(out.dir.p, "/temp")
  dir.create(path = out.dir.w, recursive = TRUE)
  setwd(out.dir.p)
  
  cfile <- paste0(out.dir.p, "/", samplename, "_CELfile.txt")
  affycel.df <- data.frame(cel_files = CEL, stringsAsFactors = FALSE)
  write.table(affycel.df, file = cfile, quote = FALSE, sep = "\t", row.names = FALSE)

  apt.cmd <- c(paste0(bin.dir, "/apt-copynumber-cyto-ssa_", my.os, " "),
               "--allele-peak-count 24 ",
               "--nearest-power-of-two 512 ",
               "--keep-intermediate-data false ",
               paste0("--dual-channel-normalization ", tolower(as.character(dual.norm)), " "),
               "--mangle-probeset-names false ",
               "--waviness-block-size 50 ",
               "--waviness-genomic-span 0 ",
               "--igender-male-threshold 1.5 ",
               "--igender-female-threshold 0.9 ",
               "--sig-gender-xx-cutoff 0.61 ",
               "--sig-gender-xx-cutoff-high 0.95 ",
               "--sig-gender-y-cutoff 0.58 ",
               "--sig-gender-reference-chromosome 2 ",
               "--force false ",
               "--text-output true ",
               "--oschp-output true ",
               "--cychp-output false ",
               paste0("--set-analysis-name ", samplename, " "),
               paste0("--snp-qc-snp-list ", res.dir, "/", apt.files$snplist, " "),
               paste0("--x-probes-file ", res.dir, "/", apt.files$xprobes, " "),
               paste0("--y-probes-file ", res.dir, "/", apt.files$yprobes, " "),
               paste0("--annotation-file ", res.dir, "/", apt.files$annotdb, " "),
               paste0("--reference-file ", res.dir, "/", apt.files$refmodel, " "),
               "--y-target -0.50 ",
               paste0("--out-dir ", out.dir.w, " "),
               paste0("--log-file ", out.dir.w, "/apt-copynumber-cyhd-ssa.log "),
               paste0("--temp-dir ", out.dir.w, "/temp "),
               "--alpha-cn-calibrate 0.564278 ",
               "--alpha-X-cn-calibrate 0.619453 ",
               "--alpha-Y-cn-calibrate 0.494620 ",
               "--beta-cn-calibrate 1 ",
               "--beta-X-cn-calibrate 1 ",
               "--beta-Y-cn-calibrate 1 ",
               "--brlmmp-HARD 3 ",
               "--brlmmp-SB 0.450000 ",
               "--brlmmp-CM 2 ",
               "--brlmmp-bins 100 ",
               "--brlmmp-mix 1 ",
               "--brlmmp-bic 2 ",
               "--brlmmp-CSepPen 0.0 ",
               "--brlmmp-CSepThr 16 ",
               "--brlmmp-lambda 1.0 ",
               "--brlmmp-wobble 0.05 ",
               "--brlmmp-copyqc 0.000000 ",
               "--brlmmp-copytype 0 ",
               "--brlmmp-ocean 0.0 ",
               "--brlmmp-clustertype 1 ",
               "--brlmmp-transform mva ",
               "--brlmmp-MS 0.05 ",
               paste0("--log2-ratio-adjuster-image-node:enable=", tolower(as.character(!normal.diploid))," "),
               "--hpf-mini-block-rows 8 ",
               "--hpf-mini-block-cols 8 ",
               "--hpf-global-smooth-weight 256 ",
               "--hpf-local-smooth-weight 64 ",
               "--hpf-converged 0.000100 ",
               "--log2-ratio-adjuster-wave-node:enable=true ",
               "--wave-bandwidth 101 ",
               "--wave-bin-count 25 ",
               "--wave-count 6 ",
               "--wave-smooth true ",
               "--median-autosome-median-normalization true ",
               "--ad-outlier-trim 3.0 ",
               "--ap-ad-step 20 ",
               "--ap-ad-window 100 ",
               "--ap-ad-point-count 129 ",
               "--ap-ad-bandwidth 0.25 ",
               "--ap-ad-cutoff 0.05 ",
               "--ap-ad-threshold 0.35 ",
               "--ap-ad-symmetry true ",
               "--ap-ad-height-threshold 0.000000 ",
               "--ap-ad-height-threshold-bound 1024.000000 ",
               paste0("--normal-diploid-normalization ", tolower(as.character(normal.diploid)), " "),
               paste0("--fld-file ", res.dir, "/", apt.files$fld, " "),
               "--ndd-weights-type FLD ",
               "--ndd-point-count 128 ",
               "--ndd-bandwidth 0.25 ",
               "--ndd-cutoff 0.05 ",
               "--ndd-step 40 ",
               "--ndd-window 275 ",
               "--ndd-min-markers 10 ",
               "--ndd-log2-threshold 0.28 ",
               "--ndd-th1 0,0.4 ",
               "--ndd-th2 0.7,1.3 ",
               "--ndd-th3 0.7,1.3 ",
               "--ndd-th4 0,0.4 ",
               "--ndd-min-percent-nd 0 ",
               "--ndd-min-number-nd 8000 ",
               "--ndd-max-number-nd 1500000 ",
               "--ndd-med-l2r-low-percentile-ND 0.015 ",
               "--ndd-med-l2r-threshold-ND 0.20 ",
               "--ndd-iqr-factor-coarse 5.0 ",
               "--ndd-iqr-factor-fine 3.0 ",
               "--nddlrcov-min-nd-marker-count-per-bin 30 ",
               "--nddlrcov-use-nd-lr-markers on ",
               "--nddqc-min-count 2000 ",
               "--nddicov-min-nd-probe-count-per-bin 30 ",
               "--nddicov-use-nd-signal-markers off ",
               "--kernel-sigma-span 50 ",
               "--hmmCN_state 0,1,2,3,4 ",
               "--hmmCN_mu -2,-0.5,0,0.33,0.54 ",
               "--hmmCN_sigma 0.3,0.3,0.3,0.3,0.3 ",
               "--hmmCN_state-X 0,1,2,3,4 ",
               "--hmmCN_mu-X -2,-0.55,0,0.36,0.57 ",
               "--hmmCN_sigma-X 0.3,0.3,0.3,0.3,0.3 ",
               "--hmmCN_state-Y 0,1,2,3,4 ",
               "--hmmCN_mu-Y -1.4,-0.5,0,0.3,0.51 ",
               "--hmmCN_sigma-Y 0.3,0.3,0.3,0.3,0.3 ",
               "--diagonal-weight-Y 0.995 ",
               "--mapd-weight-Y 0.2 ",
               "--min-segment-size-Y 5 ",
               "--hmm-confidence-weight-Y 0.6 ",
               "--diagonal-weight 0.995 ",
               "--mapd-weight 0.2 ",
               "--min-segment-size 5 ",
               "--hmm-confidence-weight 0.6 ",
               "--hmm-shrink true ",
               "--diagonal-weight-X 0.995 ",
               "--mapd-weight-X 0.2 ",
               "--min-segment-size-X 5 ",
               "--hmm-confidence-weight-X 0.6 ",
               "--agr-denominator-source CN ",
               "--agr-denominator-percentile 0.5 ",
               "--cn-gender-cutoff 0.5 ",
               "--loh-error-rate 0.05 ",
               "--loh-beta 0.001 ",
               "--loh-alpha 0.01 ",
               "--loh-separation 1000000 ",
               "--loh-min-marker-count 10 ",
               "--loh-no-call-threshold 0.05 ",
               "--loh-min-genomic-span 1000000 ",
               "--gains-boundries 0.0820,0.1588,0.2318,0.3130 ",
               "--losses-boundries -0.1049,-0.2263,-0.3477,-0.4691 ",
               "--marker-bandwidth 6000 ",
               "--confidence-window 251 ",
               "--run-y-chromosome false ",
               paste0("--cel-files ", cfile, " "),
               paste0("--ap-baf-fld-weights ", res.dir, "/", apt.files$fld, " "),
               "--cn-neutral-loh-node:enable=true"
               )
  tmsg("Running APT ...")
  cmdres <- try(system(command = paste0(apt.cmd, collapse = ""), intern = TRUE))
  
  oscf <- list.files(path = out.dir.w, pattern = "\\.oschp$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  
  if (!file.exists(oscf)) return(cmdres)

  logf <- list.files(path = out.dir.w, pattern = "\\.log$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  qcf <- list.files(path = out.dir.w, pattern = "\\.CopyNumberReport\\.txt$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)
  
  ## Renaming files
  tmsg("Renaming OSCHP ...")
  new.oscf <- paste0(out.dir.p, "/", samplename, "_", tool.version, "_", apt.build, ".oschp")
  new.logf <- paste0(out.dir.p, "/", samplename, "_", tool.version, "_", apt.build, ".log")
  new.qcf <- paste0(out.dir.p, "/", samplename, "_", tool.version, "_", apt.build, ".qc.txt")
  file.rename(from = oscf[1], to = new.oscf)
  file.rename(from = logf[1], to = new.logf)
  file.rename(from = qcf[1], to = new.qcf)
  
  ## Cleaning
  if(!temp.files.keep) {
    tmsg("Removing temporary files ...")
    unlink(x = out.dir.w, recursive = TRUE, force = TRUE)
  }
  setwd(oridir)
  
  tmsg("Done.")
  return(new.oscf)
}

apt.cytoscan.process.batch <- function(CEL.list.file = NULL, nthread = 1, cluster.type = "PSOCK", ...) {
  ## Checking the CEL.list.file
  if (is.null(CEL.list.file)) stop("A CEL.list.file is required !")
  if (!file.exists(CEL.list.file)) stop("Could not find CEL.list.file !")
  message("Reading and checking CEL.list.file ...")
  myCELs <- read.table(file = CEL.list.file, header = TRUE, sep="\t", check.names = FALSE, as.is = TRUE)
  head.ok <- c("cel_files", "SampleName")
  head.chk <- all(colnames(CEL.list.file) == head.ok)
  if (!head.chk) {
    message("Invalid header in CEL.list.file !")
    message(paste0("EXPECTED : ", head.ok))
    message(paste0("FOUND : ", colnames(myCELs)))
    stop("Invalid header.")
  }
  sn.chk <- duplicated(myCELs$SampleName)
  if (any(sn.chk)) {
    message("CEL.list.file contains duplicated SampleNames !")
    message(myCELs$SampleName[which(duplicated(myCELs$SampleName))])
    stop("Duplicated SampleNames.")
  }
  fecheck <- !vapply(myCELs$cel_files, file.exists, TRUE)
  fecheck.pos <- which(fecheck)
  if (length(fecheck.pos) > 0) stop(paste0("\n", "CEL file could not be found : ", myCELs$cel_files[fecheck.pos], collapse = ""))
  
  ## Adjusting cores/threads
  message("Adjusting number of threads if needed ...")
  avail.cores <- parallel::detectCores(logical = TRUE)
  if (is.null(nthread)) { nthread <- avail.cores -1; message(paste0("Reset nthread to ", nthread)) }
  if (nrow(myCELs) < nthread) { nthread <- nrow(myCELs); message(paste0("Reset nthread to ", nthread)) }
  if (avail.cores <= nthread) message(paste0(" WARNING : nthread set to ", nthread, " while available logical threads number is ", avail.cores, " !"))
  
  ## Building cluster
  `%dopar%` <- foreach::"%dopar%"
  cl <- parallel::makeCluster(spec = nthread, type = cluster.type, outfile = "")
  doParallel::registerDoParallel(cl)
  
  csres <- foreach::foreach(p = seq_len(nrow(myCELs)), .inorder = FALSE, .errorhandling = "stop") %dopar% {
    apt.cytoscan.process(CEL = myCELs$cel_files[p], samplename = myCELs$SampleName[p], ...)
  }
  
  ## Stopping cluster
  message("Stopping cluster ...")
  parallel::stopCluster(cl)
  
  message("Done.")
}

## Print thread-tagged message
tmsg <- function(text = NULL) { message(paste0(" [", Sys.info()[['nodename']], ":", Sys.getpid(), "] ", text)) }

## A more robust way to get machine OS type
get.os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  return(tolower(os))
}

