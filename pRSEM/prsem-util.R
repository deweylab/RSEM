#
#  pliu 20150510
#
#  utility module for pRSEM
#

Util <- new.env()

Util$checkInstallCRAN <- function(pkg_name) {
  if ( ! pkg_name %in% rownames(installed.packages())){
    install.packages(pkg_name)
  }
}

Util$checkInstallBioc <- function(pkg_name) {
  if ( ! pkg_name %in% rownames(installed.packages() ) ){
    source("http://bioconductor.org/biocLite.R")
    biocLite(pkg_name)
  }
}

Util$checkIfGzipByExt <- function(file_name) {
  library(tools)
  ext <- file_ext(file_name)
  is_gz <- FALSE
  if ( ext %in% c('gz', 'gzip') ) {
    is_gz <- TRUE
  }
  return(is_gz)
}

Util$getFileNameSansExt <- function(infname) {
  fname <- basename(infname) 
  fname_sans_ext <- strsplit(fname, '.', fixed=T)[[1]][1]
  return(fname_sans_ext)
}
