#
#  pliu 20160911
#
#  install Bioconductor and pRSEM-required libraries to this directory
#
#  require R-3.3.1
#
#  CRAN: data.table, caTools
#  Local: spp
#  BioC v3.3: ShortRead, GenomicRanges
#
#  install devtools first and use its' install_version to install packages in
#  particular version.
#
#  packages repos are obtained from
#  source("http://bioconductor.org/biocLite.R")
#  biocinstallRepos()
#

main <- function() {
  param <- list(
    lib_loc = './',
    repos = list(
      BioCsoft  = "http://bioconductor.org/packages/3.3/bioc/",
      BioCann   = "http://bioconductor.org/packages/3.3/data/annotation/",
      BioCexp   = "http://bioconductor.org/packages/3.3/data/experiment/",
      BioCextra = "http://bioconductor.org/packages/3.3/extra/",
      CRAN      = "http://cran.us.r-project.org"
    ),

    pkg_spp  = '../phantompeakqualtools/spp_1.10.1_on_R3.3/',

    pkg2ver = list(
      ## name          version
      caTools       = '1.17.1',  ## for spp
      data.table    = '1.9.6',
      GenomicRanges = '1.24.3',
      ShortRead     = '1.30.0'
    )
  )

  options(repos=structure(c(CRAN=param$repos$CRAN)))
  installRLib(param)
}


installRLib <- function(param) {
  prsem_installed_pkgs <- rownames(installed.packages(lib.loc=param$lib_loc))

  if ( ! 'devtools' %in% prsem_installed_pkgs ) {
    install.packages('devtools', lib=param$lib_loc, type='source')
  }

  .libPaths(c(param$lib_loc, .libPaths()))
  library(devtools)

  for ( pkg_name in names(param$pkg2ver) ) {
    pkg_version <- param$pkg2ver[[pkg_name]]
    if ( ! pkg_name %in% prsem_installed_pkgs ) {
      install_version(pkg_name, version=pkg_version, repos=param$repos,
                      lib=param$lib_loc, type='source')
    }
  }

  if ( ! 'spp' %in% prsem_installed_pkgs ) {
    install.packages(param$pkg_spp, lib=param$lib_loc, repos=NULL,
                     type='source')
  }
}

main()
