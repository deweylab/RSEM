#
#  pliu 20150618
#
#  install Bioconductor and pRSEM-required libraries to this directory
#
#  require R-3.2.0,
#
#  CRAN: data.table
#  Local: spp
#  BioC: ShortRead, GenomicRanges
#

main <- function() {
  param <- list(
    lib_loc = './',
    typ2urlprefix = list(
      spp  = '../phantompeakqualtools/',
      cran = 'http://cran.r-project.org/src/contrib/',
      bioc = 'http://www.bioconductor.org/packages/release/bioc/src/contrib/'
    ),

    typ2pkgs = list(
      spp  = c( 'spp_1.10.1' ),

      cran = c( 'caTools_1.17.1',
                'futile.logger_1.4.1',
                'futile.options_1.0.0',
                'snow_0.3-13',
                'lambda.r_1.1.7',
                'bitops_1.0-6',
                'hwriter_1.3.2',
                'latticeExtra_0.6-26',
                'RColorBrewer_1.1-2',
                'reshape2_1.4.1',
                'plyr_1.8.3',
                'Rcpp_0.11.6',
                'stringr_1.0.0',
                'stringi_0.4-1',
                'magrittr_1.5',
                'chron_2.3-46',
                'data.table_1.9.4' 
              ),

      bioc = c( ## do not change package order
                'BiocGenerics_0.14.0',
                'BiocParallel_1.2.3',
                'Biostrings_2.36.1',
                'Biobase_2.28.0',
                'S4Vectors_0.6.0',
                'IRanges_2.2.4',
                'XVector_0.8.0',
                'zlibbioc_1.14.0',

                'Rsamtools_1.20.4',
                'GenomeInfoDb_1.4.1',

                'GenomicRanges_1.20.5',
                'GenomicAlignments_1.4.1',
                'ShortRead_1.26.0'
              )
    )
  )

  prsem_installed_pkgs <- rownames(installed.packages(lib.loc=param$lib_loc))

  for ( typ in names(param$typ2pkgs) ) {
    urlprefix <- param$typ2url[[typ]]
    for ( pkg_name in param$typ2pkgs[[typ]] ) {
      words <- strsplit(pkg_name, '_')[[1]]
      if ( ! words[1] %in% prsem_installed_pkgs ) {
        pkg_url <- paste0(urlprefix, pkg_name, '.tar.gz')
        install.packages(pkg_url, lib=param$lib_loc, repos=NULL, type='source')
      }
    }
  }
}

main()
