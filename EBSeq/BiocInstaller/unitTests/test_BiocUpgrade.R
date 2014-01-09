test_useDevel <- function()
{
    if (!BiocInstaller:::IS_END_OF_LIFE) {
        checkException(useDevel(), silent=TRUE)
    } else if (!BiocInstaller:::IS_DOWNGRADEABLE) {
        checkException(useDevel(FALSE), silent=TRUE)
    }
    if (!BiocInstaller:::IS_UPGRADEABLE) {
        checkException(useDevel(), silent=TRUE)
        opts <- options(warn=2); on.exit(options(opts))
        checkException(biocLite("BiocUpgrade"))
    }
}

test_getContribUrl_exist <- function()
{
    fun <- BiocInstaller:::.getContribUrl
    
    vers <- BiocInstaller:::BIOC_VERSION
    checkTrue(grepl(vers, fun(vers)))
    if (BiocInstaller:::IS_UPGRADEABLE) {
        vers <- BiocInstaller:::UPGRADE_VERSION
        checkTrue(grepl(vers, fun(vers)))
    }
    if (BiocInstaller:::IS_DOWNGRADEABLE) {
        vers <- BiocInstaller:::DOWNGRADE_VERSION
        checkTrue(grepl(vers, fun(vers)))
    }
}
