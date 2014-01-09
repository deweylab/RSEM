repos <- biocinstallRepos()

test_biocinstallRepos_named_repositories <- function()
{

    allOS <- c("BioCsoft", "CRAN", "BioCann", "BioCexp", "BioCextra")
    windowsOnly <- "CRANextra"

    checkTrue(all(allOS %in% names(repos)))
    if (.Platform$OS.type == "windows")
    {
        checkTrue(windowsOnly %in% names(repos))
    } else {
        checkTrue(!windowsOnly %in% names(repos))
    }
   
}

test_biocinstallRepos_noNA_repositories <- function()
{
    checkTrue(!any(is.na(repos)))
}

test_biocinstallRepos_order <- function()
{
    checkIdentical("BioCsoft", names(repos)[[1]])
}
