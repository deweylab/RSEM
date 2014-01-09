## Mirrors: uncomment the following and change to your favorite CRAN mirror
## if you don't want to use the default (cran.fhcrc.org, Seattle, USA).
## options("repos" = c(CRAN="http://cran.fhcrc.org"))

## Mirrors: uncomment the following and change to your favorite Bioconductor
## mirror, if you don't want to use the default (www.bioconductor.org,
## Seattle, USA)
## options("BioC_mirror" = "http://www.bioconductor.org")

local({
    currBiocVers <- 
        package_version(readLines("http://bioconductor.org/bioc-version",
        warn=FALSE))
    vers <- getRversion()
    biocVers <- tryCatch({
        BiocInstaller::biocVersion() # recent BiocInstaller
    }, error=function(...) {         # no / older BiocInstaller
        tools:::.BioC_version_associated_with_R_version
    })

    if (biocVers < currBiocVers) {
        txt <- strwrap(sprintf("Your Bioconductor is out-of-date, upgrade
            to version %s by following instructions at
            http://bioconductor.org/install.", currBiocVers))
        message(paste(txt, collapse="\n"))
    }

    if (vers > "2.13" && biocVers > "2.8") {

        if (exists("biocLite", .GlobalEnv, inherits=FALSE)) {
            txt <- strwrap("There is an outdated biocLite() function in the
                global environment; run 'rm(biocLite)' and try again.")
            stop("\n", paste(txt, collapse="\n"))
        }
          
        if (!suppressWarnings(require("BiocInstaller", quietly=TRUE))) {
            a <- NULL
            p <- file.path(Sys.getenv("HOME"), ".R", "repositories")
            if (file.exists(p)) {
                a <- tools:::.read_repositories(p)
                if (!"BioCsoft" %in% rownames(a)) 
                    a <- NULL
            }
            if (is.null(a)) {
                p <- file.path(R.home("etc"), "repositories")
                a <- tools:::.read_repositories(p)
            }
            if (!"package:utils" %in% search()) {
                url <- "http://bioconductor.org/biocLite.R"
                txt <- sprintf("use 'source(\"%s\")' to update 'BiocInstaller'
                                after 'utils' package is attached",
                               url)
                message(paste(strwrap(txt), collapse="\n  "))
            } else {
                ## add a conditional for Bioc releases occuring WITHIN
                ## a single R minor version
                if (vers >= "2.15" && vers < "2.16") {
                    a["BioCsoft", "URL"] <- sub(as.character(biocVers), "2.11",
                      a["BioCsoft", "URL"])
                    biocVers <- numeric_version("2.11")
                }
                install.packages("BiocInstaller", repos=a["BioCsoft", "URL"])
                if (!suppressWarnings(require("BiocInstaller",
                                              quietly=TRUE))) {
                    url0 <- "http://www.bioconductor.org/packages"
                    url <- sprintf("%s/%s/bioc",
                                   url0, as.character(biocVers))
                    txt0 <- "'biocLite.R' failed to install 'BiocInstaller',
                            use 'install.packages(\"%s\", repos=\"%s\")'"
                    txt <- sprintf(txt0, "BiocInstaller", url)
                    message(paste(strwrap(txt), collapse="\n  "))
                }
            }
        }
    } else {
        source("http://bioconductor.org/getBioC.R")
        biocLite <<-
            function(pkgs, groupName="lite", ...)
            {
                if (missing(pkgs))
                    biocinstall(groupName=groupName, ...)
                else
                    biocinstall(pkgs=pkgs, groupName=groupName, ...)
            }
    }
})
