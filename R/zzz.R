#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    # pkgVersion <- packageDescription(pkgname, fields="Version")
    # msg <- paste0(pkgname, " v", pkgVersion, "  ",
    #               "For help: https://yulab-smu.top/biomedical-knowledge-mining-book/", "\n\n")
    
    # Define a cache directory
    options(clusterProfiler_cache_dir = tempdir())

    packageStartupMessage(yulab.utils::yulab_msg(pkgname))
}


