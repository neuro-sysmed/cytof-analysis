
load_libraries <- function() {

    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    library(BiocManager)

    if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
    library(devtools)

    if (!require("CATALYST", quietly = TRUE))
        install_github("neuro-sysmed/CATALYST")
    library(CATALYST)

    if (!require("flowCore", quietly = TRUE))
    BiocManager::install("flowCore")
    library(flowCore)

    if (!require("stringi", quietly = TRUE))
    install.packages("stringi")
    library(stringi)

    if (!require("readr", quietly = TRUE))
    install.packages("readr")
    library(readr)


}

