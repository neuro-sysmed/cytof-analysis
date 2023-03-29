
load_libraries <- function(force_install=FALSE) {

    if (force_install || !require("BiocManager", quietly = TRUE))
        install.packages("BiocManager", force=force_install)
    suppressMessages(suppressWarnings(library(BiocManager, quietly = TRUE)))

    if (force_install || !require("devtools", quietly = TRUE))
        install.packages("devtools", force=force_install)
    library(devtools, quietly = TRUE)

    if (force_install || !require("CATALYST", quietly = TRUE))
        install_github("neuro-sysmed/CATALYST", force=T)
    suppressPackageStartupMessages(library(CATALYST, quietly = TRUE))

    if (force_install || !require("flowCore", quietly = TRUE))
        BiocManager::install("flowCore", force=force_install)
    suppressMessages(suppressWarnings(library(flowCore, quietly = TRUE)))

    if (force_install || !require("stringi", quietly = TRUE))
        install.packages("stringi", force=force_install)
    library(stringi, quietly = TRUE)

    if (force_install || !require("readr", quietly = TRUE))
        install.packages("readr", force=force_install)
    library(readr, quietly = TRUE)


}


mk_outdirs <- function(output_dir ) {

 

    cat(fcs_dir, log_dir, "\n")

    return(list(fcs_dir, log_dir))

}

read_fcs_and_est_cutoffs <- function( directory ) {

    sce <- prepData( directory )
    sce <- assignPrelim(sce, sample_key)
#    cat(table(sce$bc_id))
    cat('ExtCutoffs...\n')
    sce <- estCutoffs(sce)

    cutoffs <- metadata(sce)$sep_cutoffs

#    table(sce$bc_id)


    # should be stored for later as well!
    cutoffs <- metadata(sce)$sep_cutoffs
    write.table(table(sce$bc_id), paste(log_dir, "barcode_counts.csv", sep="/"), sep="\t", col.names = NA, row.names = TRUE)
    write.table(cutoffs, paste(log_dir, "barcode_cutoffs.csv", sep="/"), sep="\t", col.names = NA, row.names = TRUE)


    plotYields(sce, which = c(0),
               out_path = directory,
               out_name = "barcode_separation")


    return(sce)
}

save_debarcoded_fcs_files <- function( sce, file_fstring ) {

    (fs <- sce2fcs(sce,split_by = "bc_id"))

    all(c(fsApply(fs, nrow)) == table(sce$bc_id))
    ids <- fsApply(fs, identifier)
    for (id in ids) {
        cat(sprintf("Writing %s barcode to %s file\n", id, sprintf(file_fstring, id)))
        ff <- fs[[id]]
        fn <- sprintf(file_fstring, id)
        write.FCS(ff, fn)
    }

}


plot_dna_scatter <- function( sce, outfile) {
    chs <- c("193Ir", "191Ir")
    jpeg(outfile)
    plotScatter(sce, chs)
    dev.off()

}


stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}


subsample_fcs <- function(fsc_file, max_events=5500, random_seed=10) {
    sce <- flowCore::read.FCS(fsc_file)
    set.seed(random_seed)

    # If there is less than max_events, select all
    picked = sample(seq_len(nrow(sce)),size = min(nrow(sce), max_events))

    sce_picked <- sce[ picked,]

    return(sce_picked)
}


write_fcs_file <- function(sce, filename) {
    write.FCS(sce, filename)
}
