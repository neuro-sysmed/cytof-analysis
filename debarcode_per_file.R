source('cytof_utils.R')


load_libraries(F)


runs_dir <- commandArgs(trailingOnly = TRUE)
runs_dir <- "d:/data/analysis/fcs"
runs_dir <- "data/runs"

base_dir <- dirname(runs_dir)

for (folder in list.dirs(runs_dir, recursive = FALSE)) {

    run_name <- basename( folder )
    output_dir <- sprintf("%s/debarcoded/%s/", base_dir, run_name)
    cat( sprintf("Processing folder %s --> run-name %s\n", folder, run_name) )
    cat( sprintf("Storing debarcoded files in %s\n", output_dir) )

    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    log_dir <- paste(output_dir, "logs",  sep="/")
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    fcs_dir <- paste(output_dir, "fcs",  sep="/")
    dir.create(fcs_dir, recursive = TRUE, showWarnings = FALSE)


    sce <- read_fcs_and_est_cutoffs( folder )
    cutoffs <- metadata(sce)$sep_cutoffs
    cat('Plotting\n')
    plot_dna_scatter(sce, paste(log_dir, "DNA_scatter.jpg", sep='/'))

    sce_specific <- applyCutoffs(sce)
    sce_global   <- applyCutoffs(sce, sep_cutoffs = mean(cutoffs))
    sce_all      <- applyCutoffs(sce, sep_cutoffs = 0)

    seps <- data.frame (cutoff    = c(specific = -1,
                                      global   = mean(cutoffs),
                                      all      = 0),
                        yield     = c(specific = mean(sce_specific$bc_id != 0),
                                      global   = mean(sce_global$bc_id != 0),
                                      all      = mean(sce_all$bc_id != 0))
                    )

    write.table(seps, paste(log_dir, "barcode_separation.csv", sep="/"), sep="\t", col.names = NA, row.names = TRUE)

    save_debarcoded_fcs_files(sce, paste(fcs_dir, "/%s_specific.fcs", sep="/"))
}


cat('All done for now...\n')

stop_quietly()



bcodes = row.names( sample_key )

for (filename in files ) {
    write(sprintf("handling file: %s", filename), stdout())
    file_sce <- prepData( filename, transform=TRUE )

    if (metadata(file_sce)$experiment_info$n_cells == 0) {
        write(sprintf("No data in %s, skipping", filename), stdout())
        next
    }

    file_sce <- applyCutoffs(file_sce, sep_cutoffs=cutoffs)

    (fs <- sce2fcs(file_sce, split_by = "bc_id"))
    all(c(fsApply(fs, nrow)) == table(file_sce$bc_id))
    ids <- fsApply(fs, identifier)
    for (id in ids) {

        dir.create( sprintf("%s/%s", tmp_dir, id), recursive = TRUE)

        ff <- fs[[id]]                     # subset 'flowFrame'
        fn <- sprintf("%s/%s/%s", tmp_dir, id, basename(filename))
#        fn <- sprintf("processed/%s/fcs/%s_specific.fcs", run_dir, id)
        write(sprintf("writing barode %s to file %s", id, fn), stdout())
        write.FCS(ff, fn)                  # write frame to FCS
    }

#    break
}

# Now concatenate the files into a single file per barcode.

for (barcode in bcodes) {
    write(sprintf('Merging files for barcode %s', barcode), stdout())
    bc_files <- list.files(path = sprintf("%s/%s", tmp_dir, barcode),full.names = TRUE, pattern = ".fcs")
    outfile <- sprintf("%s/%s.fcs", fcs_dir,barcode)
    premessa::concatenate_fcs_files(bc_files, output.file = outfile)
}







#for (bcode in bcodes ) {
#    plotEvents(sce, which = c(0, bcode), n = 25,
#               out_path = img_dir,
#               out_name = paste("events_", bcode, sep=''))
#
#    pdf_file <- paste(img_dir, "/", sep='')
#    pdf_file <- paste(pdf_file, "mahal_", bcode, ".pdf", sep='')
#    print(pdf_file)
#
#    pdf(pdf_file)
#    plotMahal(sce, which = bcode)
#    dev.off()
#}




(fs <- sce2fcs(sce_global, split_by = "bc_id"))
(c(fsApply(fs, nrow)) == table(sce_global$bc_id))
ids <- fsApply(fs, identifier)
for (id in ids) {
    ff <- fs[[id]]                     # subset 'flowFrame'
    fn <- sprintf("processed/%s/fcs/%s_global.fcs", run_dir, id)

#    fn <- file.path(outputdir, fn)
#    flowCore:::updateTransformKeywords(ff)
    write.FCS(ff, fn)                  # write frame to FCS
}

(fs <- sce2fcs(sce_specific, split_by = "bc_id"))
all(c(fsApply(fs, nrow)) == table(sce_specific$bc_id))
ids <- fsApply(fs, identifier)
for (id in ids) {
    ff <- fs[[id]]                     # subset 'flowFrame'
    fn <- sprintf("processed/%s/fcs/%s_specific.fcs", run_dir, id)
#    fn <- file.path(outputdir, fn)
    write.FCS(ff, fn)                  # write frame to FCS
}
