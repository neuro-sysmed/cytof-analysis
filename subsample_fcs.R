library( flowCore )

source('cytof_utils.R')



data_dir <- 'data/debarcoded/'

file_pattern <- '.fcs' 
max_events <- 5500

subsample_dir <- sprintf("%s/subsampled_%d", data_dir, max_events)


for (file in list.files(data_dir, recursive = T, pattern=file_pattern)) {
    run_dir <- dirname(file)

    file <- sprintf("%s/%s", data_dir, file)
#    cat(sprintf("%s\n", basename(file)))
    output_dir <- sprintf("%s/%s", subsample_dir, run_dir)
#    cat(output_dir)
    outfile = sprintf("%s/%s", output_dir, basename(file))
#    cat(outfile)
#    break




    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)



    cat("Subsampling from:", file, "to", outfile, "\n")

    sce_sub <- subsample_fcs(file, max_events = max_events)
    write_fcs_file(sce_sub, outfile)

}

