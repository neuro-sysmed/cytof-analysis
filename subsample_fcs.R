library( flowCore )

source('cytof_utils.R')



data_dir <- 'data/debarcoded/'
file_pattern <- 'all.fcs' 
max_events = 5500




for (file in list.files(data_dir, recursive = T, pattern=file_pattern)) {
    file = sprintf("%s/%s", data_dir, file)

    output_dir = sprintf("%s/subsampled_%d/", dirname(dirname(file)), max_events)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


    outfile = sprintf("%s/%s", output_dir, basename(file))

    cat( "Subsampling from:", file, "to", outfile, "\n")

    sce_sub <- subsample_fcs(file)
    write_fcs_file(sce_sub, outfile)

}

