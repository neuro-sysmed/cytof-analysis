library(premessa)
library( stringr )


file_pattern <- 'all.fcs' 
events <- c(5500, 55000)
input_dir <- 'data/debarcoded/subsampled_'

for (event in events) {

    data_dir <- sprintf('%s%d/', input_dir, event)
    files <- paste(data_dir, list.files(data_dir, recursive = T, pattern=file_pattern), sep="")
    cat(sprintf("Processing %s\n", data_dir))

#    cat(files )

    ctrl_files  <- files[ str_detect(files, "C5")]
    smple_files <- files[! str_detect(files, "C5")]

#    cat(ctrl_files)
#    cat(smple_files)
    concatenate_fcs_files(files, sprintf("%s/all.fcs",data_dir))
    concatenate_fcs_files(smple_files, sprintf("%s/samples.fcs",data_dir))
    concatenate_fcs_files(ctrl_files, sprintf("%s/ctrls.fcs",data_dir))

}



#outfile <- files[1]
#infiles <- files[-1]


#