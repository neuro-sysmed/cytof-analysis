library(premessa)


files  <- commandArgs(trailingOnly = TRUE)
#cat( files )

outfile <- files[1]
infiles <- files[-1]


concatenate_fcs_files(infiles, outfile)