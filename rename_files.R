ss <- read.delim('Barcode_sample_info.csv')

data_dir <- "data_full/debarcoded/subsampled_5500/"

for(j in 1:nrow(ss)) {
  for(k in 2:ncol(ss))  {

    run_name = colnames(ss)[k]
    file_name = ss[j,1]
    sample_name = ss[j,k]

    exts <- c('all', 'global', 'specific')

    for (ext in exts) {
        ff <- read.FCS( sprintf("%s/%s/fcs/%s_%s.fcs", data_dir, run_name, file_name, ext) )
        identifier(ff) <- sample_name
        write.FCS(ff, sprintf("%s/%s_%s.fcs", data_dir, sample_name, ext))
    }
  }
}
