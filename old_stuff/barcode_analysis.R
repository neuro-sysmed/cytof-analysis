library(CATALYST)
#library(flowCore)
#library(readr)

run_dir <- commandArgs(trailingOnly=TRUE)

# create directories to store the various data in
#for (process_dir in c('global','specific', 'all')) {
##  for (res_dir in c('fcs','images','output')) {
#    output_dir <- paste("processed", run_dir, process_dir, sep="/")
#    dir.create(output_dir,recursive = TRUE)
##  }
#}

img_dir <- paste("processed",run_dir, "images", sep="/")
dir.create(img_dir, recursive = TRUE)
log_dir <- paste("processed", run_dir, "logs",  sep="/")
dir.create(log_dir, recursive = TRUE)

fcs_dir <- paste("processed", run_dir, "fcs",  sep="/")
dir.create(fcs_dir, recursive = TRUE)

tmp_dir <- paste("tmp",run_dir, "images", sep="/")
dir.create(tmp_dir, recursive = TRUE)


raw_dir <- paste("raw", run_dir, sep="/")


files <- list.files(path = raw_dir,full.names = TRUE, pattern = ".fcs")

sce <- prepData( raw_dir, transform=TRUE )
# should be written to a log file! Nr of data points per raw fcs file
#table(sce$sample_id)
#names(int_colData(sce))

# debarcoding! sample_key comes from catalyst
sce <- assignPrelim(sce, sample_key)
#rowData(sce)$is_bc
#rownames(sce)[rowData(sce)$is_bc]
# Should be written to a log file! Nr of data points per barcode!
table(sce$bc_id)
write.table(table(sce$bc_id), paste(log_dir, "barcode_counts.csv", sep="/"), sep="\t", col.names = NA, row.names = TRUE)

sce <- estCutoffs(sce)
# should be stored for later as well!
cutoffs <- metadata(sce)$sep_cutoffs
write.table(cutoffs, paste(log_dir, "barcode_cutoffs.csv", sep="/"), sep="\t", col.names = NA, row.names = TRUE)


plotYields(sce, which = c(0),
           out_path = img_dir,
           out_name = "barcode_separation")

bcodes = row.names( sample_key )



for (bcode in bcodes ) {
  # Filename be changed to reflect the name instead of the barcode,
  # as it oddly does not look to be possible within the function.
  #
  plotYields(sce, which = c(bcode),
             out_path = img_dir,
             out_name = paste("barcode_", bcode, sep=''))
}

for (file in files ) {
    file_sce <- prepData( file, transform=TRUE )
    file_sce <- applyCutoffs(file_sce, sep_cutoffs=cutoffs)

    fn <- sprintf("%s/%s", tmp_dir, basename(files[1]))

    write.FCS(ff, fn)                  # write frame to FCS


}


sce_specific <- applyCutoffs(sce)
sce_global   <- applyCutoffs(sce, sep_cutoffs = mean(cutoffs))
sce_all      <- applyCutoffs(sce, sep_cutoffs = 0)

seps <- data.frame (cutoff    = c(specific = -1,
                                global     = mean(cutoffs),
                                all        = 0),
                    yield = c(specific = mean(sce_specific$bc_id != 0),
                                global     = mean(sce_global$bc_id != 0),
                                all        = mean(sce_all$bc_id != 0))
                  )

write.table(seps, paste(log_dir, "barcode_separation.csv", sep="/"), sep="\t", col.names = NA, row.names = TRUE)





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


chs <- c("193Ir", "191Ir")
chs <- c("Time", "191Ir")
jpg_file <- paste(img_dir, "DNA_scatter_all.jpg", sep='/')
jpeg(jpg_file)
plotScatter(sce_all, chs)
dev.off()

jpg_file <- paste(img_dir, "DNA_scatter_global.jpg", sep='/')
jpeg(jpg_file)
plotScatter(sce_global, chs)
dev.off()

jpg_file <- paste(img_dir, "DNA_scatter_specific.jpg", sep='/')
jpeg(jpg_file)
plotScatter(sce_specific, chs)
dev.off()



# chs <- c('CD3', 'CD4', 'CD8a', 'CD90')
# #chs <- grep("^CD", rownames(sce), value = TRUE)
# for (bcode in bcodes ) {
# 
#   pdf_file <- paste(img_dir, "/", sep='')
#   pdf_file <- paste(pdf_file, "cell_separation_", bcode, ".pdf", sep='')
#   pdf(pdf_file)
#   p <- plotScatter(sce, chs)
#   p$facet$params$ncol <- 3;
#   p;
#   dev.off()
# }

#ids <- sample(rownames(sample_key), 3)
#ids
#sce <- sce[chs, sce$bc_id %in% ids]
#plotScatter(sce, chs, color_by = "bc_id")
#sce <- sce2

(fs <- sce2fcs(sce_all, split_by = "bc_id"))
all(c(fsApply(fs, nrow)) == table(sce_all$bc_id))
ids <- fsApply(fs, identifier)
for (id in ids) {
    ff <- fs[[id]]                     # subset 'flowFrame'
    fn <- sprintf("processed/%s/fcs/%s_all.fcs", run_dir, id)
#    fn <- file.path(outputdir, fn)
    write.FCS(ff, fn)                  # write frame to FCS
}

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
    write(id, stdout())
    ff <- fs[[id]]                     # subset 'flowFrame'
    fn <- sprintf("processed/%s/fcs/%s_specific.fcs", run_dir, id)
#    fn <- file.path(outputdir, fn)
    write.FCS(ff, fn)                  # write frame to FCS
}
