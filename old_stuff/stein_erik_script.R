# Libraries
library(CATALYST)
#BiocManager::install("CATALYST")
library(flowCore)
#BiocManager::install("flowCore")
library(stringi)
#install.packages("stringi")
library(cowplot)
#install.packages("cowplot")
library(devtools)
#install.packages("devtools")
library(premessa)
#install_github("ParkerICI/premessa")
library(readr)



# Make a new Project in R somewhere on your computher that you can easily get to; documents. 
# In the folder with the R project file (.Rproj), make a folder called data and put your data in 
# that folder. To avoid many errors introduced by various algorithms, please use data that has 
# note been through other algorithms; directly from the XT (_prosessed.fcs). However, doing 
# gaussian gating in cytobank may be a good idea. If you do, please dont use a very strict DNA vs
# time gate, its best to leave, most of the doublet removal to the debarcoding. 

# Find data
###########
data_dir <- '~/Dropbox/Rstudio/Friday Meetings/data' 
files <- list.files(path = data_dir,full.names = TRUE, pattern = ".fcs")
files # Check that you have some file paths

# Look for double delimiters and fix
####################################
ff <- read.FCS(files[1], emptyValue= FALSE) #read the first file
stri_detect_fixed(unlist(ff@description), "\\") #check the whole header for double delims

# If you have, use the code below. I the code below I had double delims at ff@description$FolderName
for (file in files){
  print(file)
  ff <- read.FCS(file, emptyValue= FALSE)
  file_name <- file
  ff@description$FolderName = "file"
  new_file_name <- gsub(".fcs", "_edited.fcs", file_name)
  write.FCS(ff, filename = new_file_name, what="numeric", delimiter = "\\")
}

# Put all edit files in a separate folder!


# Fix channel names and remove the BCKG channel
###############################################
#remove all special signs from the channel name: - , . () [] ? % 
premessa::paneleditor_GUI() 
# Show Premessa to the folder with the edited files. 
# Edt the panels as neccecary and click prosess; the new files with edited pabels will 
# be made in a new folder. When it finishes click Dismiss and close the browser window. 
# In the consol press esc until a > appears again


# Concatenate data
##################
data_dir <- "~/Dropbox/Rstudio/Friday Meetings/data/renamed"
files <- list.files(path = data_dir,full.names = TRUE, pattern = ".fcs")
premessa::concatenate_fcs_files(files, output.file = "data/Concact_.fcs")

# This will make a concatenated file .fcs in the data folder. The name of the file 
# will here be Concact_.fcs. You might want to change that with some something more 
# meaningful or with better spelling.


# Debarcode
###########
file <- "~/Dropbox/Rstudio/Friday Meetings/data/Concact_.fcs"
sce <- prepData(file, transform = TRUE) #reads the data in sce
table(sce$sample_id)
key <- read.csv2('Bckey.csv', # Read debarcoding key
                 check.names=FALSE, 
                 row.names=1, 
                 header = TRUE)

key # We need integer data here. 

sce <- assignPrelim(sce, key) # Does preliminary debarcoding 
rownames(sce)[rowData(sce)$is_bc]
table(sce$bc_id)

sce <- estCutoffs(sce) # Estimate the threshold values for each sample. 
metadata(sce)$sep_cutoffs
plotYields(sce)
plotYields(sce, out_path = getwd(), out_name = "Allsamples")
plotYields(sce, which = rownames(key), out_path = getwd(), out_name = "Eachsample")
table(sce$sample_id)

sce <- applyCutoffs(sce) # Apply estimated thresholds

#Save debarcoded data
fs <- sce2fcs(sce, split_by = "bc_id") 
write.flowSet(fs, outdir = paste0(data_dir,"/Debarcoded FCS"))

#Compensate data
################

# Remove all channels that is not used to stain the beads (without an antibody name). 
#Leave in event length, time, and the gaussian channels. 
premessa::paneleditor_GUI() 
comp_file <- "~/Dropbox/Rstudio/Fridaymeeting2/renamed/210426_007_Bend of day CBeads_03_1.fcs" 

sce_comp <- prepData(comp_file, transform = TRUE) #read the comp beads

# Debarcode the single stained comp beads
bc_ms <- c(89, 111:114, 116, 141:176, 209)
sce_comp <- assignPrelim(sce_comp, bc_ms, verbose = FALSE)
sce_comp <- applyCutoffs(estCutoffs(sce_comp))

# compute & extract spillover matrix
sce_comp <- computeSpillmat(sce_comp)
sm <- metadata(sce_comp)$spillover_matrix
sm
# do some sanity checks
chs <- channels(sce_comp)
ss_chs <- chs[rowData(sce_comp)$is_bc]
all(diag(sm[ss_chs, ss_chs]) == 1)
all(sm >= 0 & sm <= 1)
plotSpillmat(sce_comp)

# Read the debarcoded data using a meta file. Use the meta.csv as a template.
meta <- read_delim("meta.csv", delim = ";")
meta$sample_id <- paste0(meta$condition,meta$patient_id)
sce <- prepData(meta$filename, md = meta)

# compensate using NNLS-method; keep uncompensated data
sce2 <- compCytof(sce, sm, method = "nnls", overwrite = FALSE)
chs <- c("Gd157Di", "Gd158Di")
as <- c("exprs", "compexprs")
ps <- lapply(as, function(a) 
  plotScatter(sce2, chs, assay = a))
plot_grid(plotlist = ps, nrow = 1)
rm(sce2)

# Overwrite to save;
sce <- compCytof(sce, sm, method = "nnls", overwrite = TRUE)

fs <- sce2fcs(sce, split_by = 'sample_id')
write.flowSet(fs, outdir = paste0(data_dir,"/comp"))


# Read in debarcoded data
#########################

source("~/Dropbox/Rstudio/Friday Meetings/utilityFunc.R")

# Read the first file to look at the channel names
pathtosamples <- "~/Dropbox/Rstudio/IngaDBC/data/Edited for 190/renamed/Debarcoded FCS/"
frame <- PanelToCSV(pathtosamples =  pathtosamples, nameoffile = "panel_FM.csv")

FuncInds <- c(53,	52,	29,	59,	43,	49,	51,	40,	38,	39,	54,	34,	37,	46,	41) # Functional markers
ToUse <- c(71,	30,	61,	17,	48,	47,	35,	56,	62,	36,	16,	45,	18,	42,	28,	32,	8,	31,	50,	33,	19,	63,	55,	44,	58,	57,	21,	60) # Clustering/DR

pData(parameters(frame))[,2][FuncInds] 
pData(parameters(frame))[,2][ToUse]

panel <- constPanel(frame, FuncInds, ToUse) # make panel and meta for prepData() function

# Read the meta file
meta <- read_delim("meta.csv", delim = ";")

# Threshold subsampling (optional)
set.seed(0) 
sampling.ceiling <- 10000 # max number of cells

fullframe <- read.flowSet(meta$filename, transformation = FALSE, truncate_max_range = FALSE)
fs.ds <- fsApply(fullframe, function(ff) {
  idx <- sample.int(nrow(ff), min(sampling.ceiling, nrow(ff)))
  ff[idx,]  # alt. ff[order(idx),]
})

sce <– prepData(fs.ds, panel, meta, md_cols = list(factors = colnames(meta)))
table(sce$sample_id)
sce

plotCounts(sce, 
           group_by = "sample_id", 
           color_by = "condition")

rm(fs.ds,fullframe,panel) 


# FlowSOM ands UMAP
###################
set.seed(1601)
sce <- cluster(sce, features = "type", 
               xdim = 10, ydim = 10, maxK = 30, 
               verbose = FALSE, seed = 1) 
sce <- runDR(sce, dr = "UMAP", cells = 5000, features = "type") # cells is per sample!

# Save/Open SCE
saveRDS(sce, file ="IngaDBC.RData") # Save the data for late use

Plotter(sce, name="Plotter", toremove = NULL)
clustPlotter(sce, name="Plotter", clust = "meta30")
             
ncells <- table(sce$sample_id, cluster_ids(sce, k = "mPARC"))