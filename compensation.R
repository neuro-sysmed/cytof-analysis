source('d:/analysis/cytof-analysis-main/cytof_utils.R')
library(cowplot)


source('CATALYST_overloads.R')

load_libraries(F)
datafile   <- "D:/data/analysis/debarcoded/subsampled_55000/samples.fcs"
datafile   <- 'raw/all_5500.fcs'

beadsfile  <-  "D:/data/analysis/Beads/Beads/" 
beadsfile  <- "raw/Beads170322_Processed_2.fcs"
#170322_Processed_2.fcs"

bc_ms <- c(89, 102, 104:106, 108, 110, 112:114, 141:156, 158:176, 209)

sce_beads <- prepData(beadsfile)
#sce <- assignPrelim(sce, sample_key, verbose = FALSE)
sce_beads <- assignPrelim(sce_beads, bc_ms, verbose = FALSE)
sce_beads <- applyCutoffs(estCutoffs(sce_beads))

# compute & extract spillover matrix

sce_spill <- computeSpillmat(sce_beads)
chs <- channels(sce_spill)
ss_chs <- chs[rowData(sce_spill)$is_bc]
sm <- metadata(sce_spill)$spillover_matrix


# Only use the channels with antibodies in them
sm_bc <- sm[,  colnames(sm) %in% ss_chs]

# do some sanity checks
all(diag(sm_bc[ss_chs, ss_chs]) == 1)
## [1] TRUE
all(sm_bc >= 0 & sm_bc <= 1)
## [1] TRUE

plotSpillmat(sce_spill, sm_bc) 

sce_data <- prepData(datafile)


sce_comped <- compCytof(sce_data, sm_bc, method = "nnls", overwrite = FALSE)


chs <- c("Dy162Di", "Gd160Di")
as <- c("exprs", "compexprs")
ps <- lapply(as, function(a) 
  plotScatter(sce_comped, chs, assay = a))
plot_grid(plotlist = ps, nrow = 1)

fcs_comped <- sce2fcs(sce_comped)
write_fcs_file(fcs_comped, 'all_compensated.fcs')


