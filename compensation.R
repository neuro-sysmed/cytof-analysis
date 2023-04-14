source('d:/analysis/cytof-analysis-main/cytof_utils.R')
library(cowplot)


load_libraries(F)
datafile   <- "D:/data/analysis/debarcoded/subsampled_55000/samples.fcs"
beadsfile  <-  "D:/data/analysis/Beads/Beads/" 
beadsfile  <- 
#170322_Processed_2.fcs"

bc_ms <- c(89, 102, 104:106, 108, 110, 112:114, 141:156, 158:176, 209)

sce <- prepData(beadsfile)
#sce <- assignPrelim(sce, sample_key, verbose = FALSE)
sce <- assignPrelim(sce, bc_ms, verbose = FALSE)
sce <- applyCutoffs(estCutoffs(sce))

# compute & extract spillover matrix

sce_spill <- computeSpillmat(sce)
sm <- metadata(sce_spill)$spillover_matrix

# values larger than 1
which(  sm > 1, arr.ind=T )
# manually set them to 1 for now.
sm[46,65] = 1


# do some sanity checks
chs <- channels(sce_spill)
ss_chs <- chs[rowData(sce_spill)$is_bc]
all(diag(sm[ss_chs, ss_chs]) == 1)
## [1] TRUE
all(sm >= 0 & sm <= 1)
## [1] TRUE

# drop the background and some empty channels
sm2 = sm[, colnames(sm) != "BCKG190Di"]

sm2 <- sm[, !colnames(sm) %in% c("BCKG190Di", "Cd111Di", "Cd116Di", "Ce140Di", "Gd157Di", "Pb208Di")]
sm3 <- sm2[, !colnames(sm2) %in% c("Pd102Di", "Pd104Di", "Pd105Di", "Pd106Di", "Pd108Di", "Pd110Di")]

#sm3 <- sm2[!rownames(sm2) %in% c("BCKG190Di", "Rh103Di","Cd111Di", "Cd116Di", "Ce140Di", "Gd157Di", "Pb208Di"), ] 
sm3 <- sm3[!rownames(sm2) %in% c("Pd102Di", "Pd104Di", "Pd105Di", "Pd106Di", "Pd108Di", "Pd110Di"),]

plotSpillmat(sce_spill, sm3) 

sce_comped <- compCytof(sce_spill, sm3, method = "nnls", overwrite = FALSE)


chs <- c("Dy162Di", "Gd160Di")
as <- c("exprs", "compexprs")
ps <- lapply(as, function(a) 
  plotScatter(sce_comped, chs, assay = a))
plot_grid(plotlist = ps, nrow = 1)


