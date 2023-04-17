plotSpillmat <- function(x, sm = NULL, anno = TRUE, 
    isotope_list = CATALYST::isotope_list,
    hm_pal = c("white", "lightcoral", "red2", "darkred"), anno_col = "black") {
    
    # check validity of input spillover matrix
    chs <- colnames(sm)
    bc_chs <- chs[rowData(x)$is_bc]
    bc_chs <- chs[rowData(x[rowData(x)$is_bc])$is_bc]
    bc_idx <- which(bc_chs %in% colnames(sm))
    bc_rng <- seq(min(bc_idx), max(bc_idx))
    sm <- .make_symetric(sm)[bc_rng, bc_rng]
    
    df <- melt(sm * 100)
    
    colnames(df) <- c("emitting", "receiving", "spill")
    df$spillover <- paste0(sprintf("%2.3f", df$spill), "%")
    max <- ceiling(max(100 * sm[row(sm) != col(sm)]) / 0.25) * 0.25
    total <- paste0(sprintf("%2.2f", colSums(sm) * 100 - 100), "%")
    
    labs <- chs[bc_rng]
    ex <- !labs %in% chs[bc_idx]
    lab_cols <- rep("black", nrow(sm))
    lab_cols[ex] <- "grey"
    total[ex] <- NA
    
    if (hm_pal[1] != "white")
        hm_pal <- c("white", hm_pal)

    p <- ggplot(df, 
        aes_string("receiving", "emitting", group = "spillover")) +
        geom_tile(aes_string(fill = "spill"), col = "lightgrey") +
        scale_fill_gradientn(
            guide = FALSE, limits = c(0, max),
            colors = hm_pal, na.value = "lightgrey") +
        scale_x_discrete(limits = colnames(sm), expand = c(0, 0)) +
        scale_y_discrete(limits = rev(rownames(sm)), expand = c(0, 0)) +
        coord_fixed() + theme_bw() + theme(
            axis.title = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
            panel.border = element_blank(), panel.grid = element_blank())

    if (anno) {
        anno <- sprintf("%.1f", df$spill)
        anno[df$spill == 0 | df$spill == 100] <- ""
        p <- p + geom_text(
            aes_string(label = "anno"), 
            col = anno_col, size = 2)
    }

    return(p)
}
