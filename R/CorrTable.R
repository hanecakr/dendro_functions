CorrTable <-
        function(x,
                 y = NULL,
                 overlap = 50,
                 remove.duplicates = TRUE,
                 output = "table" #c("matrix", "table")
                 ) {
        ### checks
                if (is.null(y)) {
                        y = x
                        noRef = TRUE}
                else {
                        noRef = FALSE
                }
                if (!("rwl" %in% class(x))) {
                        warning("'x' is not class rwl")
                }
                if (!all(diff(as.numeric(row.names(x))) == 1)) {
                        stop(
                        "The tree-ring series 'x' have/has no consecutive years in increasing order as rownames"
                        )
                }
                if (!all(diff(as.numeric(row.names(y))) == 1)) {
                        stop(
                        "The master series 'y' have/has no consecutive years in increasing order as rownames"
                        )
                }
                if (any(length(overlap) != 1 |
                        !is.numeric(overlap) |
                        overlap %% 1 != 0 |
                        overlap < 3)) {
                        stop("'overlap' should be a single integer >=3")
                }
                if (overlap < 50) {
                        warning(
                        "The minimum number of overlap is lower than 50. This might lead to statistically insignificant matches."
                        )
                }

        ### reduce x and y to common overlap, based on rownames = years
                interval_x <- rownames(x)
                interval_y <- rownames(y)
                interval <- intersect(interval_x, interval_y)
                x <- x[rownames(x) %in% interval, ]
                y <- y[rownames(y) %in% interval, ]
                
        ### compute overlap
                n <- dim(x)[2]
                m <- dim(y)[2]
                overlap_i <- matrix(NA_real_, nrow = n, ncol = m)
                rownames(overlap_i) <- names(x)
                colnames(overlap_i) <- names(y)
                
                for (i in 1:n) {
                        # substract single column from each matrix column => NA when no overlap!
                        OVL <- x[, i] - y
                        OVL <- colSums(!is.na(OVL))
                        OVL[OVL==0] <- NA
                        OVL[OVL<overlap] <- NA
                        overlap_i[i, ] <- OVL
                }
                
        ### parallel variation (%PV | GLK)
                GLK_mat <- matrix(NA_real_, nrow = n, ncol = m)
                rownames(GLK_mat) <- names(x)
                colnames(GLK_mat) <- names(y)
                treering_sign_x <- apply(x, 2, diff)
                treering_sign_x <- sign(treering_sign_x)
                treering_sign_y <- apply(y, 2, diff)
                treering_sign_y <- sign(treering_sign_y)
                for (i in 1:n) {
                        # substract single column from each matrix column => NA when no overlap!
                        treering_GC <- abs(treering_sign_x[, i] - treering_sign_y)
                        GLK_values <-
                                1 - (colSums(treering_GC, na.rm = TRUE) / (2 * (overlap_i[i, ]-1)))
                        GLK_mat[i,] <- GLK_values
                }
                
                if (noRef == TRUE) {
                        diag(GLK_mat) <- 1
                }
        ### probability associated with %PV | GLK
                s_df <- 1 / (2 * sqrt(overlap_i))
                z_df <- (GLK_mat - .5) / s_df
                z_normcdf <-
                        apply(z_df, 2, function(z)
                                pnorm(z, mean = 0, sd = 1))
                GLK_p <- 2 * (1 - z_normcdf)

        ### compute t-values according to the Hollstein 1980 algorithm
                tHo_mat <- matrix(NA_real_, nrow = n, ncol = m)
                rownames(tHo_mat) <- names(x)
                colnames(tHo_mat) <- names(y)
                wuch_x <- apply(x, 2, function(x) {x/lag(x)})
                wuch_x <- 100*log10(wuch_x)
                if (noRef == FALSE){
                        wuch_y <- apply(y, 2, function(x) {x/lag(x)})
                        wuch_y <- 100*log10(wuch_y)
                        r <- cor(wuch_x, wuch_y, method = "pearson", use = "pairwise.complete.obs")
                        
                } else {
                        r <- cor(wuch_x, wuch_x, method = "pearson", use = "pairwise.complete.obs")
                }
                # if r <0, r is set to zero
                r[r<0] <- 0
                tHo_mat <- round(r * sqrt(overlap_i-2)/sqrt(1-r^2), 2)
        
        ### compute t-values according to the Baillie-Pilcher algotihm
                tBP_mat <- matrix(NA_real_, nrow = n, ncol = m)
                rownames(tBP_mat) <- names(x)
                colnames(tBP_mat) <- names(y)

                movav5_x <- apply(x, 2, function(x) {MovAv(x, w = 5)})
                rownames(movav5_x) <- rownames(x)
                movav5_x <- 100*x/movav5_x
                movav5_x <- log(movav5_x)
                if (noRef == FALSE){
                        movav5_y <- apply(y, 2, function(x) {MovAv(x, w = 5)})
                        rownames(movav5_y) <- rownames(y)
                        movav5_y <- 100*y/movav5_y
                        movav5_y <- log(movav5_y)
                        # calculate Pearson r between pairs of series for common overlap
                        r <- cor(movav5_x, movav5_y, method = "pearson", use = "pairwise.complete.obs")
                } else {
                        r <- cor(movav5_x, movav5_x, method = "pearson", use = "pairwise.complete.obs")
                        
                }
                # if r <0, r is set to zero
                r[r<0] <- 0
                tBP_mat <- round(r * sqrt(overlap_i-2)/sqrt(1-r^2), 2)   
        
        ### r-Pearson
                r_pearson <- cor(x, y, method = "pearson", use = "pairwise.complete.obs")
                overlap_r <- is.na(overlap_i)
                overlap_r <- ifelse(overlap_r == TRUE, NA, 1)
                r_pearson <- r_pearson * overlap_r
                
        ### output
        corr_table <-
                list(overlap = overlap_i,
                     glk = round(100 * GLK_mat, 1),
                     glk_p = GLK_p,
                     tHo = tHo_mat,
                     tBP = tBP_mat,
                     r_pearson = r_pearson)

        if (output == "table"){
                
                descr <- dplR::rwl.stats(x) %>% 
                        select(series, first, last, year)
                
                if(noRef == FALSE){
                descr_y <- dplR::rwl.stats(y) %>%
                        select(series, first, last, year)
                descr <- rbind(descr, descr_y)
                }
                
                for (i in names(corr_table)){
                        if (i == names(corr_table)[1]){
                                corr_table_i <- corr_table[[i]] %>%
                                        as.data.frame() %>%
                                        rownames_to_column("series") %>%
                                        pivot_longer(., cols = !series, names_to = "reference", values_to = i)
                        } else {
                                corr_table_i <- corr_table[[i]] %>%
                                        as.data.frame() %>%
                                        rownames_to_column("series") %>%
                                        pivot_longer(., cols = !series, names_to = "reference", values_to = i) %>% 
                                        left_join(corr_table_i, ., by = c("series", "reference"))
                        }
                }
        corr_table <-
                corr_table_i %>%
                left_join(descr %>% select(-year), by = c("reference" = "series")) %>%
                rename(ref_end = last,
                       ref_start = first) %>% 
                left_join(descr %>% select(-first), by = "series") %>%
                rename(series_length = year,
                       series_end = last)

        col_order <- c("series",
                       "series_length",
                       "series_end",
                       "reference",
                       "ref_start",
                       "ref_end",
                       "overlap",
                       "glk",
                       "glk_p",
                       "tBP",
                       "tHo",
                       "r_pearson")
        corr_table <- corr_table[, col_order]
        corr_table <- corr_table[complete.cases(corr_table[, 6:ncol(corr_table)]), ]

        ### remove duplicates
        if (remove.duplicates == TRUE && noRef == TRUE){
                corr_table <- subset(corr_table, series < reference)
        }

        ### sort by series, then by tHo
        corr_table <-
                corr_table[order(corr_table$series, -corr_table$tHo), ]
        }

return(corr_table)

}
