#' @name get_SMDs
#' @discription Compute (standardized) mean differences
#' @param
#' @return
#' @noRd

get_SMDs <- function(ps.dat,
                     sample = "sample",
                     s.wt = "s.wt",
                     f.wt = "f.wt",
                     c.vars,
                     m.vars = NULL,
                     standardize,
                     full.label = FALSE) {


    ps.dat$sample <- ps.dat[, sample]
    ps.dat$s.wt <- ps.dat[, s.wt]
    ps.dat$f.wt <- ps.dat[, f.wt]

    yes.p10 <- (sum(ps.dat$sample=="p10") > 0)
    yes.p01 <- (sum(ps.dat$sample=="p01") > 0)


    if (!is.null(m.vars) && !yes.p10 && !yes.p01)
        stop("With m.vars specified, data must include one (or both) of pseudo samples p10 and p01.")

    ps.dat <- ps.dat[, names(ps.dat) %in% c("sample", "s.wt", "f.wt",
                                            c.vars, m.vars)]

    p11 <- ps.dat[ps.dat$sample=="p11", ]
    p00 <- ps.dat[ps.dat$sample=="p00", ]


    if (yes.p10) p10 <- ps.dat[ps.dat$sample=="p10", ]
    if (yes.p01) p01 <- ps.dat[ps.dat$sample=="p01", ]

    rm(ps.dat)

    full <- rbind(p11, p00)
    full$f.wt <- full$s.wt




    size1 <- sum(p11$s.wt)
    size0 <- sum(p00$s.wt)

    sd.pooled <- sapply(c(c.vars, m.vars), function(z) {

        if (z %in% standardize) {

            p11.s2 <- weighted_var(p11[, z], p11$s.wt)
            p00.s2 <- weighted_var(p00[, z], p00$s.wt)

            pooled.s2 <- (p11.s2 * size1 + p00.s2 * size0) / (size1 + size0)

            sd <- sqrt(pooled.s2)

        } else {

            sd <- 1
        }

        sd

    })



    smds <- lapply(list(unw = "s.wt", wtd = "f.wt"), function(w) {

        c.fu  <- sapply(c.vars, function(z) weighted.mean(full[, z], full[, w]))
        c.p11 <- sapply(c.vars, function(z) weighted.mean(p11[, z],  p11[, w]))
        c.p00 <- sapply(c.vars, function(z) weighted.mean(p00[, z],  p00[, w]))

        diff <- cbind(p11.full = c.p11 - c.fu,
                      p00.full = c.p00 - c.fu,
                      p11.p00  = c.p11 - c.p00)

        rownames(diff) <- c.vars



        if (!is.null(m.vars)) {
            m.diff <- matrix(NA, ncol = 3, nrow = length(m.vars))
            rownames(m.diff) <- m.vars
        }



        if (yes.p10) {

            c.p10 <- sapply(c.vars, function(z) weighted.mean(p10[, z], p10[, w]))

            diff <- cbind(diff,
                          p10.full = c.p10 - c.fu,
                          p11.p10  = c.p11  - c.p10,
                          p10.p00  = c.p10 - c.p00);  rm(c.p10)


            if (!is.null(m.vars)) {
                m.p10 <- sapply(m.vars, function(z) weighted.mean(p10[, z], p10[, w]))
                m.p00 <- sapply(m.vars, function(z) weighted.mean(p00[, z], p00[, w]))

                m.diff <- cbind(m.diff,
                                p10.full = NA,
                                p11.p10  = NA,
                                p10.p00  = m.p10 - m.p00);  rm(m.p10, m.p00)
            }
        }



        if (yes.p01) {

            c.p01 <- sapply(c.vars, function(z) weighted.mean(p01[, z], p01[, w]))

            diff <- cbind(diff,
                          p01.full = c.p01 - c.fu,
                          p11.p01  = c.p11 - c.p01,
                          p01.p00  = c.p01 - c.p00);  rm(c.p01)

            if (!is.null(m.vars)) {
                m.p01 <- sapply(m.vars, function(z) weighted.mean(p01[, z], p01[, w]))
                m.p11 <- sapply(m.vars, function(z) weighted.mean(p11[, z], p11[, w]))

                m.diff <- cbind(m.diff,
                                p01.full = NA,
                                p11.p01  = m.p11 - m.p01,
                                p01.p00  = NA);           rm(m.p01, m.p11)
            }
        }



        if (exists("m.diff")) { diff <- rbind(diff, m.diff);  rm(m.diff)
        }



        sd.pooled <- matrix(sd.pooled, ncol = 1)
        sd.pooled <- sd.pooled[, rep(1, ncol(diff))]
        diff <- diff / sd.pooled

        diff <- data.frame(diff)
        diff.names <- names(diff)

        diff$variable <- c(c.vars, m.vars)

        diff <- reshape_gather(data = diff,
                               columns = diff.names,
                               key = "contrast",
                               value = "smd")

        diff <- diff[, names(diff)!="wide.row"]

        diff <- diff[!is.na(diff$smd), ]

        diff$var.type <- ifelse(diff$variable %in% c.vars, "covariate", "mediator")

        if (w=="f.wt") diff$type <- "weighted\n(pseudo samples)"
        if (w=="s.wt") diff$type <- "unweighted"

        diff
    })

    tmp <- smds$unw
    tmp <- tmp[!(tmp$var.type=="covariate" &
                     tmp$contrast %in% c("p11.p10", "p01.p00")), ]
    smds$unw <- tmp



    smds <- do.call("rbind", smds)


    if (yes.p10 && yes.p01) {
        smds$contrast <-
            factor(smds$contrast,
                   levels = c("p11.full", "p00.full", "p11.p00",
                              "p10.full", "p11.p10", "p10.p00",
                              "p01.full", "p11.p01", "p01.p00"))

        if (full.label)
            smds$contrast <-
                factor(smds$contrast,
                       labels = c("p11 - full  (anchor)", "p00 - full  (anchor)",
                                  "p11 - p00  (for TE)",
                                  "p10 - full  (anchor)",
                                  "p11 - p10  (for NIE1)", "p10 - p00  (for NDE0)",
                                  "p01 - full  (anchor)",
                                  "p11 - p01  (for NDE1)", "p01 - p00  (for NIE0)"))

    } else if (yes.p10) {
        smds$contrast <-
            factor(smds$contrast,
                   levels = c("p11.full", "p00.full", "p11.p00",
                              "p10.full", "p11.p10", "p10.p00"))

        if (full.label)
            smds$contrast <-
                factor(smds$contrast,
                       labels = c("p11 - full  (anchor)", "p00 - full  (anchor)",
                                  "p11 - p00  (for TE)",
                                  "p10 - full  (anchor)",
                                  "p11 - p10  (for NIE1)", "p10 - p00  (for NDE0)"))

    } else if (yes.p01) {
        smds$contrast <-
            factor(smds$contrast,
                   levels = c("p11.full", "p00.full", "p11.p00",
                              "p01.full", "p11.p01", "p01.p00"))

        if (full.label)
            smds$contrast <-
                factor(smds$contrast,
                       labels = c("p11 - full  (anchor)", "p00 - full  (anchor)",
                                  "p11 - p00  (for TE)",
                                  "p01 - full  (anchor)",
                                  "p11 - p01  (for NDE1)", "p01 - p00  (for NIE0)"))

    } else {
        smds$contrast <-
            factor(smds$contrast,
                   levels = c("p11.full", "p00.full", "p11.p00"))

        if (full.label)
            smds$contrast <-
                factor(smds$contrast,
                       labels = c("p11 - full  (anchor)", "p00 - full  (anchor)",
                                  "p11 - p00  (for TE)"))

    }

    smds

}
