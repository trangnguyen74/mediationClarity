#### .boot_ci.se ###############################################################

#'
#'
#' @importFrom stats quantile sd

.boot_ci.se <- function(
    data,
    stratify,
    boot.num,
    seed,
    method,
    FUN,
    FUN.inputs
) {

    if (stratify) { stratify <- data$.a
    } else        { stratify <- NULL
    }

    boot.wts <- .get_boot.wts(s.weights = data$.s.wt,
                              stratify = stratify,
                              boot.num = boot.num,
                              seed = seed,
                              method = method)


    ests <- sapply(1:boot.num, function(z) {

        dat <- data
        dat$.s.wt <- boot.wts[z, ]

        do.call(FUN, c(FUN.inputs,
                       list(data        = dat,
                            output.data = FALSE)))
    })

    ci <- apply(ests, 1,
                function(z) quantile(z, prob = c(.025, .975), na.rm = TRUE))
    se <- apply(ests, 1,
                function(z) sd(z, na.rm = TRUE))

    t(rbind(ci, se = se))
}




#### .get_boot.wts ##############################################################

# TODO: add multinomial boot weights for the resampling method
# TODO: restore RNG state


#' Get weights for bootstrapping.
#'
#' @param s.weights Vector of sampling weights. If no sampling weights, use vector of 1s.
#' @param stratify Whether to draw bootstrap samples in treated and control units separately.
#' @param boot.num How many sets of bootstrap weights (or number of bootstrap samples)
#' @param seed Seed for reproducibility
#' @param method Bootstrap method. Currently allows "dirichlet", i.e., continuous weights bootstrapping where there weights are drawn from a Dirichlet distribution. TOADD: "multinmial" (or "resample")
#' @return A matrix of bootstrap weights where number of column is number of observations, and each row contains the set of bootstrap weights that forms a bootstrap sample.
#' @noRd

.get_boot.wts <- function(s.weights,
                          stratify,
                          boot.num,
                          seed,
                          method) {


    draw_boot.wts <- function(s.wts,
                              boot.num,
                              method) {

        if (method=="cont-wt") {
            draw <- gtools::rdirichlet(n = boot.num, alpha = s.wts)
        }

        draw * sum(s.wts)
    }


    set.seed(seed)

    if (is.null(stratify)) {

        out <- draw_boot.wts(s.wts = s.weights,
                             boot.num = boot.num,
                             method = method)

    } else {

        dat <- data.frame(rownum = 1:length(s.weights),
                          s.wt = s.weights,
                          stratify = factor(stratify))


        out <- lapply(levels(factor(stratify)), function(z) {

            dat.s <- dat[dat$stratify==z, ]

            wts.s <- draw_boot.wts(s.wts = dat.s$s.wt,
                                   boot.num = boot.num,
                                   method = method)

            colnames(wts.s) <- dat.s$rownum

            wts.s
        })


        out <- do.call("cbind", out)

        out <- out[, order(colnames(out))]
        colnames(out) <- NULL
    }

    out
}
