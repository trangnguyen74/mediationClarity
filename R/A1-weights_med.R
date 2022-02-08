



#### OK  weights_med ###########################################################


#' Estimate weights for natural (in)direct effects estimation
#'
#' Estimates inverse probability weights that form pseudo treated and control samples and cross-world weights that form (one or two) pseudo cross-world samples.
#' @inheritParams estimate_effects
#' @param a.c.form Formula for the P(A|C) model (the propensity score model).
#' @param a.cm.form Formula for the P(A|C,M) model.
#' @param max.stabilized.wt Max stabilized weight allowed. Larger weights are truncated to this level. Default is 30.
#' @param plot Whether to output weight distribution and balance plots. Defaults to TRUE.
#' @param c.order Order in which covariates are to be plotted. If not specify, use the order that appears in \code{a.c.form}.
#' @param m.order Order in which mediators are to be plotted. If not specify, use the order that appears in \code{a.cm.form}.
#' @param c.vars.std Covariates whose mean differences are to be standardized in balance plot. Ignore if \code{plot==FALSE}.
#' @param m.vars.std Mediators whose mean differences are to be standardized in balance plot. Ignore if \code{plot==FALSE}.
#' @return A list including\itemize{
#' \item{w.dat}{A data frame for the pseudo samples with estimated weights.}
#' \item{plot.wts}{A plot of the distributions of the weights.}
#' \item{plot.balance}{A plot of the balance in covariates and mediators of the pseudo samples.}
#' }
#' @family weighting schemes
#' @export

weights_med <- function(
    data,
    s.wt.var = NULL,
    cross.world = "10", # options: "10", "01", "both"
    a.c.form,
    a.cm.form,

    max.stabilized.wt = 30,

    plot = TRUE,
    c.order = NULL,
    m.order = NULL,
    c.vars.std = NULL,
    m.vars.std = NULL
) {

    # CLEAN INPUTS

    c.vars <- m.vars <- NULL

    .prep_med()


    # COMPUTE WEIGHTS

    w.dat <- .compute_weights.med(data              = data,
                                   cross.world       = cross.world,
                                   a.c.form          = a.c.form,
                                   a.cm.form         = a.cm.form,
                                   max.stabilized.wt = max.stabilized.wt)


    # MAKE PLOTS

    if (plot)  plots <- .plot_med(w.dat = w.dat,
                                  c.vars = c.vars,
                                  m.vars = m.vars,
                                  c.vars.std = c.vars.std,
                                  m.vars.std = m.vars.std)



    # OUTPUT

    if (plot)  return(list(w.dat = w.dat, plots  = plots))

    return(w.dat)

}



#### OK  .prep_med #############################################################

#' Internal functions: preparing inputs before main computing
#'
#' A set of internal functions to do the prep work for the main functions before the key computing takes place.
#' @name dot-prep
#' @keywords internal
NULL

#' @rdname dot-prep
#' @order 1

.prep_med <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_weights.med(top.env)

    if (top.env$plot) .check_plot.med(top.env)
}







#### OK  .clean_weights.med ####################################################

#' Internal functions: clean inputs for weighting
#'
#' Functions called by \code{.prep_} functions within \code{weights_} and \code{estimate_} functions. If no weighting is required, grab treatment variable.
#' @inheritParams env-block
#' @name dot-clean_weights
#' @keywords internal
NULL



#' @rdname dot-clean_weights
#' @order 1
#' @details \code{.clean_weights.med()} is used by \code{.prep_med()}.

.clean_weights.med <- function(env) {

    if (!is.numeric(env$max.stabilized.wt))
        stop("max.stabilized.wt must be a numeric value.")


    a.c  <- env$a.c.form
    a.cm <- env$a.cm.form


    if (!formula(a.c)[[2]]==formula(a.cm)[[2]])
        stop("Treatment variable not the same in a.c.form and a.cm.form.")

    if (!all(all.vars(formula(a.c)[[3]]) %in%
             all.vars(formula(a.cm)[[3]])))
        stop("Some C variable(s) in a.c.form but not in a.cm.form")

    a.var  <- all.vars(formula(a.c)[[2]])
    c.vars <- all.vars(formula(a.c)[[3]])
    m.vars <- setdiff(all.vars(formula(a.cm)[[3]]),
                      all.vars(formula(a.c)[[3]]))

    stray.vars <- setdiff(c(a.var, c.vars, m.vars), names(env$data))

    if (length(stray.vars)>0)
        stop(paste("Variable(s)", paste(stray.vars, collapse = ", "), "not found in dataset."))

    if (!is_binary01(env$data[, a.var]))
        stop(paste("Treatment variable (", a.var, ") must be numeric and in binary 0/1 form."))

    env$data$.a <- env$data[, a.var]

    env$c.vars <- c.vars
    env$m.vars <- m.vars

}




#### OK  .check_plot.med #######################################################

#' Internal functions: check inputs for plotting
#'
#' Functions called by \code{.prep_} functions within \code{weights_} an \code{estimate_} functions.
#' @inheritParams env-block
#' @name dot-check_plot
#' @keywords internal
NULL

#' @rdname dot-check_plot
#' @order 1
#' @details \code{.check_plot.med()} is called by \code{.prep_med()}.

.check_plot.med <- function(env) {

    c.vars <- env$c.vars
    m.vars <- env$m.vars
    c.order <- env$c.order
    m.order <- env$m.order
    c.vars.std <- env$c.vars.std
    m.vars.std <- env$m.vars.std


    if (!is.null(c.order)) {

        if (!setequal(c.vars, c.order)) {
            warning("Variables in c.order do not match covariates from a.c.form. Ignoring c.order.")
        } else
            env$c.vars <- c.vars <- c.order
    }


    if (!is.null(m.order)) {

        if (!setequal(m.vars, m.order)) {
            warning("Variables in m.order do not match mediators obtained from the combination of a.c.form and a.cm.form. Ignoring m.order.")
        } else
            env$m.vars <- m.vars <- m.order
    }



    if (is.null(c.vars.std) && is.null(m.vars.std)) {

        maybe.cont <- sapply(c(c.vars, m.vars),
                             function(z) maybe_continuous(env$data[, z]))

        if (any(maybe.cont))
            message(paste("Consider whether the balance plot should use standardized mean differences for numeric covariate/mediators",
                          paste(c(c.vars, m.vars)[which(maybe.cont)],
                                collapse = ", "),
                          "(if they are continuous variables).",
                          "To turn off this message, specify c.vars.std=\"\", m.vars.std=\"\"."))

        return()
    }


    if (length(c.vars.std==1) && c.vars.std=="" &&
        length(m.vars.std==1) && m.vars.std=="")
        return()




    c.vars.std <- setdiff(c.vars.std, "")
    m.vars.std <- setdiff(m.vars.std, "")

    if (length(setdiff(c.vars.std, c.vars))>0)
        stop("Variables specified in c.vars.std are not all contained in model formula a.c.form.")

    if (length(setdiff(m.vars.std, m.vars))>0)
        stop("Variables specified in m.vars.std are not all part of M variables based on formulas a.c.form and a.cm.form.")



    vars.std <- c(c.vars.std, m.vars.std)

    ok.std <- sapply(vars.std, function(z) maybe_continuous(env$data[, z]))

    if (!all(ok.std))
        stop(paste("Check variable(s)",
                   paste(vars.std[which(!ok.std)], collapse = ", "),
                   "before proceeding. Only include continuous variables in c.vars.std and m.vars.std."))

}





#### OK  .compute_weights.med ##################################################

#' (For maintainer) Compute weights
#'
#' Internal functions that compute weights using cleaned inputs.
#' @param data Prepared internal dataset (e.g., result of \code{.pred_data()}).
#' @param cross.world Cleaned cross.world (e.g., from \code{.clean_cross.world()}).
#' @param a.c.form,a.cm.form Cleaned formulas (e.g., from \code{.clean_a.forms()}).
#' @param max.stabilized.wt A scalar.
#' @return A data frame for the pseudo samples with estimated weights, in stacked format.
#' @name dot-compute_weights
#' @keywords internal
NULL

#' @rdname dot-compute_weights
#' @order 1

.compute_weights.med <- function(
    data,
    cross.world,
    a.c.form,
    a.cm.form,
    max.stabilized.wt

) {

    a.c.fu <- glm(formula = a.c.form,
                  data    = data,
                  weights = data$.s.wt,
                  family  = quasibinomial)
    a.cm.fu <- glm(formula = a.cm.form,
                   data    = data,
                   weights = data$.s.wt,
                   family  = quasibinomial)



    max.wt <- lapply(list(control=0, treat=1), function(z) {
        max.stabilized.wt * (sum(data$.s.wt) / sum(data$.s.wt * (data$.a==z)))
    })




    p00 <- data[data$.a==0, ];  p00$.samp <- "p00"
    p11 <- data[data$.a==1, ];  p11$.samp <- "p11"

    p00$.w.wt <- 1 + exp( predict(a.c.fu, newdata = p00, type = "link"))
    p11$.w.wt <- 1 + exp(-predict(a.c.fu, newdata = p11, type = "link"))

    p00$.w.wt <- .trunc_right(p00$.w.wt, max.wt$control)
    p11$.w.wt <- .trunc_right(p11$.w.wt, max.wt$treat)

    p00$.f.wt <- p00$.s.wt * p00$.w.wt
    p11$.f.wt <- p11$.s.wt * p11$.w.wt

    out <- rbind(p00, p11)



    if ("10" %in% cross.world) {

        p10 <- data[data$.a==1, ];  p10$.samp <- "p10"

        p10$.w.wt <-
            exp(-predict(a.cm.fu, newdata = p10, type = "link")) *
            (1 + exp(predict(a.c.fu, newdata = p10, type = "link")))

        p10$.w.wt <- .trunc_right(p10$.w.wt, max.wt$treat)

        p10$.f.wt <- p10$.s.wt * p10$.w.wt

        out <- rbind(out, p10)
    }


    if ("01" %in% cross.world) {

        p01 <- data[data$.a==0, ];  p01$.samp <- "p01"

        p01$.w.wt <-
            exp(predict(a.cm.fu, newdata = p01, type = "link")) *
            (1 + exp(-predict(a.c.fu, newdata = p01, type = "link")))

        p01$.w.wt <- .trunc_right(p01$.w.wt, max.wt$control)

        p01$.f.wt <- p01$.s.wt * p01$.w.wt

        out <- rbind(out, p01)
    }

    out

}




#### OK  .plot_med #############################################################

#' (For maintainer) Plot weight distributions and mean balance
#'
#' @name dot-plot_w.dat
#' @keywords internal
NULL


#' @rdname dot-plot_w.dat
#' @order 1
#' @param w.dat Weighted data
#' @param c.vars Names of covariates, already cleaned.
#' @param m.vars Names of mediators, already cleaned.
#' @param c.vars.std Names of covariates whose mean differences are to be standardized, already cleaned.
#' @param m.vars.std Names of mediators whose mean differences are to be standardized, already cleaned.
#' @param key.balance Whether all balance components are crucial to the estimator using the weighting method, defaults to FALSE. See relevance in \bold{Value} section.

.plot_med <- function(w.dat,
                      c.vars,
                      m.vars,
                      c.vars.std,
                      m.vars.std,
                      key.balance = FALSE) {


    out <- .plot_wt_dist(w.dat)


    if (key.balance) { bal.name <- "key.balance"
    } else           { bal.name <- "balance"
    }

    out[[bal.name]] <- .plot_balance.med(w.dat = w.dat,
                                         c.vars = c.vars,
                                         m.vars = m.vars,
                                         c.vars.std = c.vars.std,
                                         m.vars.std = m.vars.std,
                                         key.balance = key.balance)

    out
}





#### OK  .plot_balance.med #####################################################

#' (For maintainer) Plot balance
#'
#' Internal functions that makes balance plots for weighted data.
#' @inheritParams .plot_med
#' @importFrom ggplot2 ggplot aes geom_vline geom_point scale_color_manual scale_shape_manual labs theme_bw facet_wrap xlim
#' @importFrom rlang .data
#' @return Plot of balance on the means of covariates and mediators between relevant pseudo samples and full sample. If \code{estimate_wtd==TRUE}, add "(anchor)" and "(for <effect>)" notes to plot labels to draw attention to how each balance matters to the estimator.
#' @name dot-plot_balance
#' @keywords internal
NULL

#' @rdname dot-plot_balance
#' @order 1

.plot_balance.med <- function(w.dat,
                              c.vars,
                              m.vars,
                              c.vars.std,
                              m.vars.std,
                              key.balance = FALSE) {

    smd.dat <- .get_smd.med(w.dat = w.dat,
                            c.vars = c.vars,
                            m.vars = m.vars,
                            c.vars.std = c.vars.std,
                            m.vars.std = m.vars.std)

    if (key.balance)
        smd.dat$contrast <-
        factor(smd.dat$contrast,
               levels = c("p11 - full", "p00 - full", "p11 - p00",
                          "p10 - full", "p11 - p10", "p10 - p00",
                          "p01 - full", "p11 - p01", "p01 - p00"),
               labels = c("p11 - full  (anchor)",
                          "p00 - full  (anchor)",
                          "p11 - p00  (for TE)",
                          "p10 - full  (anchor)",
                          "p11 - p10  (for NIE1/NDE0)",
                          "p10 - p00  (for NIE1/NDE0)",
                          "p01 - full  (anchor)",
                          "p11 - p01  (for NIE0/NDE1)",
                          "p01 - p00  (for NIE0/NDE1)"))

    ggplot(data = smd.dat,
           aes(x = .data$mean.diff,
               y = factor(.data$variable,
                          levels = rev(levels(.data$variable))))) +
        geom_vline(xintercept = 0,
                   color = "gray60") +
        geom_point(aes(color = .data$var.type,
                       shape = .data$contrast.type),
                   fill = "white",
                   size = 1.5,
                   stroke = .5) +
        labs(x = "differences in means",
             y = "") +
        scale_color_manual(name = "", values = c("black", "magenta")) +
        scale_shape_manual(name = "", values = c(21, 19)) +
        theme_bw() +
        xlim(min(c(-.3, smd.dat$mean.diff)),
             max(c( .3, smd.dat$mean.diff))) +
        facet_wrap(~.data$contrast, ncol = 3)

}



#### OK  .get_smd.med ##########################################################

#' (For maintainer) Compute  (standardized) mean differences
#'
#' Internal functions called by \code{.plot_balance.} functions to compute (standardized) mean differences to be plotted.
#' @return A data frame containing the (standardized) mean differences.
#' @keywords internal
#' @name dot-get_smd
#' @keywords internal
NULL

#' @rdname dot-get_smd
#' @order 1
#' @param w.dat Dataset for pseudo samples, already cleaned and dummy coded.
#' @param c.vars Names of covariates, checked and dummied.
#' @param m.vars Names of mediators, checked and dummied.
#' @param c.vars.std Covariates to be standardized, already checked.
#' @param m.vars.std Mediators to be standardized, already checked.

.get_smd.med <- function(w.dat,
                         c.vars,
                         m.vars,
                         c.vars.std,
                         m.vars.std) {

    tmp <- .dummies_2sets(data = w.dat,
                          columns1 = c.vars,
                          columns2 = m.vars)

    w.dat <- tmp$data
    c.vars <- tmp$columns1
    m.vars <- tmp$columns2;  rm(tmp)

    c.smd <- .get_c.smd(w.dat      = w.dat,
                        vars        = c.vars,
                        standardize = c.vars.std)

    m.smd <- .get_m.smd(w.dat      = w.dat,
                        vars        = m.vars,
                        standardize = m.vars.std)

    smd <- rbind(c.smd, m.smd)

    c.vars <- ifelse(c.vars %in% c.vars.std, paste0("*", c.vars), c.vars)
    m.vars <- ifelse(m.vars %in% m.vars.std, paste0("*", m.vars), m.vars)

    smd$variable <- factor(smd$variable, levels = c(c.vars, m.vars))

    smd
}




#### .get_c.smd & .get_m.smd ###################################################

#' Covariate or mediator (standardized) mean differences
#'
#' Internal functions called by \code{.plot_balance()} or \code{.plot_balance.Ypred()} to compute covariate/mediator (standardized) mean differences required for balance plotting.
#' @param w.dat Data for pseudo samples, e.g., output of \code{.compute_weights.med()} or \code{.compute_weights.te()} that has had categorical variables dummy coded.
#' @param vars Names of numeric variables on which to compute mean differences.
#' @param w.dat For \code{.get_smd.Ypred()}: data for (pseudo) subsamples, e.g., output of \code{.compute_weights.Ypred()} that has had categorical variables dummy coded.
#' @param m.vars For \code{.get_smd.Ypred()}: names of mediators.
#' @param standardize Names of variables for which the mean differences are to be standardized.
#' @return A data frame holding the (standardized) mean differences.
#' @return With \code{.get_c.smd()}, these are covariate mean differences between pseudo samples and full sample and across pseudo samples.
#' @return With \code{.get_m.smd()}, these are mediator mean differences between relevant pseudo samples (i.e., between \code{p10} and \code{p00} and/or between \code{p01} and \code{p11}).
#' @return With \code{.get_smd.Ypred()}, these are covariate and mediator mean differences between pseudo cross-world SUBsamples and the original subsamples they mimic (i.e., between \code{s10} and \code{s00} and/or between \code{s01} and \code{s11}).
#' @noRd


.get_c.smd <- function(w.dat,
                       vars,
                       standardize) {

    w.dat <- w.dat[, c(".samp", ".s.wt", ".f.wt", vars)]

    yes.p10 <- (sum(w.dat$.samp=="p10") > 0)
    yes.p01 <- (sum(w.dat$.samp=="p01") > 0)



    # separate pseudo samples and full sample
    p11 <- w.dat[w.dat$.samp=="p11", ]
    p00 <- w.dat[w.dat$.samp=="p00", ]

    if (yes.p10) p10 <- w.dat[w.dat$.samp=="p10", ]
    if (yes.p01) p01 <- w.dat[w.dat$.samp=="p01", ]

    rm(w.dat)

    full <- rbind(p11, p00)
    full$.f.wt <- full$.s.wt



    # denominator for standardization
    std.denom <- sapply(vars, function(z) {

        if (!z %in% standardize) return(1)

        .get_sd.pooled(variable = z, dat1 = p11, dat0 = p00)
    })



    # compute (standardized) mean differences
    smd <- lapply(list(unw = ".s.wt", wtd = ".f.wt"), function(w) {

        means.fu  <- sapply(vars, function(z) .wtd_mean(full[, z], full[, w]))
        means.p11 <- sapply(vars, function(z) .wtd_mean(p11[, z],  p11[, w]))
        means.p00 <- sapply(vars, function(z) .wtd_mean(p00[, z],  p00[, w]))

        diff <- cbind(p11.full = (means.p11 - means.fu) / std.denom,
                      p00.full = (means.p00 - means.fu) / std.denom,
                      p11.p00  = (means.p11 - means.p00) / std.denom)


        if (yes.p10) {

            means.p10 <- sapply(vars,
                                function(z) .wtd_mean(p10[, z], p10[, w]))

            diff <- cbind(diff,
                          p10.full = (means.p10 - means.fu) / std.denom,
                          p11.p10  = (means.p11 - means.p10) / std.denom,
                          p10.p00  = (means.p10 - means.p00) / std.denom)
            rm(means.p10)
        }

        if (yes.p01) {

            means.p01 <- sapply(vars,
                                function(z) .wtd_mean(p01[, z], p01[, w]))

            diff <- cbind(diff,
                          p01.full = (means.p01 - means.fu) / std.denom,
                          p11.p01  = (means.p11 - means.p01) / std.denom,
                          p01.p00  = (means.p01 - means.p00) / std.denom)
            rm(means.p01)
        }

        diff.names <- colnames(diff)

        diff <- data.frame(diff, row.names = NULL)
        diff$variable <- vars

        diff <- reshape_gather(data     = diff,
                               columns  = diff.names,
                               key      = "contrast",
                               value    = "mean.diff",
                               wide.row = FALSE)

        if (w==".f.wt") diff$contrast.type <- "weighted"
        if (w==".s.wt") diff$contrast.type <- "pre-weighting"

        diff
    })

    smd$unw <- smd$unw[!smd$unw$contrast %in% c("p11.p10", "p01.p00"), ]

    smd <- do.call("rbind", smd)

    rownames(smd) <- NULL


    smd$contrast <-
        factor(smd$contrast,
               levels = c("p11.full", "p00.full", "p11.p00",
                          "p10.full", "p11.p10", "p10.p00",
                          "p01.full", "p11.p01", "p01.p00"),
               labels = c("p11 - full", "p00 - full", "p11 - p00",
                          "p10 - full", "p11 - p10",  "p10 - p00",
                          "p01 - full", "p11 - p01",  "p01 - p00"))

    smd$var.type <- "covariate"

    smd <- smd[,c("var.type", "variable", "contrast.type", "contrast", "mean.diff")]

    smd$variable <- ifelse(smd$variable %in% standardize,
                           paste0("*", smd$variable),
                           smd$variable)

    smd

}


#' @noRd
.get_m.smd <- function(w.dat,
                       vars,
                       standardize) {

    w.dat <- w.dat[, c(".samp", ".s.wt", ".f.wt", vars)]

    yes.p10 <- (sum(w.dat$.samp=="p10") > 0)
    yes.p01 <- (sum(w.dat$.samp=="p01") > 0)


    # separate pseudo samples
    p00 <- w.dat[w.dat$.samp=="p00", ]
    p11 <- w.dat[w.dat$.samp=="p11", ]

    if (yes.p10) p10 <- w.dat[w.dat$.samp=="p10", ]
    if (yes.p01) p01 <- w.dat[w.dat$.samp=="p01", ]

    rm(w.dat)


    # denominator for standardization
    std.denom <- sapply(vars, function(z) {

        if (!z %in% standardize) return(1)

        .get_sd.pooled(variable = z, dat1 = p11, dat0 = p00)
    })



    # compute mean differences
    smd <- lapply(list(unw = ".s.wt", wtd = ".f.wt"), function(w) {

        diff <- NULL

        if (yes.p10) {

            means.p10 <- sapply(vars,
                                function(z) .wtd_mean(p10[, z], p10[, w]))
            means.p00 <- sapply(vars,
                                function(z) .wtd_mean(p00[, z],  p00[, w]))

            diff <- cbind(diff,
                          p10.p00 = (means.p10 - means.p00) / std.denom)
        }

        if (yes.p01) {

            means.p11 <- sapply(vars,
                                function(z) .wtd_mean(p11[, z],  p11[, w]))
            means.p01 <- sapply(vars,
                                function(z) .wtd_mean(p01[, z], p01[, w]))

            diff <- cbind(diff,
                          p11.p01 = (means.p11 - means.p01) / std.denom)
        }

        diff.names <- colnames(diff)

        diff <- data.frame(diff, row.names = NULL)
        diff$variable <- vars

        diff <- reshape_gather(data     = diff,
                               columns  = diff.names,
                               key      = "contrast",
                               value    = "mean.diff",
                               wide.row = FALSE)

        if (w==".f.wt") diff$contrast.type <- "weighted"
        if (w==".s.wt") diff$contrast.type <- "pre-weighting"

        diff
    })

    smd <- do.call("rbind", smd)

    rownames(smd) <- NULL


    smd$contrast <-
        factor(smd$contrast,
               levels = c("p10.p00", "p11.p01"),
               labels = c("p10 - p00", "p11 - p01"))

    smd$var.type <- "mediator"

    smd <- smd[,c("var.type", "variable", "contrast.type", "contrast", "mean.diff")]

    smd$variable <- ifelse(smd$variable %in% standardize,
                           paste0("*", smd$variable),
                           smd$variable)

    smd

}




#### OK  param env #############################################################

#' @param env An environment
#' @name env-block
#' @keywords internal
NULL
