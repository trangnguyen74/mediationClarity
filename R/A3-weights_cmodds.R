



#### weights_cm.odds ###########################################################


#' Estimates odds weights for the Ypred estimator
#' @inheritParams weights_med
#' @param cm.vars.std Names of covariates and mediators whose mean differences are to be standardized in balance plotting. Ignore if \code{plot==FALSE}.
#' @param cm.order Order in which covariates and mediators are to be plotted. If not specify, use order that appears in \code{a.cm.form}.
#' @return \tabular{ll}{
#' \code{w.dat} \tab A data frame for the pseudo and actual subsamples (pseudo ones with estimated distribution-morphing odds weights) \cr
#' \code{wt.dist.plot} \tab If \code{plot==TRUE}, a density plot of the distribution-morphing weights, and (if sampling weights vary) a density plot of the final weights \cr
#' \code{balance.plot} \tab If \code{plot==TRUE}, a plot of the balance in covariates and mediators of the pseudo subsample(s).
#' }
#' @family weighting schemes
#' @export


weights_cm.odds <- function(
    data,
    s.wt.var = NULL,
    cross.world = "10", # options: "10", "01", "both"
    a.cm.form,

    max.stabilized.wt = 30,

    plot = TRUE,
    cm.order = NULL,
    cm.vars.std = NULL
) {


    # CLEAN INPUTS

    cm.vars <- NULL

    .prep_cm.odds()


    # COMPUTE WEIGHTS

    w.dat <- .compute_weights.cm.odds(data              = data,
                                      cross.world       = cross.world,
                                      a.cm.form         = a.cm.form,
                                      max.stabilized.wt = max.stabilized.wt)


    # MAKE PLOTS

    if (plot) plots <- .plot_cm.odds(w.dat = w.dat,
                                     cm.vars = cm.vars,
                                     cm.vars.std = cm.vars.std)




    # OUTPUT

    if (plot)
        return(mget(c("w.dat", "plots")))

    return(w.dat)

}




#### .prep_cm.odds #########################################################

#' @rdname dot-prep
#' @order 3
.prep_cm.odds <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_weights.cm.odds(top.env)

    if (top.env$plot) .check_plot.cm.odds(top.env)
}




#### .clean_weights.cm.odds ####################################################

#' @rdname dot-clean_weights
#' @order 3

.clean_weights.cm.odds <- function(env) {

    if (!is.numeric(env$max.stabilized.wt))
        stop("max.stabilized.wt must be a numeric value.")


    a.var   <- all.vars(formula(env$a.cm.form)[[2]])
    cm.vars <- all.vars(formula(env$a.cm.form)[[3]])

    stray.vars <- setdiff(c(a.var, cm.vars), names(env$data))

    if (length(stray.vars)>0)
        stop(paste("Variable(s)",
                   paste(stray.vars, collapse = ", "),
                   "in a.cm.form not found in dataset."))

    if (!is_binary01(env$data[, a.var]))
        stop(paste("Treatment variable (",
                   a.var,
                   ") must be numeric and in binary 0/1 form."))

    env$data$.a <- env$data[, a.var]


    env$cm.vars <- cm.vars
}



#### .check_plot.cm.odds ###################################################

#' @rdname dot-check_plot
#' @order 3

.check_plot.cm.odds <- function(env) {

    cm.vars     <- env$cm.vars
    cm.order    <- env$cm.order
    cm.vars.std <- env$cm.vars.std


    if (!is.null(cm.order)) {

        if (!setequal(cm.vars, cm.order)) {
            warning("Variables in cm.order do not match covariates from a.c.form. Ignoring c.order.")
        } else
            env$cm.vars <- cm.vars <- cm.order
    }




    if (is.null(cm.vars.std)) {

        maybe.cont <- sapply(cm.vars.std,
                             function(z) maybe_continuous(env$data[, z]))

        if (any(maybe.cont))
            message(paste("Consider whether the balance plot should use standardized mean differences for numeric covariate/mediators",
                          paste(cm.vars.std[which(maybe.cont)],
                                collapse = ", "),
                          "(if they are continuous variables).",
                          "To turn off this message, specify cm.vars.std=\"\"."))

        return()
    }


    if (length(cm.vars.std)==1 && cm.vars.std=="")  return()



    cm.vars.std <- setdiff(cm.vars.std, "")

    if (length(setdiff(cm.vars.std, cm.vars))>0)
        stop("Variables specified in cm.vars.std are not all contained in model formula a.cm.form.")



    ok.std <- sapply(cm.vars.std, function(z) maybe_continuous(env$data[, z]))

    if (!all(ok.std))
        stop(paste("Check variable(s)",
                   paste(cm.vars.std[which(!ok.std)], collapse = ", "),
                   "before proceeding. Only include continuous variables in c.cm.vars.std and m.cm.vars.std."))
}




#### .compute_weights.cm.odds ##################################################

#' @rdname dot-compute_weights
#' @order 3

.compute_weights.cm.odds <- function(
    data,
    cross.world,
    a.cm.form,
    max.stabilized.wt

) {

    a.cm.fu <- glm(formula = a.cm.form,
                   data    = data,
                   weights = data$.s.wt,
                   family  = quasibinomial)



    max.wt <- lapply(list(control=0, treat=1), function(z) {
        max.stabilized.wt * (sum(data$.s.wt) / sum(data$.s.wt * (data$.a==z)))
    })




    s00 <- data[data$.a==0, ];  s00$.samp <- "s00"
    s11 <- data[data$.a==1, ];  s11$.samp <- "s11"

    s00$.w.wt <- 1
    s11$.w.wt <- 1

    s00$.f.wt <- s00$.s.wt
    s11$.f.wt <- s11$.s.wt

    out <- rbind(s00, s11)




    if ("10" %in% cross.world) {

        s10 <- data[data$.a==1, ];  s10$.samp <- "s10"

        s10$.w.wt <- exp(-predict(a.cm.fu, newdata = s10, type = "link"))
        s10$.w.wt <- .trunc_right(s10$.w.wt, max.wt$treat)

        s10$.f.wt <- s10$.s.wt * s10$.w.wt

        out <- rbind(out, s10)
    }


    if ("01" %in% cross.world) {

        s01 <- data[data$.a==0, ];  s01$.samp <- "s01"

        s01$.w.wt <- exp(predict(a.cm.fu, newdata = s01, type = "link"))
        s01$.w.wt <- .trunc_right(s01$.w.wt, max.wt$control)

        s01$.f.wt <- s01$.s.wt * s01$.w.wt

        out <- rbind(out, s01)
    }

    out

}



#### .plot_cm.odds #############################################################

#' @rdname dot-plot_w.dat
#' @order 3
#' @param cm.vars Names of covariates and mediators.
#' @param cm.vars.std Names of covariates and mediators for which to use standardized mean differences in balance plots.

.plot_cm.odds <- function(w.dat,
                          cm.vars,
                          cm.vars.std) {

    out <- .plot_wt_dist.cm.odds(w.dat)
    out$balance <- .plot_balance.cm.odds(w.dat = w.dat,
                                         cm.vars = cm.vars,
                                         cm.vars.std = cm.vars.std)

    out
}




#### .plot_wt_dist.cm.odds #####################################################

#' @rdname dot-plot_wt_dist
#' @order 2

.plot_wt_dist.cm.odds <- function(
    w.dat,
    point.alpha = .1,
    jitter.width = .3
) {


    w.dat$.w.wt <- stabilize_weight(weight   = w.dat$.w.wt,
                                    group    = w.dat$.samp,
                                    s.weight = w.dat$.s.wt)


    if (is_constant(w.dat$.s.wt)) {

        dat <- w.dat[, c(".samp", ".w.wt")]
        dat <- dat[dat$.samp %in% c("s10", "s01"), ]


        p <- ggplot(data = dat,
                    aes(x = .data$.samp,
                        y = .data$.w.wt)) +
            geom_jitter(height = 0,
                        width = jitter.width,
                        alpha = point.alpha) +
            geom_violin(color = "red", fill = NA) +
            labs(x = "",
                 y = "distribution morphing weights (stabilized)") +
            theme_bw()

        return(list(w.wt.distribution = p))
    }




    w.dat$.s.wt <- stabilize_weight(weight   = w.dat$.s.wt,
                                    group    = w.dat$.samp,
                                    s.weight = rep(1, nrow(w.dat)))

    dat <- w.dat[, c(".samp", ".s.wt", ".w.wt")]
    dat <- dat[dat$.samp %in% c("s10", "s01"), ]

    w.wt.distribution <-
        ggplot(data = dat,
               aes(x = .data$.samp,
                   y = .data$.w.wt,
                   weight = .data$.s.wt)) +
        geom_jitter(height = 0,
                    width = jitter.width,
                    alpha = point.alpha,
                    aes(size = .data$.s.wt)) +
        geom_violin(color = "red", fill = NA) +
        labs(x = "",
             y = "distribution morphing weights (stabilized)",
             size = "sampling weight (stabilized)") +
        theme_bw() +
        theme(legend.position = "bottom")




    w.dat$.f.wt <- stabilize_weight(weight   = w.dat$.f.wt,
                                    group    = w.dat$.samp,
                                    s.weight = rep(1, nrow(w.dat)))

    dat <- w.dat[, c(".samp", ".f.wt")]

    f.wt.distribution <-
        ggplot(data = dat,
               aes(x = .data$.samp,
                   y = .data$.f.wt)) +
        geom_jitter(height = 0,
                    width = jitter.width,
                    alpha = point.alpha) +
        geom_violin(color = "red", fill = NA) +
        labs(x = "",
             y = "final weights (stabilized)\n(combining sampling and distribution morphing)") +
        theme_bw()


    mget(c("w.wt.distribution", "f.wt.distribution"))
}




#### .plot_balance.cm.odds #####################################################


#' Plot balance
#'
#' Internal function that makes balance plot for mediation weights and TE weights.
#' @param w.dat Data for (pseudo) subsamples, e.g. \code{output of .compute_weights.Ypred()}
#' @param cm.vars Names of covariates and mediators, already cleaned.
#' @param cm.vars.std Names of mediators whose mean differences are to be standardized, already cleaned
#' @importFrom ggplot2 ggplot aes geom_vline geom_point scale_color_manual scale_shape_manual labs theme_bw facet_wrap xlim
#' @importFrom rlang .data
#' @return Plot of balance on the means of covariates and mediators between pseudo subsample(s) and relevant subsample(s).
#' @rdname dot-plot_balance
#' @order 3

.plot_balance.cm.odds <- function(w.dat,
                                  cm.vars,
                                  cm.vars.std) {


    smd.dat <- .get_smd.cm.odds(w.dat   = w.dat,
                              vars     = cm.vars,
                              standardize = cm.vars.std)


    ggplot(data = smd.dat,
           aes(x = .data$mean.diff,
               y = factor(.data$variable,
                          levels = rev(levels(.data$variable))))) +
        geom_vline(xintercept = 0,
                   color = "gray60") +
        geom_point(aes(shape = .data$contrast.type),
                   fill = "white",
                   size = 1.5,
                   stroke = .5) +
        labs(x = "differences in means", y = "") +
        scale_shape_manual(name = "", values = c(21, 19)) +
        theme_bw() +
        xlim(min(-.3, min(smd.dat$mean.diff)),
             max( .3, max(smd.dat$mean.diff))) +
        facet_wrap(~ .data$contrast, ncol = 3)
}




#### .get_smd.cm.odds ##########################################################

#' @rdname dot-get_smd
#' @order 3
#' @param cm.vars Names of covariates and mediators, checked and dummied.
#' @param cm.vars.std Covariates and mediators to be standardized, already checked.

.get_smd.cm.odds <- function(w.dat,
                             vars,
                             standardize) {


    tmp <- .make_dummies(data = w.dat,
                         columns = vars,
                         output.names = TRUE,
                         warning = FALSE)

    vars <- tmp$columns
    w.dat  <- tmp$data;   rm(tmp)





    w.dat <- w.dat[, c(".samp", ".s.wt", ".f.wt", vars)]

    yes.s10 <- (sum(w.dat$.samp=="s10") > 0)
    yes.s01 <- (sum(w.dat$.samp=="s01") > 0)

    s00 <- w.dat[w.dat$.samp=="s00", ]
    s11 <- w.dat[w.dat$.samp=="s11", ]

    if (yes.s10) s10 <- w.dat[w.dat$.samp=="s10", ]
    if (yes.s01) s01 <- w.dat[w.dat$.samp=="s01", ]

    rm(w.dat)



    # denominator for standardization
    std.denom <- sapply(vars, function(z) {

        if (!z %in% standardize) return(1)

        .get_sd.pooled(variable = z, dat1 = s11, dat0 = s00)
    })



    # compute (standardized) mean differences
    smd <- lapply(list(unw = ".s.wt", wtd = ".f.wt"), function(w) {

        if (yes.s10) {

            means.s10 <- sapply(vars,
                                function(z) .wtd_mean(s10[, z], s10[, w]))
            means.s00 <- sapply(vars,
                                function(z) .wtd_mean(s00[, z], s00[, w]))

            diff <- cbind(s10.s00  = (means.s10 - means.s00) / std.denom)
            rm(means.s10, means.s00)
        }

        if (yes.s01) {

            means.s01 <- sapply(vars,
                                function(z) .wtd_mean(s01[, z], s01[, w]))
            means.s11 <- sapply(vars,
                                function(z) .wtd_mean(s11[, z], s11[, w]))

            diff <- cbind(diff,
                          s01.s11  = (means.s01 - means.s11) / std.denom)
            rm(means.s01, means.s11)
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


    smd$contrast <-factor(smd$contrast,
                          levels = c("s10.s00", "s01.s11"),
                          labels = c("s10 - s00", "s01 - s11"))


    smd$variable <- ifelse(smd$variable %in% standardize,
                           paste0("*", smd$variable),
                           smd$variable)

    vars <- ifelse(vars %in% standardize,
                   paste0("*", vars),
                   vars)

    smd$variable <- factor(smd$variable, levels = vars)



    smd[,c("variable", "contrast.type", "contrast", "mean.diff")]

}









