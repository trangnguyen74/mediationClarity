



#### OK  weights_ipw ###########################################################



#' Estimate inverse probability weights
#'
#' This function estimates weights that form two pseudo treated and pseudo control samples (as if for total effect estimation)
#' @inheritParams weights_med
#' @importFrom stats formula glm quasibinomial predict
#' @return A list including\tabular{ll}{
#' \code{w.dat} \tab A data frame for the pseudo samples with estimated weights. \cr
#' \code{plot.wts} \tab A plot of the distributions of the weights. \cr
#' \code{plot.balance} \tab A plot of the balance in covariates and mediators of the pseudo samples. \cr
#' }
#' @family weighting schemes
#' @export

weights_ipw <- function(
    data,
    s.wt.var = NULL,
    a.form,

    max.stabilized.wt = 30,

    plot = TRUE,
    vars.std = NULL,
    vars.order = NULL
) {


    # CLEAN INPUTS

    vars <- NULL

    .prep_ipw()


    # COMPUTE WEIGHTS

    w.dat <- .compute_weights.ipw(data              = data,
                                  a.form            = a.form,
                                  max.stabilized.wt = max.stabilized.wt)


    # MAKE PLOTS

    if (plot) {

        plots <- .plot_ipw(w.dat    = w.dat,
                           vars     = vars,
                           vars.std = vars.std)
    }


    # OUTPUT

    if (!plot) return(w.dat)

    return(list(w.dat = w.dat,
                plots  = plots))

}




#### OK  .prep_ipw #############################################################

#' @rdname dot-prep
#' @order 2
.prep_ipw <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_weights.ipw(top.env)

    if (top.env$plot) .check_plot.ipw(top.env)
}




#### OK  .clean_weights.ipw ####################################################


#' @rdname dot-clean_weights
#' @order 2

.clean_weights.ipw <- function(env) {

    if (!is.numeric(env$max.stabilized.wt))
        stop("max.stabilized.wt must be a numeric value.")

    a.form <- env$a.form

    a.var  <- all.vars(formula(a.form)[[2]])
    c.vars <- all.vars(formula(a.form)[[3]])

    stray.vars <- setdiff(c(a.var, c.vars), names(env$data))

    if (length(stray.vars)>0)
        stop(paste("Variable(s)", paste(stray.vars, collapse = ", "), "not found in dataset."))

    if (!is_binary01(env$data[, a.var]))
        stop(paste("Treatment variable (", a.var, ") must be numeric and in binary 0/1 form."))

    env$data$.a <- env$data[, a.var]

    env$vars <- c.vars

}




#### OK  .check_plot.ipw #######################################################

#' @rdname dot-check_plot
#' @order 2

.check_plot.ipw <- function(env) {

    vars       <- env$vars
    vars.order <- env$vars.order
    vars.std   <- env$vars.std


    if (!is.null(vars.order)) {

        if (!setequal(vars, vars.order)) {
            warning("Variables in c.order do not match covariates from a.c.form. Ignoring c.order.")
        } else
            env$vars <- vars <- vars.order
    }





    if (is.null(vars.std)) {

        maybe.cont <- sapply(vars,
                             function(z) maybe_continuous(env$data[, z]))

        if (any(maybe.cont))
            message(paste("Consider whether the balance plot should use standardized mean differences for numeric covariates",
                          paste(c.vars[which(maybe.cont)],
                                collapse = ", "),
                          "(if they are continuous variables).",
                          "To turn off this message, specify c.vars.std=\"\"."))

        return()
    }


    if (length(vars.std)==1 && vars.std=="")  return()


    vars.std <- setdiff(vars.std, "")

    if (length(setdiff(vars.std, vars))>0)
        stop("Variables specified in c.vars.std are not all contained in model formula a.c.form.")


    ok.std <- sapply(vars.std, function(z) maybe_continuous(env$data[, z]))

    if (!all(ok.std))
        stop(paste("Check variable(s)",
                   paste(vars.std[which(!ok.std)], collapse = ", "),
                   "before proceeding. Only include continuous variables in c.vars.std."))

}




#### OK  .compute_weights.ipw ##################################################

#' @rdname dot-compute_weights
#' @order 2

.compute_weights.ipw <- function(
    data,
    a.form,
    max.stabilized.wt

) {

    a.fu <- glm(formula = a.form,
                  data    = data,
                  weights = data$.s.wt,
                  family  = quasibinomial)



    max.wt <- lapply(list(control=0, treat=1), function(z) {
        max.stabilized.wt * (sum(data$.s.wt) / sum(data$.s.wt * (data$.a==z)))
    })



    p00 <- data[data$.a==0, ];  p00$.samp <- "p00"
    p11 <- data[data$.a==1, ];  p11$.samp <- "p11"

    p00$.w.wt <- 1 + exp( predict(a.fu, newdata = p00, type = "link"))
    p11$.w.wt <- 1 + exp(-predict(a.fu, newdata = p11, type = "link"))

    p00$.w.wt <- .trunc_right(p00$.w.wt, max.wt$control)
    p11$.w.wt <- .trunc_right(p11$.w.wt, max.wt$treat)

    p00$.f.wt <- p00$.s.wt * p00$.w.wt
    p11$.f.wt <- p11$.s.wt * p11$.w.wt

    rbind(p00, p11)

}




#### OK  .plot_ipw #############################################################

#' @rdname dot-plot_w.dat
#' @order 2

.plot_ipw <- function(w.dat,
                      vars,
                      vars.std) {

    c(.plot_wt_dist(w.dat),
      .plot_balance.ipw(w.dat = w.dat,
                        vars = vars,
                        vars.std = vars.std))
}




#### OK  .plot_balance.ipw #####################################################

#' @rdname dot-plot_balance
#' @order 2

.plot_balance.ipw <- function(w.dat,
                              vars,
                              vars.std) {


    smd.dat <- .get_smd.ipw(w.dat = w.dat,
                            vars = vars,
                            standardize = vars.std)


    p <-
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
        labs(x = "differences in means",
             y = "") +
        scale_shape_manual(name = "", values = c(21, 19)) +
        theme_bw() +
        xlim(min(c(-.3, smd.dat$mean.diff)),
             max(c( .3, smd.dat$mean.diff))) +
        facet_wrap(~.data$contrast, ncol = 3)


    list(balance = p)

}




#### OK  .get_smd.ipw ##########################################################

#' @rdname dot-get_smd
#' @order 2

.get_smd.ipw <- function(w.dat,
                         vars,
                         standardize) {


    tmp <- .make_dummies(data = w.dat,
                         columns = vars,
                         output.names = TRUE,
                         warning = FALSE)

    w.dat <- tmp$data
    vars  <- tmp$columns; rm(tmp)


    w.dat <- w.dat[, c(".samp", ".s.wt", ".f.wt", vars)]

    p11 <- w.dat[w.dat$.samp=="p11", ]
    p00 <- w.dat[w.dat$.samp=="p00", ];  rm(w.dat)

    full <- rbind(p11, p00)
    full$.w.wt <- 1
    full$.f.wt <- full$.s.wt


    std.denom <- sapply(vars, function(z) {

        if (!z %in% standardize) return(1)

        .get_sd.pooled(variable = z, dat1 = p11, dat0 = p00)
    })



    smd <- lapply(list(unw = ".s.wt", wtd = ".f.wt"), function(w) {

        means.fu  <- sapply(vars, function(z) .wtd_mean(full[, z], full[, w]))
        means.p11 <- sapply(vars, function(z) .wtd_mean(p11[, z],  p11[, w]))
        means.p00 <- sapply(vars, function(z) .wtd_mean(p00[, z],  p00[, w]))

        diff <- cbind(p11.full = (means.p11 - means.fu) / std.denom,
                      p00.full = (means.p00 - means.fu) / std.denom,
                      p11.p00  = (means.p11 - means.p00) / std.denom)

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
               levels = c("p11.full", "p00.full", "p11.p00"),
               labels = c("p11 - full", "p00 - full", "p11 - p00"))


    smd$variable <- ifelse(smd$variable %in% standardize,
                           paste0("*", smd$variable),
                           smd$variable)

    vars <- ifelse(vars %in% standardize,
                   paste0("*", vars),
                   vars)

    smd$variable <- factor(smd$variable, levels = vars)


    smd[,c("variable", "contrast.type", "contrast", "mean.diff")]

}























































