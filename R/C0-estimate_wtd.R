



#### OK  estimate_wtd ##########################################################

#' Pure weighting estimator
#'
#' Function that implements the pure weighting estimator of natural (in)direct effects.
#' @inheritParams estimate_effects
#' @inheritParams weights_med
#' @param y.var Name of the outcome variable.
#' @return If \code{plot==FALSE} and \code{boot.num==0}, a point estimate vector including the potential outcome means and the effects.
#' @return Otherwise, a list of objects including:\tabular{ll}{
#' \code{estimates} \tab point estimates (and if \code{boot.num>0}) 95\% quantile intervals and standard errors \cr
#' \cr
#' \code{boot.seed} \tab the seed used for the bootstrap (if \code{boot.num>0}) \cr
#' \cr
#' \code{wt.dist.plot} \tab a \code{ggplot2} plot of weight distributions (if \code{plot==TRUE}) \cr
#' \cr
#' \code{balance.plot} \tab a \code{ggplot2} plot of balance on the means of covariate and mediator for the pseudo samples (if \code{plot==TRUE}) \cr
#' }
#' @family estimators
#' @export

estimate_wtd <- function(
    data,
    s.wt.var     = NULL,
    cross.world  = "10",
    effect.scale = "MD",

    boot.num      = 999,
    boot.seed     = NULL,
    boot.method   = "cont-wt",
    boot.stratify = TRUE,

    a.c.form,
    a.cm.form,
    y.var,
    max.stabilized.wt = 30,

    plot = TRUE,
    c.std = NULL,
    m.std = NULL,
    c.order = NULL,
    m.order = NULL) {


    # CLEAN INPUTS

    c.vars <- m.vars <- NULL

    .prep_estimate_wtd()


    key.inputs <- mget(c("cross.world",
                         "effect.scale",
                         "a.c.form",
                         "a.cm.form",
                         "max.stabilized.wt"))




    # POINT ESTIMATION

    if (!plot) {

        estimates <- do.call(".point_est.wtd", c(key.inputs,
                                                 list(data        = data,
                                                      output.data = FALSE)))
    } else {

        tmp <- do.call(".point_est.wtd", c(key.inputs,
                                           list(data        = data,
                                                output.data = TRUE)))
        estimates <- tmp$estimates

        plots <- .plot_med(w.dat = tmp$w.dat,
                           c.vars = c.vars,
                           m.vars = m.vars,
                           c.std = c.std,
                           m.std = m.std,
                           key.balance = TRUE);     rm(tmp)
    }


    # BOOTSTRAP

    if (boot.num > 0) {
        ci.se <- .boot_ci.se(data       = data,
                             stratify   = boot.stratify,
                             boot.num   = boot.num,
                             seed       = boot.seed,
                             method     = boot.method,
                             FUN        = ".point_est.wtd",
                             FUN.inputs = key.inputs)

        estimates <- cbind(estimate = estimates,
                           ci.se)
        rm(ci.se)
    }



    # OUTPUT

    if (!plot && boot.num==0)  return(estimates)


    out <- list(estimates = estimates)

    if (boot.num > 0)  out$boot.seed <- boot.seed
    if (plot)          out$plots     <- plots

    out
}






#### OK  .prep_estimate_wtd ####################################################

#' @rdname dot-prep
#' @order 4
.prep_estimate_wtd <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_effect.scale(top.env)

    .clean_boot(top.env)

    .clean_weights.med(top.env)

    .clean_y.wtd(top.env)

    if (top.env$plot) .check_plot.med(top.env)
}



#### OK  .clean_y.wtd ##########################################################

#' (For maintainer) Clean outcome-related inputs
#'
#' Internal functions called by \code{estimate_} functions to clean inputs related to the outcome including \code{y.var}, \code{y.c.form}, \code{y.cm.form}, etc.
#' @name dot-clean_y
#' @keywords internal
NULL

#' @rdname dot-clean_y
#' @order 1
#' @inheritParams env-block

.clean_y.wtd <- function(env) {

    y.var <- env$y.var

    if (!y.var %in% names(env$data))
        stop(paste("Variable", y.var, "(y.var) not found in dataset."))

    env$data$.y <- env$data[, y.var]
}




#### OK  .point_est.wtd ########################################################

#' (For maintainer) Point estimation
#'
#' Internal functions to be called within \code{estimate_\*} to obtain the point estimate of effects.
#' @param data A dataset that has been prepared to have the sampling with variable named s.wt (e.g., using function \code{.clean_inputs_generic()}.
#' @param cross.world Three options: "10", "01" or c("10", "01").
#' @param effect.scale Three options: "MD", "MR", "OR".
#' @inheritParams weights_med
#' @param output.data Whether to output the weighted data in addition to the estimated potential outcome means and effects. Defaults to FALSE.
#' @return If \code{plot==FALSE}, effect estimates (numeric vector). If \code{plot==TRUE}, a list including this vector of effect estimates and the weight distribution and mean balance plots.
#' @name dot-point_est
#' @keywords internal
NULL

#' @rdname dot-point_est
#' @order 1

.point_est.wtd <- function(
    data,
    cross.world,
    effect.scale,
    a.c.form,
    a.cm.form,
    max.stabilized.wt = 30,
    output.data = FALSE
) {

    w.dat <- .compute_weights.med(
        data              = data,
        cross.world       = cross.world,
        a.c.form          = a.c.form,
        a.cm.form         = a.cm.form,
        max.stabilized.wt = max.stabilized.wt
    )

    estimates <- .get_means.and.effects(w.dat = w.dat,
                                        effect.scale = effect.scale)


    if (!output.data) return(estimates)

    list(estimates = estimates,
         w.dat    = w.dat)

}






