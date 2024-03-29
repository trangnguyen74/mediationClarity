



#### OK  estimate_Ypred ########################################################

#' Estimator Ypred
#'
#' Function that implements estimator Ypred
#' @inheritParams estimate_wtd
#' @inheritParams weights_med
#' @inheritParams estimate_psYpred
#' @param cm.std blah
#' @param cm.order blah
#' @param y10.c.form Model formula for E[Y(1,M0)|C]
#' @param y01.c.form Model formula for E[Y(0,M1)|C]
#' @family estimators
#' @seealso [estimate_YpredMR()] which is closely related but more robust
#' @export

estimate_Ypred <- function(
    data,
    s.wt.var     = NULL,
    cross.world  = "10",
    effect.scale = "MD",

    boot.num      = 999,
    boot.seed     = NULL,
    boot.method   = "cont-wt",
    boot.stratify = TRUE,

    a.cm.form,
    max.stabilized.wt = 30,

    plot        = TRUE,
    cm.std = NULL,
    cm.order    = NULL,

    y.c.form   = NULL,
    y.c1.form  = NULL,
    y.c0.form  = NULL,
    y10.c.form = NULL,
    y01.c.form = NULL,
    y.link     = "identity"
) {


    # CLEAN INPUTS

    cm.vars <- c.vars <- m.vars <- y.family <- NULL

    .prep_Ypred()

    key.inputs <- mget(c("cross.world",
                         "effect.scale",
                         "a.cm.form",
                         "max.stabilized.wt",
                         "y.c1.form",
                         "y.c0.form",
                         "y10.c.form",
                         "y01.c.form",
                         "y.family"))






    # POINT ESTIMATION

    if (!plot) {

        estimates <- do.call(".point_est.Ypred", c(key.inputs,
                                                   list(data        = data,
                                                        output.data = FALSE)))
    } else {

        tmp <- do.call(".point_est.Ypred", c(key.inputs,
                                             list(data        = data,
                                                  output.data = TRUE)))
        estimates <- tmp$estimates

        plots <- .plot_odds(w.dat          = tmp$w.dat,
                            vars           = cm.vars,
                            vars.std       = cm.std,
                            estimate.Ypred = TRUE);  rm(tmp)
    }


    # BOOTSTRAP

    if (boot.num > 0) {
        ci.se <- .boot_ci.se(data       = data,
                             stratify   = boot.stratify,
                             boot.num   = boot.num,
                             seed       = boot.seed,
                             method     = boot.method,
                             FUN        = ".point_est.Ypred",
                             FUN.inputs = key.inputs)

        estimates <- cbind(estimate = estimates,
                           ci.se)
        rm(ci.se)
    }



    # OUTPUT

    if (!plot && boot.num==0) return(estimates)

    out <- list(estimates = estimates)

    if (boot.num > 0)  out$boot.seed <- boot.seed
    if (plot)          out$plots     <- plots

    out

}




#### OK  .prep_Ypred ###########################################################

#' @rdname dot-prep
#' @order 9

.prep_Ypred <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_effect.scale(top.env)

    .clean_boot(top.env)

    .clean_weights.Ypred(top.env)

    .clean_y.Ypred(top.env)

    if (top.env$plot) .check_plot.Ypred(top.env)
}




#### .clean_weights.Ypred ####################################################

#' @rdname dot-clean_weights
#' @order 6
#' @details \code{.clean_weights.Ypred()} is used by \code{.prep_Ypred()}.

.clean_weights.Ypred <- function(env) {

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



#### OK  .clean_y.Ypred ##################################################

#' @rdname dot-clean_y
#' @order 5

.clean_y.Ypred <- function(env) {

    yes10 <- ("10" %in% env$cross.world)
    yes01 <- ("01" %in% env$cross.world)

    cm.vars <- env$cm.vars

    y.c <- env$y.c.form
    y.c1 <- env$y.c1.form
    y.c0 <- env$y.c0.form
    y10.c <- env$y10.c.form
    y01.c <- env$y01.c.form


    if (is.null(y.c)) {

        if (is.null(y.c1))
            stop("Must specify either y.c1.form or y.c.form.")

        if (is.null(y.c0))
            stop("Must specify either y.c0.form or y.c.form.")

        if (yes10 && is.null(y10.c))
            stop("For cross.world==\"10\", must specify either y10.c.form or y.c.form.")

        if (yes01 && is.null(y01.c))
            stop("For cross.world==\"01\", must specify either y01.c.form or y.c.form.")

    } else {

        if (is.null(env$y.c1.form))  y.c1 <- y.c

        if (is.null(env$y.c0.form))  y.c0 <- y.c

        if (yes10 && is.null(y10.c))  y10.c <- y.c

        if (yes01 && is.null(y01.c))  y01.c <- y.c
    }

    env$y.c1.form <- y.c1
    env$y.c0.form <- y.c0
    env$y10.c.form <- y10.c
    env$y01.c.form <- y01.c






    y.var  <- unique(c(all.vars(formula(y.c1)[[2]]),
                       all.vars(formula(y.c0)[[2]])))

    if (yes10)  y.var <- unique(c(y.var, all.vars(formula(y10.c)[[2]])))
    if (yes01)  y.var <- unique(c(y.var, all.vars(formula(y01.c)[[2]])))

    if (length(y.var)>1)
        stop("Outcome variable is not unique across outcome models.")

    env$data$.y <- env$data[, y.var]



    c.vars <- unique(c(all.vars(formula(y.c1)[[3]]),
                       all.vars(formula(y.c0)[[3]])))

    if (yes10)  c.vars <- unique(c(c.vars, all.vars(formula(y10.c)[[3]])))
    if (yes01)  c.vars <- unique(c(c.vars, all.vars(formula(y01.c)[[3]])))

    stray.c <- setdiff(c.vars, cm.vars)

    if (length(stray.c)>0)
        stop(paste("Covariate(s)", paste(stray.c, collapse = ", "), "in outcome model(s) not found in a.c.form."))



    m.vars <- setdiff(cm.vars, c.vars)
    c.vars <- setdiff(cm.vars, m.vars)

    if (length(m.vars)==0)
        stop("Based on outcome given covariates model formula(s), all predictors in a.cm.form are covariates, and none are mediators.")

    env$c.vars <- c.vars
    env$m.vars <- m.vars





    y.link <- env$y.link

    if (!(y.link %in% c("identity", "logit", "logistic", "log")))
        stop("y.link not recognized or supported. Supported options include: \"identity\" (for linear model with numeric outcome), \"logit\" (for binary outcome or outcome bounded in (0,1) interval), and \"log\" (for non-negative outcome).")

    if (length(unique(env$data$.y))==2) {

        if (!(y.link %in% c("logit", "logistic")))
            warning("The outcome is binary. Logit model is used.")

        env$y.family <- "quasibinomial"

    } else if (is.numeric(env$data$.y) && all(env$data$.y>=0)) {

        if (y.link=="log") {
            env$y.family <- "quasipoisson"

        } else if (y.link=="identity") {
            env$y.family <- "gaussian"

        } else if (all(env$data$.y<=1)) {
            env$y.family <- "quasibinomial"

        } else {
            warning("Logit link not allowed for outcome that is non-binary and not bounded in the (0,1) interval. Identity link is used instead. May also consider log link.")
            env$y.family <- "gaussian"
        }
    } else if (is.numeric(env$data$.y)) {
        if (!y.link=="identity")
            warning("Outcome is numeric variable with negative values. Identity link is used.")
        env$y.family <- "gaussian"

    } else
        stop("Outcome type not supported.")
}




#### OK  .check_plot.Ypred #####################################################

#' @rdname dot-check_plot
#' @order 4
#' @details \code{.check_plot.Ypred()} is called by \code{.prep_Ypred()}.

.check_plot.Ypred <- function(env) {

    if (is.null(env$cm.order)) env$cm.order <- c(env$c.vars, env$m.vars)


    cm.vars  <- env$cm.vars
    cm.order <- env$cm.order
    cm.std   <- env$cm.std


    if (!is.null(cm.order)) {

        if (!setequal(cm.vars, cm.order)) {
            warning("Variables in cm.order do not match covariates from a.c.form. Ignoring c.order.")
        } else
            env$cm.vars <- cm.vars <- cm.order
    }




    if (is.null(cm.std)) {

        maybe.cont <- sapply(cm.std,
                             function(z) maybe_continuous(env$data[, z]))

        if (any(maybe.cont))
            message(paste("Consider whether the balance plot should use standardized mean differences for numeric covariate/mediators",
                          paste(cm.std[which(maybe.cont)],
                                collapse = ", "),
                          "(if they are continuous variables).",
                          "To turn off this message, specify cm.std=\"\"."))

        return()
    }


    if (length(cm.std)==1 && cm.std=="")  return()



    cm.std <- setdiff(cm.std, "")

    if (length(setdiff(cm.std, cm.vars))>0)
        stop("Variables specified in cm.std are not all contained in model formula a.cm.form.")



    ok.std <- sapply(cm.std, function(z) maybe_continuous(env$data[, z]))

    if (!all(ok.std))
        stop(paste("Check variable(s)",
                   paste(cm.std[which(!ok.std)], collapse = ", "),
                   "before proceeding. Only include continuous variables in cm.std."))
    cm.vars  <- env$cm.vars
    cm.order <- env$cm.order
    cm.std   <- env$cm.std


    if (!is.null(cm.order)) {

        if (!setequal(cm.vars, cm.order)) {
            warning("Variables in cm.order do not match covariates from a.c.form. Ignoring c.order.")
        } else
            env$cm.vars <- cm.vars <- cm.order
    }




    if (is.null(cm.std)) {

        maybe.cont <- sapply(cm.std,
                             function(z) maybe_continuous(env$data[, z]))

        if (any(maybe.cont))
            message(paste("Consider whether the balance plot should use standardized mean differences for numeric covariate/mediators",
                          paste(cm.std[which(maybe.cont)],
                                collapse = ", "),
                          "(if they are continuous variables).",
                          "To turn off this message, specify cm.std=\"\"."))

        return()
    }


    if (length(cm.std)==1 && cm.std=="")  return()



    cm.std <- setdiff(cm.std, "")

    if (length(setdiff(cm.std, cm.vars))>0)
        stop("Variables specified in cm.std are not all contained in model formula a.cm.form.")



    ok.std <- sapply(cm.std, function(z) maybe_continuous(env$data[, z]))

    if (!all(ok.std))
        stop(paste("Check variable(s)",
                   paste(cm.std[which(!ok.std)], collapse = ", "),
                   "before proceeding. Only include continuous variables in cm.std."))

}




#### OK  .point_est.Ypred ######################################################

#' @rdname dot-point_est
#' @order 4
#' @param y10.c.form,y01.c.form Model formulas for E[Y(1,M0)|C] and E[Y(0,M1)|C].

.point_est.Ypred <- function(
    data,
    cross.world,
    effect.scale,
    a.cm.form,
    max.stabilized.wt = 30,
    output.data = FALSE,
    y.c1.form,
    y.c0.form,
    y10.c.form = NULL,
    y01.c.form = NULL,
    y.family
) {

    w.dat <- .compute_weights.odds(data              = data,
                                      cross.world       = cross.world,
                                      a.form         = a.cm.form,
                                      max.stabilized.wt = max.stabilized.wt
    )

    s11 <- w.dat[w.dat$.samp=="s11", ]
    s00 <- w.dat[w.dat$.samp=="s00", ]

    full <- rbind(s11, s00)

    y.c1.s11 <- glm(formula = y.c1.form,
                    data = s11,
                    weights = data$.s.wt,
                    family = y.family)

    y.c0.s00 <- glm(formula = y.c0.form,
                    data = s00,
                    weights = data$.s.wt,
                    family = y.family)

    pred <- full[".f.wt"]
    pred$p00 <- predict(y.c0.s00, newdata = full, type = "response")
    pred$p11 <- predict(y.c1.s11, newdata = full, type = "response")


    if ("10" %in% cross.world) {

        y10.c.s10 <- glm(formula = y10.c.form,
                         data = w.dat[w.dat$.samp=="s10", ],
                         weights = data$.f.wt,
                         family = y.family)

        pred$p10 <- predict(y10.c.s10, newdata = full, type = "response")
    }


    if ("01" %in% cross.world) {

        y01.c.s01 <- glm(formula = y01.c.form,
                         data = w.dat[w.dat$.samp=="s01", ],
                         weights = data$.f.wt,
                         family = y.family)

        pred$p01 <- predict(y01.c.s01, newdata = full, type = "response")
    }


    pred <- .reshape_gather(pred,
                           columns = setdiff(colnames(pred), ".f.wt"),
                           key = ".samp",
                           value = ".y",
                           wide.row = FALSE)


    estimates <- .get_means.and.effects(w.dat = pred,
                                        effect.scale = effect.scale)


    if (!output.data) return(estimates)

    list(estimates = estimates,
         w.dat    = w.dat)

}










