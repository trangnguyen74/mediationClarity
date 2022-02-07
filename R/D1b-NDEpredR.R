#### estimate_NDEpredR ########################################################

#' Estimator NDEpredR
#'
#' Function that implements estimator NDEYpred.R
#' @inheritParams estimate_Y2predR
#' @inheritParams estimate_NDEpred
#' @family estimators
#' @family additive-effect estimators
#' @family more-robust estimators
#' @export

estimate_NDEpredR <- function(
    data,
    s.wt.var     = NULL,
    cross.world  = "10",

    boot.num      = 999,
    boot.seed     = NULL,
    boot.method   = "cont-wt",
    boot.stratify = TRUE,

    a.c.form,
    a.cm.form,
    max.stabilized.wt = 30,

    plot = TRUE,
    c.vars.std = NULL,
    m.vars.std = NULL,
    c.order    = NULL,
    m.order    = NULL,

    # detailed outcome/effect models
    y.c1.form = NULL,
    y.c0.form = NULL,

    y.cm1.form = NULL,
    y.cm0.form = NULL,

    nde0.c.form = NULL,
    nde1.c.form = NULL,

    # or shortcuts
    y.c.form  = NULL,
    y.cm.form  = NULL,
    nde.c.form  = NULL,

    y.link = "identity"

) {

    # CLEAN INPUTS

    c.vars <- m.vars <- y.family <- NULL

    .prep_NDEpredR()

    key.inputs <- mget(c("cross.world",
                         "a.c.form", "a.cm.form",
                         "max.stabilized.wt",
                         "y.c1.form", "y.c0.form",
                         "y.cm1.form", "y.cm0.form",
                         "nde0.c.form", "nde1.c.form",
                         "y.family"))


    # POINT ESTIMATION

    if (!plot) {

        estimates <- do.call(".point_est.NDEpredR",
                             c(key.inputs, list(data        = data,
                                                output.data = FALSE)))
    } else {

        tmp <- do.call(".point_est.NDEpredR",
                       c(key.inputs, list(data        = data,
                                          output.data = TRUE)))

        estimates <- tmp$estimates

        plots <- .plot_med(w.dat      = tmp$w.dat,
                           c.vars     = c.vars,
                           m.vars     = m.vars,
                           c.vars.std = c.vars.std,
                           m.vars.std = m.vars.std);  rm(tmp)
    }



    # BOOTSTRAP

    if (boot.num > 0) {
        ci.se <- .boot_ci.se(data       = data,
                             stratify   = boot.stratify,
                             boot.num   = boot.num,
                             seed       = boot.seed,
                             method     = boot.method,
                             FUN        = ".point_est.NDEpredR",
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




#### .prep_NDEpredR ##########################################################

#' @rdname dot-prep
#' @order 11

.prep_NDEpredR <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_boot(top.env)

    .clean_weights.med(top.env)

    .clean_ye.forms.NDEpred(top.env)

    if (top.env$plot) .check_plot.med(top.env)
}




#### .point_est.NDEpredR #####################################################

#' @rdname dot-point_est
#' @order 8

.point_est.NDEpredR <- function(
    data,
    cross.world,
    a.c.form,
    a.cm.form,
    max.stabilized.wt = 30,
    output.data = FALSE,
    y.c1.form,
    y.c0.form,
    y.cm1.form,
    y.cm0.form,
    nde0.c.form,
    nde1.c.form,
    y.family
) {



    w.dat <- .compute_weights.med(
        data              = data,
        cross.world       = cross.world,
        a.c.form          = a.c.form,
        a.cm.form         = a.cm.form,
        max.stabilized.wt = max.stabilized.wt
    )


    full <- w.dat[w.dat$.samp=="p11" | w.dat$.samp=="p00", ]
    full$.w.wt <- 1
    full$.f.wt <- full$.s.wt



    y.c1.p11 <- glm(formula = y.c1.form,
                    data    = w.dat[w.dat$.samp=="p11", ],
                    weights = data$.f.wt,
                    family  = y.family)

    y.c0.p00 <- glm(formula = y.c0.form,
                    data    = w.dat[w.dat$.samp=="p00", ],
                    weights = data$.f.wt,
                    family  = y.family)

    pred.te <-
        predict(y.c1.p11, newdata = full, type = "response") -
        predict(y.c0.p00, newdata = full, type = "response")


    estimates <- list(TE = .wtd_mean(pred.te, data$.s.wt))


    if ("10" %in% cross.world) {

        y.cm1.p10 <- glm(formula = y.cm1.form,
                         data    = w.dat[w.dat$.samp=="p10", ],
                         weights = data$.f.wt,
                         family  = y.family)

        p00 <- w.dat[w.dat$.samp=="p00", ]

        p00$nde0.prox <-
            predict(y.cm1.p10, newdata = p00, type = "response") - p00$.y


        nde0.c.form <- paste("nde0.prox ~",
                             paste(all.vars(formula(nde0.c.form)[[3]]),
                                   collapse = " + "))

        if (y.family=="quasibinomial")
            p00$nde0.prox <- (p00$nde0.prox + 1) / 2

        nde0.c.p00 <- glm(formula = nde0.c.form,
                          data    = p00,
                          weights = data$.f.wt,
                          family  = y.family)

        pred.nde0 <- predict(nde0.c.p00, newdata = full, type = "response")


        if (y.family=="quasibinomial")
            pred.nde0 <- pred.nde0 * 2 - 1

        estimates$NDE0 <- .wtd_mean(pred.nde0, full$.s.wt)
        estimates$NIE1 <- estimates$TE - estimates$NDE0

    }



    if ("01" %in% cross.world) {

        y.cm0.p01 <- glm(formula = y.cm0.form,
                         data    = w.dat[w.dat$.samp=="p01", ],
                         weights = data$.s.wt,
                         family  = y.family)

        p11 <- w.dat[w.dat$.samp=="p11", ]

        p11$nde1.prox <-
            p11$.y - predict(y.cm0.p01, newdata = p11, type = "response")


        nde1.c.form <- paste("nde1.prox ~",
                             paste(all.vars(formula(nde1.c.form)[[3]]),
                                   collapse = " + "))

        if (y.family=="quasibinomial")
            p11$nde1.prox <- (p11$nde1.prox + 1) / 2

        nde1.c.p11 <- glm(formula = nde1.c.form,
                          data    = p11,
                          weights = data$.f.wt,
                          family  = y.family)

        pred.nde1 <- predict(nde1.c.p11, newdata = full, type = "response")


        if (y.family=="quasibinomial")
            pred.nde1 <- pred.nde1 * 2 - 1

        estimates$NDE1 <- .wtd_mean(pred.nde1, data$.s.wt)
        estimates$NIE0 <- estimates$TE - estimates$NDE1

    }

    estimates <- unlist(estimates)


    if (!output.data) return(estimates)

    list(estimates = estimates,
         w.dat     = w.dat)
}
















































































