



#### OK  estimate_NDEpred ######################################################

#' Estimator NDEpred
#'
#' Function that implements estimator NDEYpred
#' @inheritParams estimate_Y2pred
#' @param nde0.c.form Model formula for E[NDE0|C]. Use any name for the response model in this formula, e.g., "nde.prox" or "effect".
#' @param nde1.c.form Model formula for E[NDE1|C]. Use any name for response model in this formula, e.g., "nde.prox" or "effect".
#' @param nde.c.form Shortcut to specify the same formula for both E[NDE0|C] and E[NDE1|C] models.
#' @family estimators
#' @family additive-effect estimators
#' @export

estimate_NDEpred <- function(
    data,
    s.wt.var     = NULL,
    cross.world  = "10",

    boot.num      = 999,
    boot.seed     = NULL,
    boot.method   = "cont-wt",
    boot.stratify = TRUE,

    a.var,

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

    .prep_NDEpred()

    key.inputs <- mget(c("cross.world",
                         "y.c1.form", "y.c0.form",
                         "y.cm1.form", "y.cm0.form",
                         "nde0.c.form", "nde1.c.form",
                         "y.family"))


    # POINT ESTIMATION

    estimates <- do.call(".point_est.NDEpred", c(key.inputs,
                                                 list(data = data)))


    # BOOTSTRAP

    if (boot.num==0) return(estimates)


    ci.se <- .boot_ci.se(data       = data,
                         stratify   = boot.stratify,
                         boot.num   = boot.num,
                         seed       = boot.seed,
                         method     = boot.method,
                         FUN        = ".point_est.NDEpred",
                         FUN.inputs = key.inputs)

    estimates <- cbind(estimate = estimates,
                       ci.se)


    mget(c("estimates", "boot.seed"))


}




#### .prep_NDEpred ##########################################################

#' @rdname dot-prep
#' @order 10

.prep_NDEpred <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_boot(top.env)

    .clean_a.var(top.env)

    .clean_ye.forms.NDEpred(top.env)
}




#### .clean_ye.forms.NDEpred ##############################################

#' @rdname dot-clean_y
#' @order 7

.clean_ye.forms.NDEpred <- function(env) {

    yes10 <- ("10" %in% env$cross.world)
    yes01 <- ("01" %in% env$cross.world)

    y.c    <- env$y.c.form
    y.c1   <- env$y.c1.form
    y.c0   <- env$y.c0.form
    y.cm   <- env$y.cm.form
    y.cm1  <- env$y.cm1.form
    y.cm0  <- env$y.cm0.form
    nde.c  <- env$nde.c.form
    nde0.c <- env$nde0.c.form
    nde1.c <- env$nde1.c.form
    y.link <- env$y.link



    if (is.null(y.c)) {

        if (is.null(y.c1))
            stop("Must specify either y.c1.form or y.c.form.")

        if (is.null(y.c0))
            stop("Must specify either y.c0.form or y.c.form.")

    } else {

        if (is.null(y.c1))  env$y.c1.form <- y.c1 <- y.c

        if (is.null(y.c0))  env$y.c0.form <- y.c0 <- y.c
    }



    if (is.null(y.cm)) {

        if (yes10 && is.null(env$y.cm1.form))
            stop("Must specify either y.cm1.form or y.cm.form.")

        if (yes01 && is.null(env$y.cm0.form))
            stop("Must specify either y.cm0.form or y.cm.form.")

    } else {

        if (yes10 && is.null(env$y.cm1.form))  env$y.cm1.form <- y.cm1 <- y.cm

        if (yes01 && is.null(env$y.cm0.form))  env$y.cm0.form <- y.cm0 <- y.cm
    }



    if (is.null(nde.c)) {

        if (yes10 && is.null(nde0.c))
            stop("Must specify either nde.c.form or nde0.c.form.")

        if (yes01 && is.null(nde1.c))
            stop("Must specify either nde.c.form or nde1.c.form.")

    } else {

        if (yes10 && is.null(nde0.c))  env$nde0.c.form <- nde0.c <- y.c

        if (yes01 && is.null(nde1.c))  env$nde1.c.form <- nde1.c <- y.c
    }




    y.var <- unique(all.vars(formula(y.c1)[[2]]),
                    all.vars(formula(y.c0)[[2]]))

    c.vars <- unique(all.vars(formula(y.c1)[[3]]),
                     all.vars(formula(y.c0)[[3]]))

    cm.vars <- NULL

    if (yes10) {
        y.var <- unique(c(y.var,
                          all.vars(formula(y.cm1)[[2]])))
        c.vars <- unique(c(c.vars,
                           all.vars(formula(nde0.c)[[3]])))
        cm.vars <- unique(c(cm.vars,
                            all.vars(formula(y.cm1)[[3]])))
    }

    if (yes01) {
        y.var <- unique(c(y.var,
                          all.vars(formula(y.cm0)[[2]])))
        c.vars <- unique(c(c.vars,
                           all.vars(formula(nde1.c)[[3]])))
        cm.vars <- unique(c(cm.vars,
                            all.vars(formula(y.cm0)[[3]])))
    }

    if (length(y.var)>1)
        stop("Outcome variable is not unique across outcome models.")

    env$data$.y <- env$data[, y.var]


    if (all(cm.vars %in% c.vars))
        stop("The combination of outcome and effect models imply all predictors are covariates and there are no mediators. Please double-check!")





    if (!(y.link %in% c("identity", "logit", "logistic", "log")))
        stop("y.link not recognized or supported. Supported options include: \"identity\" (for linear model with numeric outcome), \"logit\" (for binary outcome or outcome bounded in (0,1) interval), and \"log\" (for non-negative outcome).")

    if (length(unique(env$data$.y))==2) {

        if (!(y.link %in% c("logit", "logistic")))
            warning("The outcome is binary. Logit model is used.")

        env$y.family <- "quasibinomial"

    } else if (is.numeric(env$data$.y) && all(env$data$.y>=0)) {

        if (y.link=="log") {
            stop("Please use another estimator. NDEpred is only implemented with linear/logit outcome models.")

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




#### .point_est.NDEpred #####################################################

#' @rdname dot-point_est
#' @order 7

.point_est.NDEpred <- function(
    data,
    cross.world,
    y.c1.form,
    y.c0.form,
    y.cm1.form,
    y.cm0.form,
    nde0.c.form,
    nde1.c.form,
    y.family,
    output.data = FALSE # this is to work nice with boot function
    #                   # (to revisit later)
) {

    s11 <- data[data$.a==1, ]
    s00 <- data[data$.a==0, ]

    y.c1.s11 <- glm(formula = y.c1.form,
                    data    = s11,
                    weights = data$.s.wt,
                    family  = y.family)

    y.c0.s00 <- glm(formula = y.c0.form,
                    data    = s00,
                    weights = data$.s.wt,
                    family  = y.family)

    pred.te <-
        predict(y.c1.s11, newdata = data, type = "response") -
        predict(y.c0.s00, newdata = data, type = "response")


    out <- list(TE = .wtd_mean(pred.te, data$.s.wt))


    if ("10" %in% cross.world) {

        y.cm1.s11 <- glm(formula = y.cm1.form,
                         data    = s11,
                         weights = data$.s.wt,
                         family  = y.family)

        s00.tmp <- s00

        s00.tmp$nde0.prox <-
            predict(y.cm1.s11, newdata = s00, type = "response") - s00$.y


        nde0.c.form <- paste("nde0.prox ~",
                             paste(all.vars(formula(nde0.c.form)[[3]]),
                                   collapse = " + "))

        if (y.family=="quasibinomial")
            s00.tmp$nde0.prox <- (s00.tmp$nde0.prox + 1) / 2

        nde0.c.s00 <- glm(formula = nde0.c.form,
                          data    = s00.tmp,
                          weights = data$.s.wt,
                          family  = y.family)

        pred.nde0 <- predict(nde0.c.s00, newdata = data, type = "response")


        if (y.family=="quasibinomial")
            pred.nde0 <- pred.nde0 * 2 - 1

        out$NDE0 <- .wtd_mean(pred.nde0, data$.s.wt)
        out$NIE1 <- out$TE - out$NDE0

    }



    if ("01" %in% cross.world) {

        y.cm0.s00 <- glm(formula = y.cm0.form,
                         data    = s00,
                         weights = data$.s.wt,
                         family  = y.family)

        s11.tmp <- s11

        s11.tmp$nde1.prox <-
            s11$.y - predict(y.cm0.s00, newdata = s11, type = "response")


        nde1.c.form <- paste("nde1.prox ~",
                             paste(all.vars(formula(nde1.c.form)[[3]]),
                                   collapse = " + "))

        if (y.family=="quasibinomial")
            s11.tmp$nde1.prox <- (s11.tmp$nde1.prox + 1) / 2

        nde1.c.s11 <- glm(formula = nde1.c.form,
                          data    = s11.tmp,
                          weights = data$.s.wt,
                          family  = y.family)

        pred.nde1 <- predict(nde1.c.s11, newdata = data, type = "response")


        if (y.family=="quasibinomial")
            pred.nde1 <- pred.nde1 * 2 - 1

        out$NDE1 <- .wtd_mean(pred.nde1, data$.s.wt)
        out$NIE0 <- out$TE - out$NDE1

    }



    unlist(out)
}
















































































