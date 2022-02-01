



#### OK  estimate_Y2pred #######################################################

#' Estimator Y2pred
#'
#' Function that implements estimator Y2pred
#' @inheritParams estimate_psYpred
#' @inheritParams estimate_Ypred
#' @param a.var Name of treatment variable.
#' @family estimators
#' @export

estimate_Y2pred <- function(
    data,
    s.wt.var     = NULL,
    cross.world  = "10",
    effect.scale = "additive",

    boot.num      = 999,
    boot.seed     = NULL,
    boot.method   = "cont-wt",
    boot.stratify = TRUE,

    a.var,
    y.cm.form = NULL,
    y.cm1.form = NULL,
    y.cm0.form = NULL,
    y.c.form = NULL,
    y.c1.form = NULL,
    y.c0.form = NULL,
    y10.c.form = NULL,
    y01.c.form = NULL,
    y.link = "identity") {


    # CLEAN INPUTS

    y.family <- NULL

    .prep_Y2pred()

    key.inputs <- mget(c("cross.world",
                         "effect.scale",
                         "y.c1.form",
                         "y.c0.form",
                         "y.cm1.form",
                         "y.cm0.form",
                         "y10.c.form",
                         "y01.c.form",
                         "y.family"))


    # POINT ESTIMATION

    estimates<- do.call(".point_est.Y2pred", c(key.inputs,
                                               list(data = data)))


    # BOOTSTRAP

    if (boot.num==0) return(estimates)


    ci.se <- .boot_ci.se(data       = data,
                         stratify   = boot.stratify,
                         boot.num   = boot.num,
                         seed       = boot.seed,
                         method     = boot.method,
                         FUN        = ".point_est.Y2pred",
                         FUN.inputs = key.inputs)

    estimates <- cbind(estimate = estimates,
                       ci.se)


    mget(c("estimates", "boot.seed"))

}




#### OK  .prep_Y2pred ##########################################################

#' @rdname dot-prep
#' @order 9

.prep_Y2pred <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_effect.scale(top.env)

    .clean_boot(top.env)

    .clean_a.var(top.env)

    .clean_y.forms.Y2pred(top.env)
}




#### OK  .clean_a.var ##########################################################

#' @rdname dot-clean_weights
#' @order 4

.clean_a.var <- function(env) {

    a.var <- env$a.var

    if (!a.var %in% names(env$data))
        stop(paste0("Treatment variable (",
                   a.var,
                   ") not found in dataset."))

    env$data$.a <- env$data[, a.var]
}




#### OK  .clean_y.forms.Y2pred #################################################

#' @rdname dot-clean_y
#' @order 6

.clean_y.forms.Y2pred <- function(env) {

    yes10 <- ("10" %in% env$cross.world)
    yes01 <- ("01" %in% env$cross.world)

    y.c    <- env$y.c.form
    y.c1   <- env$y.c1.form
    y.c0   <- env$y.c0.form
    y10.c  <- env$y10.c.form
    y01.c  <- env$y01.c.form
    y.cm   <- env$y.cm.form
    y.cm1  <- env$y.cm1.form
    y.cm0  <- env$y.cm0.form
    y.link <- env$y.link



    if (is.null(y.cm)) {

        if (yes10 && is.null(env$y.cm1.form))
            stop("Must specify either y.cm1.form or y.cm.form.")

        if (yes01 && is.null(env$y.cm0.form))
            stop("Must specify either y.cm0.form or y.cm.form.")

    } else {

        if (yes10 && is.null(env$y.cm1.form))  env$y.cm1.form <- y.cm1 <- y.cm

        if (yes01 && is.null(env$y.cm0.form))  env$y.cm0.form <- y.cm0 <- y.cm
    }



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

        if (is.null(y.c1))  env$y.c1.form <- y.c1 <- y.c

        if (is.null(y.c0))  env$y.c0.form <- y.c0 <- y.c


        if (yes10 && is.null(y10.c))  env$y10.c.form <- y10.c <- y.c

        if (yes01 && is.null(y01.c))  env$y01.c.form <- y01.c <- y.c
    }




    y.var <- unique(all.vars(formula(y.c1)[[2]]),
                    all.vars(formula(y.c0)[[2]]))

    c.vars <- unique(all.vars(formula(y.c1)[[3]]),
                     all.vars(formula(y.c0)[[3]]))

    cm.vars <- NULL

    if (yes10) {
        y.var <- unique(c(y.var,
                          all.vars(formula(y10.c)[[2]]),
                          all.vars(formula(y.cm1)[[2]])))
        c.vars <- unique(c(c.vars,
                           all.vars(formula(y10.c)[[3]])))
        cm.vars <- unique(c(cm.vars,
                            all.vars(formula(y.cm1)[[3]])))
    }

    if (yes01) {
        y.var <- unique(c(y.var,
                          all.vars(formula(y01.c)[[2]]),
                          all.vars(formula(y.cm0)[[2]])))
        c.vars <- unique(c(c.vars,
                           all.vars(formula(y01.c)[[3]])))
        cm.vars <- unique(c(cm.vars,
                            all.vars(formula(y.cm0)[[3]])))
    }

    if (length(y.var)>1)
        stop("Outcome variable is not unique across outcome models.")

    env$data$.y <- env$data[, y.var]


    if (all(cm.vars %in% c.vars))
        stop("The combination of outcome models imply all predictors are covariates and there are no mediators. Please double-check!")





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




#### OK  .point_est.Y2pred #####################################################

#' @rdname dot-point_est
#' @order 6
#' @param a.var Name of treatment variable.

.point_est.Y2pred <- function(
    data,
    cross.world,
    effect.scale,
    a.var,
    y.c1.form,
    y.c0.form,
    y.cm1.form,
    y.cm0.form,
    y10.c.form,
    y01.c.form,
    y.family,
    output.data = FALSE
) {

    data$.w.wt <- 1
    data$.f.wt <- data$.s.wt


    pred <- data[".f.wt"]

    y.c1.s11 <- glm(formula = y.c1.form,
                    data    = data[data$.a==1, ],
                    weights = data$.s.wt,
                    family  = y.family)

    pred$p11 <- predict(y.c1.s11, newdata = data, type = "response")


    y.c0.s00 <- glm(formula = y.c0.form,
                    data    = data[data$.a==0, ],
                    weights = data$.s.wt,
                    family  = y.family)

    pred$p00 <- predict(y.c0.s00, newdata = data, type = "response")




    y.var <- all.vars(formula(y.c1.form)[[2]])

    if ("10" %in% cross.world) {

        y.cm1.s11 <- glm(formula = y.cm1.form,
                         data    = data[data$.a==1, ],
                         weights = data$.s.wt,
                         family  = y.family)

        s00 <- data[data$.a==0, ]

        s00[, y.var] <- predict(y.cm1.s11, newdata = s00, type = "response")

        y10.c.s00 <- glm(formula = y10.c.form,
                         data    = s00,
                         weights = data$.s.wt,
                         family  = y.family)

        pred$p10 <- predict(y10.c.s00, newdata = data, type = "response")
    }



    if ("01" %in% cross.world) {

        y.cm0.s00 <- glm(formula = y.cm0.form,
                         data    = data[data$.a==0, ],
                         weights = data$.s.wt,
                         family  = y.family)

        s11 <- data[data$.a==1, ]

        s11[, y.var] <- predict(y.cm0.s00, newdata = s11, type = "response")

        y01.c.s11 <- glm(formula = y01.c.form,
                         data    = s11,
                         weights = data$.s.wt,
                         family  = y.family)

        pred$p01 <- predict(y01.c.s11, newdata = data, type = "response")
    }




    pred <- reshape_gather(pred,
                           columns = c("p00", "p10", "p01", "p11"),
                           key = ".samp",
                           value = ".y",
                           wide.row = FALSE)


    .get_means.and.effects(w.dat = pred,
                           effect.scale = effect.scale)
}






