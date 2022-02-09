



#### estimate_wpCadj ############################################################

#' Estimator wp-Cadj
#'
#' Function that implements estimator wp-Cadj
#' @inheritParams estimate_psYpred
#' @inheritParams estimate_wtd
#' @family estimators
#' @export

estimate_wpCadj <- function(
    data,
    s.wt.var     = NULL,
    cross.world  = "10",
    effect.scale = "MD",

    boot.num      = 999,
    boot.seed     = NULL,
    boot.method   = "cont-wt",
    boot.stratify = TRUE,

    a.c.form,
    max.stabilized.wt = 30,

    y.cm1.form = NULL,
    y.cm0.form = NULL,
    y.cm.form  = NULL,
    y.link,

    plot    = TRUE,
    c.order = NULL,
    c.std   = NULL
) {

    c.vars <- y.family <- wkng.form <- NULL

    .prep_wpCadj()


    key.inputs <- mget(c("cross.world",
                         "effect.scale",
                         "a.c.form",
                         "max.stabilized.wt",
                         "y.cm1.form", "y.cm0.form",
                         "wkng.form",
                         "y.family"))

    # POINT ESTIMATION

    if (!plot) {

        estimates <- do.call(".point_est.wpCadj", c(key.inputs,
                                                    list(data        = data,
                                                         output.data = FALSE)))
    } else {

        tmp <- do.call(".point_est.wpCadj", c(key.inputs,
                                              list(data        = data,
                                                   output.data = TRUE)))
        estimates <- tmp$estimates

        plots <- .plot_ipw(w.dat       = tmp$w.dat,
                           vars        = c.vars,
                           vars.std    = c.std,
                           key.balance = TRUE);     rm(tmp)
    }


    # BOOTSTRAP

    if (boot.num > 0) {
        ci.se <- .boot_ci.se(data       = data,
                             stratify   = boot.stratify,
                             boot.num   = boot.num,
                             seed       = boot.seed,
                             method     = boot.method,
                             FUN        = ".point_est.wpCadj",
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




#### .prep_wpCadj #######################################################

#' @rdname dot-prep
#' @order 14

.prep_wpCadj <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_weights.ipw(top.env)

    .clean_y.wpCadj(top.env)

    if (top.env$plot) .check_plot.ipw(top.env)
}




#### .clean_y.wpCadj ####################################################

#' @rdname dot-clean_y
#' @order 10

.clean_y.wpCadj <- function(env) {

    yes10 <- ("10" %in% env$cross.world)
    yes01 <- ("01" %in% env$cross.world)

    y.cm1  <- env$y.cm1.form
    y.cm0  <- env$y.cm0.form
    y.cm   <- env$y.cm.form
    y.link <- env$y.link



    if (is.null(y.cm)) {

        if (yes10 && is.null(y.cm1))
            stop("Must specify either y.cm1.form or y.cm.form.")

        if (yes01 && is.null(y.cm0))
            stop("Must specify either y.cm0.form or y.cm.form.")

    } else {

        if (yes10 && is.null(y.cm1)) env$y.cm1.form <- y.cm1 <- y.cm
        if (yes01 && is.null(y.cm0)) env$y.cm0.form <- y.cm0 <- y.cm
    }


    y.var <- NULL

    if (yes10) y.var <- unique(c(y.var, all.vars(formula(y.cm1)[[2]])))
    if (yes01) y.var <- unique(c(y.var, all.vars(formula(y.cm0)[[2]])))

    if (length(y.var)>1)
        stop("Outcome variable is not unique across outcome models.")

    if (!y.var %in% names(env$data))
        stop(paste("Variable", y.var, "(y.var) not found in dataset."))

    env$data$.y <- env$data[, y.var]




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



    env$wkng.form <- paste(y.var,
                           "~ .a +",
                           paste(env$c.vars, collapse = " + "))

}




#### .point_est.wpCadj ####################################################

#' @rdname dot-point_est
#' @order 10
#' @param wkng.form Formula for working regression model
#' @importFrom stats coef

.point_est.wtCadj <- function(
    data,
    cross.world,
    effect.scale,
    a.c.form,
    max.stabilized.wt = 30,
    output.data = FALSE,
    y.cm1.form,
    y.cm0.form,
    wkng.form,
    y.family
) {

    w.dat <- .compute_weights.ipw(
        data              = data,
        a.form            = a.c.form,
        max.stabilized.wt = max.stabilized.wt
    )

    p00 <- w.dat[w.dat$.samp=="p00", ]
    p11 <- w.dat[w.dat$.samp=="p11", ]




    estimates <- NULL



    if ("10" %in% cross.world) {

        # p10
        y.cm1.s11 <- glm(formula = y.cm1.form,
                         data    = data[data$.a==1, ],
                         weights = data$.s.wt,
                         family  = y.family)
        p10 <- p00
        p10$.samp <- "p10"
        p10$.y <- predict(y.cm1.s11, newdata = p10, type = "response")



        # estimate NIE1
        wkng.mod <- glm(formula = wkng.form,
                        data    = rbind(p11, p10),
                        weights = data$.f.wt,
                        family  = y.family)

        if (effect.scale=="MD" && y.family=="gaussian") {

            nie1 <- unname(coef(wkng.mod)[2])

        } else if (effect.scale=="MR" && y.family=="quasipoisson") {

            nie1 <- exp(unname(coef(wkng.mod)[2]))

        } else {

            pred <- w.dat

            pred$.a <- 1
            pred$y11 <- predict(wkng.mod, newdata = pred, type = "response")

            pred$.a <- 0
            pred$y10 <- predict(wkng.mod, newdata = pred, type = "response")

            nie1 <- .get_contrast(.wtd_mean(pred$y11, pred$.s.wt),
                                  .wtd_mean(pred$y10, pred$.s.wt),
                                  type = effect.scale)

            rm(pred)
        }

        rm(wkng.mod)


        # estimate NDE0
        if (effect.scale=="MD") {

            nde0 <- .wtd_mean(p10$.y - p00$.y, p00$.f.wt)

        } else {

            p10$.a <- 1
            wkng.mod <- glm(formula = wkng.form,
                            data    = rbind(p10, p00),
                            weights = data$.f.wt,
                            family  = y.family)

            if (effect.scale=="MR" && y.family=="quasipoisson") {

                nde0 <- exp(unname(coef(wkng.mod)[2]))

            } else {

                pred <- w.dat

                pred$.a <- 1
                pred$y10 <- predict(wkng.mod, newdata = pred, type = "response")

                pred$.a <- 0
                pred$y00 <- predict(wkng.mod, newdata = pred, type = "response")

                nde0 <- .get_contrast(.wtd_mean(pred$y10, pred$.s.wt),
                                      .wtd_mean(pred$y00, pred$.s.wt),
                                      type = effect.scale)

                rm(pred)
            }

            rm(wkng.mod)
        }

        estimates <- c(estimates,
                       NDE0 = nde0,
                       NIE1 = nie1)


    }




    if ("01" %in% cross.world) {

        # p01
        y.cm0.s00 <- glm(formula = y.cm0.form,
                         data    = data[data$.a==0, ],
                         weights = data$.s.wt,
                         family  = y.family)
        p01 <- p11
        p01$.samp <- "p01"
        p01$.y <- predict(y.cm0.s00, newdata = p01, type = "response")



        # estimate NIE0
        wkng.mod <- glm(formula = wkng.form,
                        data    = rbind(p01, p00),
                        weights = data$.f.wt,
                        family  = y.family)

        if (effect.scale=="MD" && y.family=="gaussian") {

            nie0 <- unname(coef(wkng.mod)[2])

        } else if (effect.scale=="MR" && y.family=="quasipoisson") {

            nie0 <- exp(unname(coef(wkng.mod)[2]))

        } else {

            pred <- w.dat

            pred$.a <- 1
            pred$y01 <- predict(wkng.mod, newdata = pred, type = "response")

            pred$.a <- 0
            pred$y00 <- predict(wkng.mod, newdata = pred, type = "response")

            nie0 <- .get_contrast(.wtd_mean(pred$y01, pred$.s.wt),
                                  .wtd_mean(pred$y00, pred$.s.wt),
                                  type = effect.scale)

            rm(pred)
        }

        rm(wkng.mod)


        # estimate NDE1
        if (effect.scale=="MD") {

            nde1 <- .wtd_mean(p11$.y - p01$.y, p11$.f.wt)

        } else {

            p01$.a <- 0
            wkng.mod <- glm(formula = wkng.form,
                            data    = rbind(p11, p01),
                            weights = data$.f.wt,
                            family  = y.family)

            if (effect.scale=="MR" && y.family=="quasipoisson") {

                nde1 <- exp(unname(coef(wkng.mod)[2]))

            } else {

                pred <- w.dat

                pred$.a <- 1
                pred$y11 <- predict(wkng.mod, newdata = pred, type = "response")

                pred$.a <- 0
                pred$y01 <- predict(wkng.mod, newdata = pred, type = "response")

                nde1 <- .get_contrast(.wtd_mean(pred$y11, pred$.s.wt),
                                      .wtd_mean(pred$y01, pred$.s.wt),
                                      type = effect.scale)

                rm(pred)
            }

            rm(wkng.mod)
        }

        estimates <- c(estimates,
                       NIE0 = nie0,
                       NDE1 = nde1)
    }



    if (!output.data) return(estimates)

    list(estimates = estimates,
         w.dat    = w.dat)


}


























