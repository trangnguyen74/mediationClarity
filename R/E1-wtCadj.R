



#### estimate_wtCadj ############################################################

#' Estimator wt-Cadj
#'
#' Function that implements estimator wt-Cadj
#' @inheritParams estimate_psYpred
#' @inheritParams estimate_wtd
#' @family estimators
#' @export

estimate_wtCadj <- function(
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
    max.stabilized.wt = 30,

    y.var,
    y.link,

    plot       = TRUE,
    c.order    = NULL,
    m.order    = NULL,
    c.std = NULL,
    m.std = NULL
) {

    c.vars <- m.vars <- y.family <- wkng.form <- NULL

    .prep_wtCadj()


    key.inputs <- mget(c("cross.world",
                         "effect.scale",
                         "a.c.form",
                         "a.cm.form",
                         "max.stabilized.wt",
                         "wkng.form",
                         "y.family"))



    # POINT ESTIMATION

    if (!plot) {

        estimates <- do.call(".point_est.wtCadj", c(key.inputs,
                                                    list(data        = data,
                                                         output.data = FALSE)))
    } else {

        tmp <- do.call(".point_est.wtCadj", c(key.inputs,
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
                             FUN        = ".point_est.wtCadj",
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




#### .pred_wtCadj #######################################################

#' @rdname dot-prep
#' @order 13

.prep_wtCadj <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_weights.med(top.env)

    .clean_y.wtCadj(top.env)

    if (top.env$plot) .check_plot.med(top.env)
}




#### .clean_y.wtCadj ####################################################

#' @rdname dot-clean_y
#' @order 9

.clean_y.wtCadj <- function(env) {

    y.var <- env$y.var

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
                           "~ .a0 + .a1 +",
                           paste(env$c.vars, collapse = " + "))
}




#### .point_est.wtCadj ####################################################

#' @rdname dot-point_est
#' @order 10
#' @param wkng.form Formula for working regression model
#' @importFrom stats coef

.point_est.wtCadj <- function(
    data,
    cross.world,
    effect.scale,
    a.c.form,
    a.cm.form,
    max.stabilized.wt = 30,
    output.data = FALSE,
    wkng.form,
    y.family
) {

    w.dat <- .compute_weights.med(
        data              = data,
        cross.world       = cross.world,
        a.c.form          = a.c.form,
        a.cm.form         = a.cm.form,
        max.stabilized.wt = max.stabilized.wt
    )

    w.dat$.a0 <- ifelse(w.dat$.samp=="p00", -1, 0)
    w.dat$.a1 <- ifelse(w.dat$.samp=="p11",  1, 0)


    estimates <- NULL

    if ("10" %in% cross.world) {

        wkng.mod <- glm(formula = wkng.form,
                        data    = w.dat[w.dat$.samp!="p01", ],
                        weights = data$.f.wt,
                        family  = y.family)


        if (effect.scale=="MD" && y.family=="gaussian") {

            estimates <- c(estimates,
                           NDE0 = unname(coef(wkng.mod)[2]),
                           NIE1 = unname(coef(wkng.mod)[3]))

        } else if (effect.scale=="MR" && y.family=="quasipoisson") {

            estimates <- c(estimates,
                           exp(c(NDE0 = unname(coef(wkng.mod)[2]),
                                 NIE1 = unname(coef(wkng.mod)[3]))))

        } else {
            pred <- w.dat[w.dat$.samp %in% c("p00", "p11"), ]

            pred$.a0 <- -1;  pred$.a1 <- 0
            pred$y00 <- predict(wkng.mod, newdata = pred, type = "response")

            pred$.a0 <-  0;  pred$.a1 <- 0
            pred$y10 <- predict(wkng.mod, newdata = pred, type = "response")

            pred$.a0 <-  0;  pred$.a1 <- 1
            pred$y11 <- predict(wkng.mod, newdata = pred, type = "response")

            estimates <-
                c(estimates,
                  NDE0 = .get_contrast(.wtd_mean(pred$y10, pred$.s.wt),
                                       .wtd_mean(pred$y00, pred$.s.wt),
                                       type = effect.scale),
                  NIE1 = .get_contrast(.wtd_mean(pred$y11, pred$.s.wt),
                                       .wtd_mean(pred$y10, pred$.s.wt),
                                       type = effect.scale))

            rm(wkng.mod, pred)
        }
    }

    if ("01" %in% cross.world) {

        wkng.mod <- glm(formula = wkng.form,
                        data    = w.dat[w.dat$.samp!="p10"],
                        weights = data$.f.wt,
                        family  = y.family)


        if (effect.scale=="MD" && y.family=="gaussian") {

            estimates <- c(estimates,
                           NIE0 = unname(coef(wkng.mod)[2]),
                           NDE1 = unname(coef(wkng.mod)[3]))

        } else if (effect.scale=="MR" && y.family=="quasipoisson") {

            estimates <- c(estimates,
                           exp(c(NIE0 = unname(coef(wkng.mod)[2]),
                                 NDE1 = unname(coef(wkng.mod)[3]))))

        } else {

            pred <- w.dat[w.dat$.samp %in% c("p00", "p11"), ]

            pred$.a0 <- -1;  pred$.a1 <- 0
            pred$y00 <- predict(wkng.mod, newdata = pred, type = "response")

            pred$.a0 <-  0;  pred$.a1 <- 0
            pred$y01 <- predict(wkng.mod, newdata = pred, type = "response")

            pred$.a0 <-  0;  pred$.a1 <- 1
            pred$y11 <- predict(wkng.mod, newdata = pred, type = "response")

            estimates <-
                c(estimates,
                  NIE0 = .get_contrast(.wtd_mean(pred$y01, pred$.s.wt),
                                       .wtd_mean(pred$y00, pred$.s.wt),
                                       type = effect.scale),
                  NDE1 = .get_contrast(.wtd_mean(pred$y11, pred$.s.wt),
                                       .wtd_mean(pred$y01, pred$.s.wt),
                                       type = effect.scale))
        }
    }

    w.dat <- w.dat[, !names(w.dat) %in% c(".a0", ".a1")]


    if (!output.data) return(estimates)

    list(estimates = estimates,
         w.dat    = w.dat)
}




#### .get_means.contrasts ## MAYBE NO NEED? #######

.get_means.contrasts <- function(dat,
                                 effect.scale,
                                 same.sample) {

    if (same.sample) {

        po.names <- setdiff(names(dat), "wt")

        if (!all(substr(po.names, 1, 1)=="y"))
            stop("Predicted outcome variable names (",
                 paste(po.names, collapse = ", "),
                 ") do not all start with \"y\".")

        po.means <- lapply(po.names, function(z)
            .wtd_mean(dat[, z], dat$wt))

        names(po.means) <- paste0(po.names,
                                  ".mean")

    } else {

        ps.names <- unique(dat$.samp)

        po.means <- lapply(ps.names, function(z) {
            dat <- dat[dat$.samp==z, ]
            .wtd_mean(dat$.y, dat$.f.wt)
        })

        names(po.means) <- paste0("y",
                                  substr(ps.names, 2, nchar(ps.names)),
                                  ".mean")

    }


    effects <- c(TE = .get_contrast(po.means$y11.mean, po.means$y00.mean,
                                    type = effect.scale))

    if ("y10" %in% names(po.means)) {
        effects <- c(effects,
                     NDE0 = .get_contrast(po.means$y10.mean, po.means$y00.mean,
                                          type = effect.scale))
        effects <- c(effects,
                     NIE1 = .get_contrast(po.means$y11.mean, po.means$y10.mean,
                                          type = effect.scale))
    }

    if ("y01" %in% names(po.means)) {
        effects <- c(effects,
                     NIE0 = .get_contrast(po.means$y01.mean, po.means$y00.mean,
                                          type = effect.scale))
        effects <- c(effects,
                     NDE1 = .get_contrast(po.means$y11.mean, po.means$y01.mean,
                                          type = effect.scale))
    }

    c(unlist(po.means), effects)
}


