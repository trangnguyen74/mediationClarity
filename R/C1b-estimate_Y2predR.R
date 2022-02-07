



#### OK  estimate_Y2predR #####################################################

#' Estimator Y2pred.R
#'
#' Function that implements estimator Y2pred.MR
#' @inheritParams estimate_Y2pred
#' @inheritParams estimate_wtd
#' @family estimators
#' @family more-robust estimators
#' @export

estimate_Y2predR <- function(
    data,
    s.wt.var     = NULL,
    cross.world  = "10",
    effect.scale = "additive",

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

    # detailed outcome models
    y.c1.form = NULL,
    y.c0.form = NULL,

    y.cm1.form = NULL,
    y.cm0.form = NULL,

    y10.c.form = NULL,
    y01.c.form = NULL,

    # or shotcuts
    y.c.form = NULL,
    y.cm.form = NULL,

    y.link = "identity") {


    # CLEAN INPUTS

    c.vars <- m.vars <- y.family <- NULL

    .prep_Y2predR()

    key.inputs <- mget(c("cross.world",
                         "effect.scale",
                         "a.c.form",
                         "a.cm.form",
                         "max.stabilized.wt",
                         "y.c1.form",
                         "y.c0.form",
                         "y.cm1.form",
                         "y.cm0.form",
                         "y10.c.form",
                         "y01.c.form",
                         "y.family"))


    # POINT ESTIMATION

    if (!plot) {

        estimates <- do.call(".point_est.Y2predR",
                             c(key.inputs, list(data        = data,
                                                output.data = FALSE)))
    } else {

        tmp <- do.call(".point_est.Y2predR",
                       c(key.inputs, list(data        = data,
                                          output.data = TRUE)))

        estimates <- tmp$estimates

        plots <- .plot_med(w.dat = tmp$w.dat,
                           c.vars = c.vars,
                           m.vars = m.vars,
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
                             FUN        = ".point_est.Y2predR",
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




#### OK  .prep_Y2predR #########################################################

#' @rdname dot-prep
#' @order 9

.prep_Y2predR <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_effect.scale(top.env)

    .clean_boot(top.env)

    .clean_weights.med(top.env)

    .clean_y.forms.Y2predR(top.env)

    if (top.env$plot) .check_plot.med(top.env)
}




#### OK  .clean_y.forms.Y2preR #################################################

#' @rdname dot-clean_y
#' @order 7

.clean_y.forms.Y2predR <- function(env) {

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

    in_c.vars <- env$c.vars
    in_m.vars <- env$m.vars



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



    if (!all(c.vars %in% in_c.vars))
        stop("Not all covariates in outcome models are contained in a.c.form.")


    m.vars <- setdiff(cm.vars, in_c.vars)

    if (length(m.vars)==0)
        stop("The combination of outcome models imply all predictors are covariates and there are no mediators. Please double-check!")


    if (!all(m.vars %in% in_m.vars))
        stop("Not all mediators in outcome models are contained in a.cm.form.")






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




#### OK  .point_est.Y2predR ####################################################

#' @rdname dot-point_est
#' @order 7

.point_est.Y2predR <- function(
    data,
    cross.world,
    effect.scale,
    a.c.form,
    a.cm.form,
    max.stabilized.wt = 30,
    output.data = FALSE,
    y.c1.form,
    y.c0.form,
    y.cm1.form,
    y.cm0.form,
    y10.c.form,
    y01.c.form,
    y.family
) {


    w.dat <- .compute_weights.med(
        data              = data,
        cross.world       = cross.world,
        a.c.form          = a.c.form,
        a.cm.form         = a.cm.form,
        max.stabilized.wt = max.stabilized.wt
    )



    p11 <- w.dat[w.dat$.samp=="p11", ]
    p00 <- w.dat[w.dat$.samp=="p00", ]

    full <- rbind(p11, p00)
    full$.w.wt <- 1
    full$.f.wt <- full$.s.wt

    y.c1.p11 <- glm(formula = y.c1.form,
                    data = p11,
                    weights = data$.f.wt,
                    family = y.family)

    y.c0.p00 <- glm(formula = y.c0.form,
                    data = p00,
                    weights = data$.f.wt,
                    family = y.family)

    pred <- full[".f.wt"]
    pred$p00 <- predict(y.c0.p00, newdata = full, type = "response")
    pred$p11 <- predict(y.c1.p11, newdata = full, type = "response")




    y.var <- all.vars(formula(y.c1.form)[[2]])

    if ("10" %in% cross.world) {

        y.cm1.p10 <- glm(formula = y.cm1.form,
                         data    = w.dat[w.dat$.samp=="p10", ],
                         weights = data$.f.wt,
                         family  = y.family)

        p00.tmp <- p00

        p00.tmp[, y.var] <- predict(y.cm1.p10, newdata = p00, type = "response")

        y10.c.p00 <- glm(formula = y10.c.form,
                         data    = p00.tmp,
                         weights = data$.f.wt,
                         family  = y.family)

        pred$p10 <- predict(y10.c.p00, newdata = full, type = "response")
    }



    if ("01" %in% cross.world) {

        y.cm0.p01 <- glm(formula = y.cm0.form,
                         data    = w.dat[w.dat$.samp=="p01", ],
                         weights = data$.s.wt,
                         family  = y.family)

        p11.tmp <- p11

        p11.tmp[, y.var] <- predict(y.cm0.p01, newdata = p11, type = "response")

        y01.c.p11 <- glm(formula = y01.c.form,
                         data    = p11.tmp,
                         weights = data$.s.wt,
                         family  = y.family)

        pred$p01 <- predict(y01.c.p11, newdata = full, type = "response")
    }




    pred <- reshape_gather(pred,
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






