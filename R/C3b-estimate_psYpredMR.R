



#### OK  estimate_psYpredMR ####################################################

#' Estimator psYpredMR
#'
#' Function that implements estimator psYpredMR
#' @inheritParams estimate_psYpred
#' @inheritParams estimate_wtd
#' @family estimators
#' @family MR-estimators
#' @export

estimate_psYpredMR <- function(
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
    c.std = NULL,
    m.std = NULL,
    c.order = NULL,
    m.order = NULL,

    y.c.form = NULL,
    y.c1.form = NULL,
    y.c0.form = NULL,
    y.cm.form = NULL,
    y.cm1.form = NULL,
    y.cm0.form = NULL,
    y.link = "gaussian") {


    # CLEAN INPUTS

    c.vars <- m.vars <- y.family <- NULL

    .prep_psYpredMR()

    key.inputs <- mget(c("cross.world",
                         "effect.scale",
                         "a.c.form",
                         "a.cm.form",
                         "max.stabilized.wt",
                         "y.c1.form",
                         "y.c0.form",
                         "y.cm1.form",
                         "y.cm0.form",
                         "y.family"))



    # POINT ESTIMATION

    if (!plot) {

        estimates <- do.call(".point_est.psYpredMR",
                             c(key.inputs, list(data        = data,
                                                output.data = FALSE)))
    } else {

        tmp <- do.call(".point_est.psYpredMR",
                       c(key.inputs, list(data        = data,
                                          output.data = TRUE)))
        estimates <- tmp$estimates

        plots <- .plot_psYpredMR(w.dat = tmp$w.dat,
                                 cross.world = cross.world,
                                 c.vars = c.vars,
                                 m.vars = m.vars,
                                 c.std = c.std,
                                 m.std = m.std);  rm(tmp)
    }


    # BOOTSTRAP

    if (boot.num > 0) {
        ci.se <- .boot_ci.se(data       = data,
                             stratify   = boot.stratify,
                             boot.num   = boot.num,
                             seed       = boot.seed,
                             method     = boot.method,
                             FUN        = ".point_est.psYpredMR",
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






#### OK  .prep_psYpredMR #######################################################

#' @rdname dot-prep
#' @order 8

.prep_psYpredMR <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_effect.scale(top.env)

    .clean_boot(top.env)

    .clean_weights.med(top.env)

    .clean_y.psYpredMR(top.env)

    if (top.env$plot) .check_plot.med(top.env)
}



#### OK  .clean_y.psYpredMR ##############################################

#' @rdname dot-clean_y
#' @order 3
#'
.clean_y.psYpredMR <- function(env) {

    yes10 <- ("10" %in% env$cross.world)
    yes01 <- ("01" %in% env$cross.world)

    in_c.vars <- env$c.vars
    in_m.vars <- env$m.vars

    y.c    <- env$y.c.form
    y.cm   <- env$y.cm.form
    y.link <- env$y.link



    # populate formulas
    if (is.null(y.c)) {

        if (yes10 && is.null(env$y.c1.form))
            stop("Must specify either y.c1.form or y.c.form.")

        if (yes01 && is.null(env$y.c0.form))
            stop("Must specify either y.c0.form or y.c.form.")

    } else {

        if (yes10 && is.null(env$y.c1.form))
            env$y.c1.form <- y.c

        if (yes01 && is.null(env$y.c0.form))
            env$y.c0.form <- y.c
    }


    if (is.null(y.cm)) {

        if (yes10 && is.null(env$y.cm1.form))
            stop("Must specify either y.cm1.form or y.cm.form.")

        if (yes01 && is.null(env$y.cm0.form))
            stop("Must specify either y.cm0.form or y.cm.form.")

    } else {

        if (yes10 && is.null(env$y.cm1.form))
            env$y.cm1.form <- y.cm

        if (yes01 && is.null(env$y.cm0.form))
            env$y.cm0.form <- y.cm
    }



    # check y.var, c.vars, m.vars
    y.var  <- NULL
    c.vars <- NULL
    m.vars <- NULL

    if (yes10) {
        y.var <- unique(c(y.var,
                          all.vars(formula(env$y.c1.form)[[2]]),
                          all.vars(formula(env$y.cm1.form)[[2]])))
        c.vars <- unique(c(c.vars,
                           all.vars(formula(env$y.c1.form)[[3]])))
    }


    if (yes01) {
        y.var <- unique(c(y.var,
                          all.vars(formula(env$y.c0.form)[[2]]),
                          all.vars(formula(env$y.cm0.form)[[2]])))
        c.vars <- unique(c(c.vars,
                           all.vars(formula(env$y.c0.form)[[3]])))
    }


    if (length(y.var)>1)
        stop("Outcome variable is not unique across outcome models.")

    env$data$.y <- env$data[, y.var]

    stray.c <- setdiff(c.vars, in_c.vars)

    if (length(stray.c)>0)
        stop(paste("Covariate(s)",
                   paste(stray.c, collapse = ", "),
                   "(that appear in outcome given covariates model(s)) are not found in a.c.form."))


    if (yes10)
        m.vars <- unique(c(m.vars,
                           setdiff(all.vars(formula(env$y.cm1.form)[[3]]),
                                   in_c.vars)))
    if (yes01)
        m.vars <- unique(c(m.vars,
                           setdiff(all.vars(formula(env$y.cm0.form)[[3]]),
                                   in_c.vars)))

    stray.m <- setdiff(m.vars, in_m.vars)

    if (length(stray.m)>0)
        stop(paste("Mediator(s)",
                   paste(stray.m, collapse = ", "),
                   "(that appear in outcome model(s)) are not part of the mediators based on a.c.form and a.cm.form."))




    # y.link -> y.family
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











#### OK  .point_est.psYpredMR ##################################################

#' @rdname dot-point_est
#' @order 3

.point_est.psYpredMR <- function(
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
    y.family
) {

    dat <- .compute_weights.med(
        data              = data,
        cross.world       = cross.world,
        a.c.form          = a.c.form,
        a.cm.form         = a.cm.form,
        max.stabilized.wt = max.stabilized.wt
    )

    rm(data)

    estimates <- NULL

    if ("10" %in% cross.world) {

        p00 <- dat[dat$.samp=="p00", ]

        y.c1.p11 <- glm(formula = y.c1.form,
                        data = dat[dat$.samp=="p11", ],
                        weights = data$.f.wt,
                        family = y.family)
        y.cm1.p10 <- glm(formula = y.cm1.form,
                         data = dat[dat$.samp=="p10", ],
                         weights = data$.f.wt,
                         family = y.family)

        pred10 <- p00[".f.wt"]
        pred10$p00 <- p00$.y
        pred10$p10 <- predict(y.cm1.p10, newdata = p00, type = "response")
        pred10$p11 <- predict(y.c1.p11,  newdata = p00, type = "response")

        pred10 <- reshape_gather(pred10, columns = c("p00", "p10", "p11"),
                                 key = ".samp",
                                 value = ".y",
                                 wide.row = FALSE)

        estimates <- c(estimates,
                       .get_means.and.effects(w.dat = pred10,
                                              effect.scale = effect.scale))
    }

    if ("01" %in% cross.world) {

        p11 <- dat[dat$.samp=="p11", ]

        y.c0.p00 <- glm(formula = y.c0.form,
                        data = dat[dat$.samp=="p00", ],
                        weights = data$.f.wt,
                        family = y.family)
        y.cm0.p01 <- glm(formula = y.cm0.form,
                         data = dat[dat$.samp=="p01", ],
                         weights = data$.f.wt,
                         family = y.family)

        pred01 <- p11[".f.wt"]
        pred01$p11 <- p11$.y
        pred01$p01 <- predict(y.cm0.p01, newdata = p11, type = "response")
        pred01$p00 <- predict(y.c0.p00,  newdata = p11, type = "response")

        pred01 <- reshape_gather(pred01, columns = c("p00", "p01", "p11"),
                                 key = ".samp",
                                 value = ".y",
                                 wide.row = FALSE)

        estimates <- c(estimates,
                       .get_means.and.effects(w.dat = pred01,
                                              effect.scale = effect.scale))
    }

    if (!output.data) return(estimates)

    list(estimates = estimates,
         w.dat    = dat)
}





#### OK  .plot_psYpredMR #######################################################

#' @rdname dot-plot_w.dat
#' @order 5
#' @param cross.world blah

.plot_psYpredMR <- function(w.dat,
                            cross.world,
                            c.vars,
                            m.vars,
                            c.std,
                            m.std) {

    c(.plot_wt_dist(w.dat),
      .plot_balance.psYpredMR(w.dat = w.dat,
                              cross.world = cross.world,
                              c.vars = c.vars,
                              m.vars = m.vars,
                              c.std = c.std,
                              m.std = m.std))
}



#### OK  .plot_balance.psYpredMR ###############################################

#' @rdname dot-plot_balance
#' @order 5
#' @param cross.world (For \code{plot_balance.psYpredMR}) blah

.plot_balance.psYpredMR <- function(w.dat,
                                    cross.world,
                                    c.vars,
                                    m.vars,
                                    c.std,
                                    m.std) {



    smd.dat <- .get_smd.med(w.dat = w.dat,
                            c.vars = c.vars,
                            m.vars = m.vars,
                            c.std = c.std,
                            m.std = m.std)



    full.balance <-
        ggplot(data = smd.dat,
               aes(x = .data$mean.diff,
                   y = factor(.data$variable,
                              levels = rev(levels(.data$variable))))) +
        geom_vline(xintercept = 0,
                   color = "gray60") +
        geom_point(aes(color = .data$var.type,
                       shape = .data$contrast.type),
                   fill = "white",
                   size = 1.5,
                   stroke = .5) +
        labs(x = "differences in means",
             y = "") +
        scale_color_manual(name = "", values = c("black", "magenta")) +
        scale_shape_manual(name = "", values = c(21, 19)) +
        theme_bw() +
        xlim(min(c(-.3, smd.dat$mean.diff)),
             max(c( .3, smd.dat$mean.diff))) +
        facet_wrap(~.data$contrast, ncol = 3)



    if (all(c("10", "01") %in% cross.world)) {
        smd.dat <- smd.dat[smd.dat$contrast %in% c("p00 - full",
                                                   "p11 - full"), ]
    } else if (cross.world=="10") {
        smd.dat <- smd.dat[smd.dat$contrast=="p00 - full", ]
    } else if (cross.world=="01") {
        smd.dat <- smd.dat[smd.dat$contrast=="p11 - full", ]
    }

    key.balance <-
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
        facet_wrap(~.data$contrast, ncol = 2)


    mget(c("key.balance", "full.balance"))
}



