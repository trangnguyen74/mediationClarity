



#### OK  estimate_psYpred ######################################################

#' Estimator psYpred
#'
#' Function that implements estimator psYpred
#' @inheritParams estimate_wtd
#' @param y.c.form,y.c1.form,y.c0.form blah
#' @param y.cm.form,y.cm1.form,y.cm0.form blah
#' @param y.link blah
#' @family estimators
#' @export

estimate_psYpred <- function(
    data,
    s.wt.var     = NULL,
    cross.world  = "10",
    effect.scale = "additive",

    boot.num      = 999,
    boot.seed     = NULL,
    boot.method   = "cont-wt",
    boot.stratify = TRUE,

    a.c.form,
    max.stabilized.wt = 30,

    plot = TRUE,
    c.vars.std = NULL,
    c.order = NULL,

    y.c.form = NULL,
    y.c1.form = NULL,
    y.c0.form = NULL,
    y.cm.form = NULL,
    y.cm1.form = NULL,
    y.cm0.form = NULL,
    y.link = "gaussian") {


    # CLEAN INPUTS

    c.vars <- y.family <- NULL

    .prep_psYpred()

    key.inputs <- mget(c("cross.world",
                         "effect.scale",
                         "a.c.form",
                         "max.stabilized.wt",
                         "y.c1.form",
                         "y.c0.form",
                         "y.cm1.form",
                         "y.cm0.form",
                         "y.family"))



    # POINT ESTIMATION

    if (!plot) {

        estimates <- do.call(".point_est.psYpred", c(key.inputs,
                                                     list(data        = data,
                                                          output.data = FALSE)))
    } else {

        tmp <- do.call(".point_est.psYpred", c(key.inputs,
                                               list(data        = data,
                                                    output.data = TRUE)))
        estimates <- tmp$estimates

        plots <- .plot_psYpred(w.dat = tmp$w.dat,
                               c.vars = c.vars,
                               c.vars.std = c.vars.std);  rm(tmp)
    }


    # BOOTSTRAP

    if (boot.num > 0) {
        ci.se <- .boot_ci.se(data       = data,
                             stratify   = boot.stratify,
                             boot.num   = boot.num,
                             seed       = boot.seed,
                             method     = boot.method,
                             FUN        = ".point_est.psYpred",
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






#### OK  .prep_psYpred #########################################################

#' @rdname dot-prep
#' @order 5

.prep_psYpred <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_effect.scale(top.env)

    .clean_boot(top.env)

    .clean_weights.ipw(top.env)

    .clean_y.forms.psYpred(top.env)

    if (top.env$plot) .check_plot.ipw(top.env)
}



#### OK  .clean_y.forms.psYpred ################################################

#' @rdname dot-clean_y
#' @order 2

.clean_y.forms.psYpred <- function(env) {

    yes10 <- ("10" %in% env$cross.world)
    yes01 <- ("01" %in% env$cross.world)

    in_c.vars <- env$c.vars

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



    # check y.var and c.vars
    y.var <- NULL
    c.vars <- NULL

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




#### OK  .point_est.psYpred ####################################################

#' @rdname dot-point_est
#' @order 2
#' @param y.c1.form,y.c0.form Model formulas for E[Y|C,A=1], E[Y|C,A=0], checked.
#' @param y.cm1.form,y.cm0.form Model formulas for E[Y|C,M,A=1], E[Y|C,M,A=0], checked.
#' @param y.family GLM (quasi-) family for outcome models.

.point_est.psYpred <- function(
    data,
    cross.world,
    effect.scale,
    a.c.form,
    max.stabilized.wt = 30,
    output.data = FALSE,
    y.c1.form,
    y.c0.form,
    y.cm1.form,
    y.cm0.form,
    y.family
) {

    dat <- .compute_weights.psYpred(
        data              = data,
        cross.world       = cross.world,
        a.c.form          = a.c.form,
        max.stabilized.wt = max.stabilized.wt
    )

    rm(data)

    estimates <- NULL

    if ("10" %in% cross.world) {

        p00 <- dat[dat$.samp=="p00", ]

        y.c1.s11 <- glm(formula = y.c1.form,
                        data = dat[dat$.samp=="s11", ],
                        weights = data$.s.wt,
                        family = y.family)
        y.cm1.s11 <- glm(formula = y.cm1.form,
                         data = dat[dat$.samp=="s11", ],
                         weights = data$.s.wt,
                         family = y.family)

        pred10 <- p00[".f.wt"]
        pred10$p00 <- p00$.y
        pred10$p10 <- predict(y.cm1.s11, newdata = p00, type = "response")
        pred10$p11 <- predict(y.c1.s11,  newdata = p00, type = "response")

        pred10 <- reshape_gather(pred10,
                                 columns = c("p00", "p10", "p11"),
                                 key = ".samp",
                                 value = ".y",
                                 wide.row = FALSE)

        estimates <- c(estimates,
                       .get_means.and.effects(w.dat = pred10,
                                              effect.scale = effect.scale))
    }

    if ("01" %in% cross.world) {

        p11 <- dat[dat$.samp=="p11", ]

        y.c0.s00 <- glm(formula = y.c0.form,
                        data = dat[dat$.samp=="s00", ],
                        weights = data$.s.wt,
                        family = y.family)
        y.cm0.s00 <- glm(formula = y.cm0.form,
                         data = dat[dat$.samp=="s00", ],
                         weights = data$.s.wt,
                         family = y.family)

        pred01 <- p11[".f.wt"]
        pred01$p11 <- p11$.y
        pred01$p01 <- predict(y.cm0.s00, newdata = p11, type = "response")
        pred01$p00 <- predict(y.c0.s00,  newdata = p11, type = "response")

        pred01 <- reshape_gather(pred01,
                                 columns = c("p00", "p01", "p11"),
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





#### OK  .compute_weights_psYpred ##############################################

#' @rdname dot-compute_weights
#' @order 4

.compute_weights.psYpred <- function(
    data,
    cross.world,
    a.c.form,
    max.stabilized.wt
) {

    tmp <- .compute_weights.ipw(data = data,
                                a.c.form = a.c.form,
                                max.stabilized.wt = max.stabilized.wt)

    out <- NULL

    if ("10" %in% cross.world) {
        p00 <- tmp[tmp$.a==0, ]

        s11 <- tmp[tmp$.a==1, ]
        s11$.samp <- "s11"
        s11$.w.wt <- 1
        s11$.f.wt <- s11$.s.wt

        out <- rbind(out, p00, s11)
    }

    if ("01" %in% cross.world) {
        p11 <- tmp[tmp$.a==1, ]

        s00 <- tmp[tmp$.a==0, ]
        s00$.samp <- "s00"
        s00$.w.wt <- 1
        s00$.f.wt <- s00$.s.wt

        out <- rbind(out, p11, s00)
    }

    out
}




#### OK  .plot_psYpred #########################################################

#' @rdname dot-plot_w.dat
#' @order 4

.plot_psYpred <- function(w.dat,
                          c.vars,
                          c.vars.std) {

    c(.plot_wt_dist.psYpred(w.dat),
      .plot_balance.psYpred(w.dat      = w.dat,
                            c.vars     = c.vars,
                            c.vars.std = c.vars.std))
}



#### OK  .plot_wt_dist.psYpred #################################################

#' @rdname dot-plot_wt_dist
#' @order 3

.plot_wt_dist.psYpred <- function(
    w.dat,
    point.alpha = .1,
    jitter.width = .3
) {


    w.dat$.w.wt <- stabilize_weight(weight   = w.dat$.w.wt,
                                    group    = w.dat$.samp,
                                    s.weight = w.dat$.s.wt)


    if (is_constant(w.dat$.s.wt)) {

        dat <- w.dat[, c(".samp", ".w.wt")]
        dat <- dat[dat$.samp %in% c("p00", "p11"), ]


        p <- ggplot(data = dat,
                    aes(x = .data$.samp,
                        y = .data$.w.wt)) +
            geom_jitter(height = 0,
                        width = jitter.width,
                        alpha = point.alpha) +
            geom_violin(color = "red", fill = NA) +
            labs(x = "",
                 y = "distribution morphing weights (stabilized)") +
            theme_bw()

        return(list(w.wt.distribution = p))
    }




    w.dat$.s.wt <- stabilize_weight(weight   = w.dat$.s.wt,
                                    group    = w.dat$.samp,
                                    s.weight = rep(1, nrow(w.dat)))

    dat <- w.dat[, c(".samp", ".s.wt", ".w.wt")]
    dat <- dat[dat$.samp %in% c("p00", "p11"), ]

    w.wt.distribution <-
        ggplot(data = dat,
               aes(x = .data$.samp,
                   y = .data$.w.wt,
                   weight = .data$.s.wt)) +
        geom_jitter(height = 0,
                    width = jitter.width,
                    alpha = point.alpha,
                    aes(size = .data$.s.wt)) +
        geom_violin(color = "red", fill = NA) +
        labs(x = "",
             y = "distribution morphing weights (stabilized)",
             size = "sampling weight (stabilized)") +
        theme_bw() +
        theme(legend.position = "bottom")




    w.dat$.f.wt <- stabilize_weight(weight   = w.dat$.f.wt,
                                    group    = w.dat$.samp,
                                    s.weight = rep(1, nrow(w.dat)))

    dat <- w.dat[, c(".samp", ".f.wt")]

    f.wt.distribution <-
        ggplot(data = dat,
               aes(x = .data$.samp,
                   y = .data$.f.wt)) +
        geom_jitter(height = 0,
                    width = jitter.width,
                    alpha = point.alpha) +
        geom_violin(color = "red", fill = NA) +
        labs(x = "",
             y = "final weights (stabilized)\n(combining sampling and distribution morphing)") +
        theme_bw()


    mget(c("w.wt.distribution", "f.wt.distribution"))

}




#### OK  .plot_balance.psYpred #################################################

#' @rdname dot-plot_balance
#' @order 4


.plot_balance.psYpred <- function(w.dat,
                                  c.vars,
                                  c.vars.std) {

    smd.dat <- .get_smd.psYpred(w.dat = w.dat,
                                c.vars = c.vars,
                                c.vars.std = c.vars.std)


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
        facet_wrap(~.data$contrast, ncol = 2)

    list(key.balance = p)

}




#### OK  .get_smd.psYpred ######################################################

#' @rdname dot-get_smd
#' @order 4

.get_smd.psYpred <- function(w.dat,
                             c.vars,
                             c.vars.std) {


    tmp <- .make_dummies(data = w.dat,
                         columns = c.vars,
                         output.names = TRUE,
                         warning = FALSE)

    w.dat <- tmp$data
    c.vars <- tmp$columns; rm(tmp)


    yes10 <- any(w.dat$.samp=="p00")
    yes01 <- any(w.dat$.samp=="p11")

    if (yes10) { full <- w.dat[w.dat$.samp %in% c("p00", "s11"), ]
    } else     { full <- w.dat[w.dat$.samp %in% c("s00", "p11"), ]

    }
    full$.w.wt <- 1
    full$.f.wt <- full$.s.wt

    if (yes10) p00 <- w.dat[w.dat$.samp=="p00", ]
    if (yes01) p11 <- w.dat[w.dat$.samp=="p11", ]


    std.denom <- sapply(c.vars, function(z) {

        if (!z %in% c.vars.std) return(1)

        .get_sd.pooled(variable = z,
                       dat1 = full[full$.a==1, ],
                       dat0 = full[full$.a==0, ])
    })


    # compute (standardized) mean differences
    smd <- lapply(list(unw = ".s.wt", wtd = ".f.wt"), function(w) {

        means.full <- sapply(c.vars,
                             function(z) .wtd_mean(full[, z], full[, w]))

        diff <- NULL

        if (yes10) {

            means.p00 <- sapply(c.vars,
                                function(z) .wtd_mean(p00[, z], p00[, w]))

            diff <- cbind(diff,
                          p00.full = (means.p00 - means.full) / std.denom)

            rm(means.p00)
        }

        if (yes01) {

            means.p11 <- sapply(c.vars,
                                function(z) .wtd_mean(p11[, z], p11[, w]))

            diff <- cbind(diff,
                          p11.full = (means.p11 - means.full) / std.denom)

            rm(means.p11)
        }

        diff.names <- colnames(diff)

        diff <- data.frame(diff, row.names = NULL)
        diff$variable <- c.vars

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




    smd$contrast <- factor(smd$contrast,
                           levels = c("p00.full", "p11.full"),
                           labels = c("p00 - full", "p11 - full"))


    smd$variable <- ifelse(smd$variable %in% c.vars.std,
                           paste0("*", smd$variable),
                           smd$variable)

    c.vars <- ifelse(c.vars %in% c.vars.std,
                     paste0("*", c.vars),
                     c.vars)

    smd$variable <- factor(smd$variable, levels = c.vars)


    smd[,c("variable", "contrast.type", "contrast", "mean.diff")]
}







































#### .clean_y.forms_c.vars ### MAY RETIRE ######################################

.clean_y.forms_c.vars <- function() {

    inputs <- mget(c("data",
                     "y.c.form", "y.cm.form",
                     "c.vars"),
                   envir = parent.frame())


    if (!all(all.vars(formula(inputs$y.c.form)[[3]]) %in% inputs$c.vars))
        stop("y.c.form contains covariates that are not in a.c.form.")


    y.var <- all.vars(formula(inputs$y.c.form)[[2]])

    if (y.var != all.vars(formula(inputs$y.cm.form)[[2]]))
        stop("Outcome variable is not the same in y.c.form and y.cm.form.")

    inputs$data$.y <- inputs$data[, y.var]


    inputs$m.vars <- setdiff(all.vars(formula(inputs$y.c.form)[[3]]),
                             inputs$c.vars)


    inputs
}





















