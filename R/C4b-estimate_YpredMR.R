



#### OK  estimate_YpredMR ######################################################

#' Estimator Ypred.MR
#'
#' Function that implements estimator Ypred.MR
#' @inheritParams estimate_wtd
#' @inheritParams estimate_Ypred
#' @family estimators
#' @family MR-estimators
#' @export

estimate_YpredMR <- function(
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

    plot = TRUE,
    c.std = NULL,
    m.std = NULL,
    c.order    = NULL,
    m.order    = NULL,

    y.c.form   = NULL,
    y.c1.form  = NULL,
    y.c0.form  = NULL,
    y10.c.form = NULL,
    y01.c.form = NULL,
    y.link     = "identity"
) {


    # CLEAN INPUTS

    c.vars <- m.vars <- y.family <- NULL

    .prep_YpredMR()


    key.inputs <- mget(c("cross.world",
                         "effect.scale",
                         "a.c.form",
                         "a.cm.form",
                         "max.stabilized.wt",
                         "y.c1.form",
                         "y.c0.form",
                         "y10.c.form",
                         "y01.c.form",
                         "y.family"))




    # POINT ESTIMATION

    if (!plot) {

        estimates <- do.call(".point_est.YpredMR",
                             c(key.inputs, list(data        = data,
                                                output.data = FALSE)))
    } else {

        tmp <- do.call(".point_est.YpredMR",
                       c(key.inputs, list(data        = data,
                                          output.data = TRUE)))

        estimates <- tmp$estimates

        plots <- .plot_YpredMR(w.dat = tmp$w.dat,
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
                             FUN        = ".point_est.YpredMR",
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




#### OK  .prep_YpredMR #########################################################

#' @rdname dot-prep
#' @order 10

.prep_YpredMR <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_effect.scale(top.env)

    .clean_boot(top.env)

    .clean_weights.med(top.env)

    .clean_y.YpredMR(top.env)

    if (top.env$plot) .check_plot.med(top.env)
}




#### OK  .clean_y.YpredMR ################################################

#' @rdname dot-clean_y
#' @order 6

.clean_y.YpredMR <- function(env) {

    yes10 <- ("10" %in% env$cross.world)
    yes01 <- ("01" %in% env$cross.world)


    in_c.vars <- env$c.vars

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

        if (is.null(y.c1)) y.c1 <- y.c

        if (is.null(y.c0)) y.c0 <- y.c

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




    c.vars  <- unique(c(all.vars(formula(y.c1)[[3]]),
                        all.vars(formula(y.c0)[[3]])))

    if (yes10)  c.vars <- unique(c(c.vars, all.vars(formula(y10.c)[[3]])))
    if (yes01)  c.vars <- unique(c(c.vars, all.vars(formula(y01.c)[[3]])))

    stray.c <- setdiff(c.vars, in_c.vars)

    if (length(stray.c)>0)
        stop(paste("Covariate(s)", paste(stray.c, collapse = ", "), "in outcome model(s) not found in a.c.form."))




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




#### OK  .point_est.YpredMR ####################################################

#' @rdname dot-point_est
#' @order 5

.point_est.YpredMR <- function(
    data,
    cross.world,
    effect.scale,
    a.c.form,
    a.cm.form,
    max.stabilized.wt = 30,
    output.data = FALSE,
    y.c1.form,
    y.c0.form,
    y10.c.form = NULL,
    y01.c.form = NULL,
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





    if ("10" %in% cross.world) {

        p10 <- w.dat[w.dat$.samp=="p10", ]

        y10.c.p10 <- glm(formula = y10.c.form,
                         data = p10,
                         weights = data$.f.wt,
                         family = y.family)

        pred$p10 <- predict(y10.c.p10, newdata = full, type = "response")
    }


    if ("01" %in% cross.world) {

        p01 <- w.dat[w.dat$.samp=="p01", ]

        y01.c.p01 <- glm(formula = y01.c.form,
                         data = p01,
                         weights = data$.f.wt,
                         family = y.family)

        pred$p01 <- predict(y01.c.p01, newdata = full, type = "response")
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




#### OK  .plot_YpredMR #########################################################

#' @rdname dot-plot_w.dat
#' @order 6

.plot_YpredMR <- function(w.dat,
                          cross.world,
                          c.vars,
                          m.vars,
                          c.std,
                          m.std) {

    c(.plot_wt_dist(w.dat),
      .plot_balance.YpredMR(w.dat = w.dat,
                            cross.world = cross.world,
                            c.vars = c.vars,
                            m.vars = m.vars,
                            c.std = c.std,
                            m.std = m.std))
}




#### OK  .plot_balance.YpredMR #################################################

#' @rdname dot-plot_balance
#' @order 6

.plot_balance.YpredMR <- function(w.dat,
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
        smd.dat <- smd.dat[smd.dat$contrast %in% c("p10 - p00",
                                                   "p11 - p01"), ]
    } else if (cross.world=="10") {
        smd.dat <- smd.dat[smd.dat$contrast=="p10 - p00", ]
    } else if (cross.world=="01") {
        smd.dat <- smd.dat[smd.dat$contrast=="p11 - p01", ]
    }

    smd.dat$contrast <- as.character(smd.dat$contrast)

    smd.dat$contrast <- ifelse(smd.dat$contrast=="p11 - p01",
                               "p01 - p11",
                               smd.dat$contrast)

    smd.dat$mean.diff <- ifelse(smd.dat$contrast=="p01 - p11",
                                -smd.dat$mean.diff,
                                smd.dat$mean.diff)

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












































