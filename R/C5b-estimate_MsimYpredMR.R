


# TODO: combine .clean_m.MsimYpred and .clean_m.MsimYpredMR into one file, with MR argument



#### estimate_MsimYpredMR #################################################

#' Estimator MsimYpredMR
#'
#' Function that implements estimator MsimYpredMR
#' @inheritParams estimate_MsimYpred
#' @inheritParams weights_med
#' @export

estimate_MsimYpredMR <- function(
    data,
    s.wt.var     = NULL,
    cross.world  = "10",
    effect.scale = "additive",

    a.c.form,
    a.cm.form,
    max.stabilized.wt = 30,

    m.c1.form = NULL,
    m.c0.form = NULL,
    m.c.form  = NULL,

    m.dist    = NULL,

    y.c1.form = NULL,
    y.c0.form = NULL,
    y.c.form  = NULL,

    y.cm1.form = NULL,
    y.cm0.form = NULL,
    y.cm.form  = NULL,

    y.link     = "identity",

    point.reps = 100,

    boot.num      = 999,
    boot.seed     = NULL,
    boot.method   = "cont-wt",
    boot.stratify = TRUE,

    plot = TRUE,
    c.std = NULL,
    m.std = NULL,
    c.order    = NULL,
    m.order    = NULL

) {

    # CLEAN INPUTS

    c.vars <- m.vars <- m.family <- y.family <- NULL

    .prep_MsimYpredMR()

    key.inputs <- mget(c("effect.scale",
                         "cross.world",
                         "a.c.form", "a.cm.form",
                         "max.stabilized.wt",
                         "m.vars",
                         "m.c1.form", "m.c0.form",
                         "m.family",
                         "y.c1.form", "y.c0.form",
                         "y.cm1.form", "y.cm0.form",
                         "y.family",
                         "point.reps"))


    # POINT ESTIMATION

    if (!plot) {

        estimates <- do.call(".point_est.MsimYpredMR",
                             c(key.inputs, list(data        = data,
                                                output.data = FALSE)))
    } else {

        tmp <- do.call(".point_est.MsimYpredMR",
                       c(key.inputs, list(data        = data,
                                          output.data = TRUE)))

        estimates <- tmp$estimates

        plots <- .plot_med(w.dat = tmp$w.dat,
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
                             FUN        = ".point_est.MsimYpredMR",
                             FUN.inputs = c(key.inputs,
                                            list(boot = TRUE)))

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




#### .prep_MsimYpredMR ####################################################

#' @rdname dot-prep
#' @order 12

.prep_MsimYpredMR <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_effect.scale(top.env)

    .clean_boot(top.env)

    .clean_weights.med(top.env)

    .clean_m.MsimYpredMR(top.env)

    .clean_y.psYpredMR(top.env)
}




#### .clean_m.MsimYpredMR ##################################################

#' @rdname dot-clean_m
#' @order 2

.clean_m.MsimYpredMR <- function(env) {

    yes10 <- ("10" %in% env$cross.world)
    yes01 <- ("01" %in% env$cross.world)

    m.c1   <- env$m.c1.form
    m.c0   <- env$m.c0.form
    m.c    <- env$m.c.form
    m.dist <- env$m.dist

    in_c.vars <- env$c.vars
    in_m.vars <- env$m.vars


    if (is.null(m.c)) {

        if (yes10 && is.null(m.c0))
            stop("Must specify either m.c0.form or m.c.form.")

        if (yes01 && is.null(m.c1))
            stop("Must specify either m.c1.form or m.c.form.")

    } else {

        if (yes10 && is.null(m.c1)) env$m.c1.form <- m.c1 <- m.c
        if (yes01 && is.null(m.c0)) env$m.c0.form <- m.c0 <- m.c
    }



    if (yes10) m0.vars <- sapply(m.c0, function(z) all.vars(formula(z)[[2]]))
    if (yes01) m1.vars <- sapply(m.c1, function(z) all.vars(formula(z)[[2]]))


    if (yes10 && yes01 && !setequal(m0.vars, m1.vars)) {
        stop("mediators in m.c1.form and m.c0.form not the same")

    } else if (yes10) { m.vars <- m0.vars
    } else if (yes01) { m.vars <- m1.vars
    }

    if (!setequal(m.vars, in_m.vars))
        stop("Mediators specified in mediator model(s) do not coincide with mediators based on a.cm.form and a.c.form.")



    c.vars <- NULL

    if (yes10)
        c.vars <- c(c.vars,
                    unlist(lapply(m.c0, function(z) all.vars(formula(z)))))
    if (yes01)
        c.vars <- c(c.vars,
                    unlist(lapply(m.c1, function(z) all.vars(formula(z)))))

    c.vars <- setdiff(unique(c.vars), m.vars)

    if (!setequal(c.vars, in_c.vars))
        stop("Covariates specified in mediator model(s) do not coincide with covariates in a.c.form.")



    if (is.null(m.dist))  stop("m.dist must be provided")


    m.family <- list()

    for (i in 1:length(m.dist)) {

        if (tolower(m.dist[[i]]) %in% c("normal", "gaussian",
                                        "continuous", "linear", "identity")) {
            m.family[[i]] <- "gaussian"

        } else if (tolower(m.dist[[i]]) %in% c("binary", "bernouli")) {

            m.family[[i]] <- "quasibinomial"

        } else
            stop("m.dist ", m.dist[[i]], " not recognized or supported. Only handles normal and binary simulation for now.")

    }

    # TODO: add both here and for MsimYpred checks whether specific M variables are binary or continuous

    env$m.family <- m.family


}




#### .point_est.MsimYpredMR ################################################

#' @rdname dot-point_est
#' @order 11

.point_est.MsimYpredMR <- function(
    data,
    cross.world,
    effect.scale,
    a.c.form,
    a.cm.form,
    max.stabilized.wt,
    m.vars,
    m.c1.form,
    m.c0.form,
    m.family,
    y.c1.form,
    y.c0.form,
    y.cm1.form,
    y.cm0.form,
    y.family,
    output.data = FALSE, # this is to work nice with boot function
    #                      (to revisit later)
    boot = FALSE,
    point.reps = 100
) {

    w.dat <- .compute_weights.med(
        data              = data,
        cross.world       = cross.world,
        a.c.form          = a.c.form,
        a.cm.form         = a.cm.form,
        max.stabilized.wt = max.stabilized.wt
    )


    po.means <-
        .reg_Ypred(dat       = w.dat[w.dat$.samp %in% c("p11", "p00"), ],
                   y.c1.form = y.c1.form,
                   y.c0.form = y.c0.form,
                   y.family  = y.family,
                   robust    = TRUE)

    if ("10" %in% cross.world)

        po.means[["y10.mean"]] <- .crw_MsimYpred(
            dat       = w.dat[w.dat$.samp %in% c("p10", "p00"), ],
            m.vars    = m.vars,
            m.c.form  = m.c0.form,
            m.family  = m.family,
            y.cm.form = y.cm1.form,
            y.family  = y.family,
            crw       = "10",
            robust    = TRUE,
            boot      = boot,
            reps      = point.reps)


    if ("01" %in% cross.world)
        po.means[["y01.mean"]] <- .crw_MsimYpred(
            dat       = w.dat[w.dat$.samp %in% c("p01", "p11"), ],
            m.vars    = m.vars,
            m.c.form  = m.c1.form,
            m.family  = m.family,
            y.cm.form = y.cm0.form,
            y.family  = y.family,
            crw       = "01",
            robust    = TRUE,
            boot      = boot,
            reps      = point.reps)


    estimates <- .po.means_to_effects(po.means)

    if (!output.data) return(estimates)

    mget(c("estimates", "w.dat"))

}















































