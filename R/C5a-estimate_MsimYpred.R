
# TODO: reorder arguments of functions: boot inputs last, y.cm after y.cm1 and y.cm0, same default for effect.scale.

# TODO: revisit output.data argument in .point_est for estimators that don't require weighting

# TODO: bundle .point and .boot in one convenient sample





#### estimate_MsimYpred ######################################################

#' Estimator MsimYpred
#'
#' Function that implements estimator MsimYpred
#' @inheritParams estimate_effects
#' @inheritParams estimate_Y2pred
#' @param m.c1.form,m.c0.form blah
#' @param m.c.form blah
#' @param m.dist blah
#' @param point.reps blah
#' @export

estimate_MsimYpred <- function(
    data,
    s.wt.var     = NULL,
    cross.world  = "10",
    effect.scale = "additive",

    a.var,

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

    boot.num      = 999,
    boot.seed     = NULL,
    boot.method   = "cont-wt",
    boot.stratify = TRUE,

    point.reps = 100

) {

    # CLEAN INPUTS

    c.vars <- m.vars <- m.family <- y.family <- NULL

    .prep_MsimYpred()

    key.inputs <- mget(c("effect.scale",
                         "cross.world",
                         "m.vars",
                         "m.c1.form", "m.c0.form",
                         "m.family",
                         "y.c1.form", "y.c0.form",
                         "y.cm1.form", "y.cm0.form",
                         "y.family",
                         "point.reps"))


    # POINT ESTIMATION

    estimates <- do.call(".point_est.MsimYpred", c(key.inputs,
                                                   list(data = data)))


    # BOOTSTRAP

    if (boot.num==0) return(estimates)


    ci.se <- .boot_ci.se(data       = data,
                         stratify   = boot.stratify,
                         boot.num   = boot.num,
                         seed       = boot.seed,
                         method     = boot.method,
                         FUN        = ".point_est.MsimYpred",
                         FUN.inputs = c(key.inputs,
                                        list(boot = TRUE)))

    estimates <- cbind(estimate = estimates,
                       ci.se)


    mget(c("estimates", "boot.seed"))

}




#### .prep_MsimYpred ######################################################

#' @rdname dot-prep
#' @order 11

.prep_MsimYpred <- function() {

    top.env <- parent.frame()

    .setup_data(top.env)

    .clean_cross.world(top.env)

    .clean_effect.scale(top.env)

    .clean_boot(top.env)

    .clean_a.var(top.env)

    .clean_m.MsimYpred(top.env)

    .clean_y.psYpredMR(top.env)
}




#### .clean_m.MsimYpred ###################################################

#' Internal: clean inputs related to mediator model(s)
#'
#' @inheritParams env-block
#' @name dot-clean_m
#' @keywords internal
NULL

#' @rdname dot-clean_m
#' @order 1

.clean_m.MsimYpred <- function(env) {

    yes10 <- ("10" %in% env$cross.world)
    yes01 <- ("01" %in% env$cross.world)

    m.c1   <- env$m.c1.form
    m.c0   <- env$m.c0.form
    m.c    <- env$m.c.form
    m.dist <- env$m.dist


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

    } else if (yes10) { env$m.vars <- m.vars <- m0.vars
    } else if (yes01) { env$m.vars <- m.vars <- m1.vars
    }


    c.vars <- NULL

    if (yes10)
        c.vars <- c(c.vars,
                    unlist(lapply(m.c0, function(z) all.vars(formula(z)))))
    if (yes01)
        c.vars <- c(c.vars,
                    unlist(lapply(m.c1, function(z) all.vars(formula(z)))))

    env$c.vars <- c.vars <- setdiff(unique(c.vars), m.vars)



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

    # TODO: add checks whether specific M variables are binary or continuous

    env$m.family <- m.family

}









#### .point_est.MsimYpred ################################################

#' @rdname dot-point_est
#' @order 10

.point_est.MsimYpred <- function(
    data,
    cross.world,
    effect.scale,
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


    po.means <- .reg_Ypred(dat       = data,
                           y.c1.form = y.c1.form,
                           y.c0.form = y.c0.form,
                           y.family  = y.family,
                           robust    = FALSE)

    if ("10" %in% cross.world)
        po.means[["y10.mean"]] <- .crw_MsimYpred(dat       = data,
                                                 m.vars    = m.vars,
                                                 m.c.form  = m.c0.form,
                                                 m.family  = m.family,
                                                 y.cm.form = y.cm1.form,
                                                 y.family  = y.family,
                                                 crw       = "10",
                                                 robust    = FALSE,
                                                 boot      = boot,
                                                 reps      = point.reps)


    if ("01" %in% cross.world)
        po.means[["y01.mean"]] <- .crw_MsimYpred(dat       = data,
                                                 m.vars    = m.vars,
                                                 m.c.form  = m.c1.form,
                                                 m.family  = m.family,
                                                 y.cm.form = y.cm0.form,
                                                 y.family  = y.family,
                                                 crw       = "01",
                                                 robust    = FALSE,
                                                 boot      = boot,
                                                 reps      = point.reps)


    .po.means_to_effects(po.means)

}




#### .simulate ############################################################

#' Internal: simulate based on a model fit
#'
#' Similar to \code{base::simulate()}, but allows a \code{newdata} argument, similar to that of \code{predict()}.
#' @importFrom stats family rbinom rnorm sigma
#' @keywords internal

.simulate <- function(object, ...) {

    if (!any(class(object) %in% c("glm", "lm", "svyglm")))
        stop("The function does not support simulating based on an object of this class: ", paste(class(object), collapse = ", "))


    newdata <- NULL


    .extract_dots(...)


    if (is.null(newdata)) { means <- predict(object, type = "response")
    } else                { means <- predict(object, type = "response",
                                             newdata = newdata)
    }



    if (family(object)[[1]] %in% c("quasibinomial", "binomial")) {

        sims <- rbinom(n = length(means), size = 1, prob = means)

    } else if (family(object)[[1]] == "gaussian") {

        sd <- sigma(object) / sqrt(mean(object$prior.weights))

        sims <- rnorm(n = length(means), mean = means, sd = sd)
    }

    sims


}


.extract_dots <- function(...) {

    # version using with() that does not pass R CMD check
    # (remove ... in function call)
    #
    # with(parent.frame(), {
    #     dots <- list(...)
    #
    #     for (i in 1:length(dots))
    #         assign(x = names(dots)[i], value = dots[[i]])
    #
    #     rm(dots)
    # })


    dots <- list(...)

    for (i in 1:length(dots))
        assign(x     = names(dots)[i],
               value = dots[[i]],
               envir = parent.frame())


}



#### .regYpred ############################################################

#' @noRd

.reg_Ypred <- function(dat,
                       y.c1.form,
                       y.c0.form,
                       y.family,
                       robust = FALSE) {

    data <- NULL

    if (robust)  dat$wt <- dat$.f.wt  else  dat$wt <- dat$.s.wt


    dat1 <- dat[dat$.a==1, ]
    dat0 <- dat[dat$.a==0, ]


    y.c1.fit <- glm(formula = y.c1.form,
                    data    = dat1,
                    weights = data$wt,
                    family  = y.family)

    y.c0.fit <- glm(formula = y.c0.form,
                    data    = dat0,
                    weights = data$wt,
                    family  = y.family)

    y11.pred <- predict(y.c1.fit, newdata = dat, type = "response")
    y00.pred <- predict(y.c0.fit, newdata = dat, type = "response")

    y11.mean <- .wtd_mean(y11.pred, dat$.s.wt)
    y00.mean <- .wtd_mean(y00.pred, dat$.s.wt)

    list(y11.mean = y11.mean,
         y00.mean = y00.mean)

}




#### .crw_msimYpred #####################################################

#' @noRd
#' @importFrom stats glm

.crw_MsimYpred <- function(dat,
                           m.vars,
                           m.c.form,
                           m.family,
                           y.cm.form,
                           y.family,
                           crw,
                           robust = FALSE,
                           boot   = FALSE,
                           reps) {

    data <- NULL


    if (robust)  dat$wt <- dat$.f.wt  else  dat$wt <- dat$.s.wt



    if (crw=="10")        { mdat <- dat[dat$.a==0,];  ydat <- dat[dat$.a==1,]
    } else if (crw=="01") { mdat <- dat[dat$.a==1,];  ydat <- dat[dat$.a==0,]
    }



    m.mod <- lapply(1:length(m.c.form), function(z) {

        glm(formula = m.c.form[[z]],
            data    = mdat,
            weights = data$wt,
            family  = m.family[[z]])
    })


    if (boot) reps <- 1

    y10.mean.reps <- sapply(1:reps, function(z) {

        for (i in 1:length(m.vars))
            dat[, m.vars[i]] <- .simulate(m.mod[[i]], newdata = dat)



        y.mod <- glm(formula = y.cm.form,
                     data    = ydat,
                     weights = data$wt,
                     family  = y.family)

        y10.pred <- predict(y.mod, newdata = dat, type = "response")

        .wtd_mean(y10.pred, dat$.s.wt)

    })

    mean(y10.mean.reps)

}



#' @noRd

.po.means_to_effects <- function(means) {

    effects <- c(TE = means[["y11.mean"]] - means[["y00.mean"]])

    if ("y10.mean" %in% names(means))
        effects <- c(effects,
                     NDE0 = means[["y10.mean"]] - means[["y00.mean"]],
                     NIE1 = means[["y11.mean"]] - means[["y10.mean"]])

    if ("y01.mean" %in% names(means))
        effects <- c(effects,
                     NIE0 = means[["y01.mean"]] - means[["y00.mean"]],
                     NDE1 = means[["y11.mean"]] - means[["y01.mean"]])

    c(unlist(means), effects)
}







































