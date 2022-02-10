
# TODO: reorder arguments of functions: boot inputs last, y.cm after y.cm1 and y.cm0, same default for effect.scale.

# TODO: rename .prep_estimate_wtd() to .prep_wtd()

# TODO: make reshape_gather() an internal function

# TODO: simplify function names .clean_y.forms to .clean_y

# TODO: correct order of .prep_ (letting MsimYpred in its place)

# TODO: correct order of .clean_y. (letting MsimYpred in its place)

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
#' @param m.link blah
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

    m.link    = "identity",

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
    boot.stratify = TRUE

) {

    # CLEAN INPUTS

    m.family <- y.family <- NULL

    .prep_MsimYpred()

    key.inputs <- mget(c("effect.scale",
                         "cross.world",
                         "m.c1.form", "m.c0.form",
                         "m.family",
                         "y.c1.form", "y.c0.form",
                         "y.cm1.form", "y.cm0.form",
                         "y.family"))


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
                         FUN.inputs = key.inputs)

    estimates <- cbind(estimate = estimates,
                       ci.se)

}




#### .prep_MsimYpred ######################################################

#' @rdname dot-prep
#' @order 13

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

.clean_m.MsimYpred <- function(env) {

    yes10 <- ("10" %in% env$cross.world)
    yes01 <- ("01" %in% env$cross.world)

    m.c1 <- env$m.c1.form
    m.c0 <- env$m.c0.form
    m.c  <- env$m.c.form


    if (is.null(m.c)) {

        if (yes10 && is.null(m.c1))
            stop("Must specify either m.c1.form or m.c.form.")

        if (yes01 && is.null(m.c0))
            stop("Must specify either m.c0.form or m.c.form.")
    }


    env$m.vars <- m.vars <-
        sapply(m.c, function(z) all.vars(formula(z)[[2]]))
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
    output.data = FALSE # this is to work nice with boot function
    # (to revisit later)
) {


    po.means <- .reg_Ypred(dat       = data,
                           y.c1.form = y.c1.form,
                           y.c0.form = y.c0.form,
                           y.family  = y.family,
                           robust    = FALSE)

    if ("10" %in% cross.world)
        po.means[["y10.mean"]] <- .crw_MsimYpred(dat = data,
                                                 m.vars    = m.vars,
                                                 m.c.form  = m.c0.form,
                                                 m.family  = m.family,
                                                 y.cm.form = y.cm1.form,
                                                 y.family  = y.family,
                                                 crw       = "10",
                                                 robust    = FALSE)


    if ("01" %in% cross.world)
        po.means[["y01.mean"]] <- .crw_MsimYpred(dat = data,
                                                 m.vars    = m.vars,
                                                 m.c.form  = m.c1.form,
                                                 m.family  = m.family,
                                                 y.cm.form = y.cm0.form,
                                                 y.family  = y.family,
                                                 crw       = "01",
                                                 robust    = FALSE)


    .po.means_to_effects(po.means)

}




#### .simulate ############################################################

#' @noRd
#' @importFrom stats family rbinom rnorm sigma

.simulate <- function(object, ...) {

    if (!any(class(object) %in% c("glm", "lm", "svyglm")))
        stop("The function does not support simulating based on an object of this class: ", paste(class(object), collapse = ", "))



    if ("newdata" %in% ls(...)) {
        means <- predict(object, newdata = get("newdata"), type = "response")
    } else {
        means <- predict(object, type = "response")
    }



    if (family(object)[[1]] %in% c("quasibinomial", "binomial")) {

        sims <- rbinom(n = length(means), size = 1, prob = means)

    } else if (family(object)[[1]] == "gaussian") {

        sd <- sigma(object) / sqrt(mean(object$prior.weights))

        sims <- rnorm(n = length(means), mean = means, sd = sd)
    }

    sims
}



#### .regYpred ############################################################

#' @noRd

.reg_Ypred <- function(dat,
                       y.c1.form,
                       y.c0.form,
                       y.family,
                       robust = FALSE) {

    if (robust)  dat$wt <- dat$.f.wt  else  dat$wt <- dat$.s.wt


    dat1 <- dat[dat$.a==1, ]
    dat0 <- dat[dat$.a==0, ]


    y.c1.fit <- glm(formula = y.c1.form,
                    data    = dat1,
                    weights = dat1$wt,
                    family  = y.family)

    y.c0.fit <- glm(formula = y.c0.form,
                    data    = dat0,
                    weights = dat0$wt,
                    family  = y.family)

    y11.pred <- predict(y.c1.fit, newdata = dat, type = "response")
    y00.pred <- predict(y.c0.fit, newdata = dat, type = "response")

    y11.mean <- .wtd_mean(y11.pred, dat$.s.wt)
    y00.mean <- .wtd_mean(y00.pred, dat$.s.wt)

    list(y11.mean = y11.mean,
         y00.mean = y00.mean)

}


#' @noRd

.crw_MsimYpred <- function(dat,
                           m.vars,
                           m.c.form,
                           m.family,
                           y.cm.form,
                           y.family,
                           crw,
                           robust = FALSE) {


    if (robust)  dat$wt <- dat$.f.wt  else  dat$wt <- dat$.s.wt



    if (crw=="10")        { mdat <- dat[dat$.a==0,];  ydat <- dat[dat$.a==1,]
    } else if (crw=="01") { mdat <- dat[dat$.a==1,];  ydat <- dat[dat$.a==0,]
    }



    m.mod <- lapply(1:length(m.c.form), function(z) {

        glm(formula = m.c.form[z],
            data    = mdat,
            weights = mdat$wt,
            family  = m.family[z])
    })

    for (i in 1:length(m.vars))
        dat[, m.vars[i]] <- .simulate(m.mod[[i]], newdata = dat)



    y.mod <- glm(formula = y.cm.form,
                 data    = ydat,
                 weights = ydat$wt,
                 family  = y.family)

    y10.pred <- predict(y.mod, newdata = dat, type = "response")

    .wtd_mean(y10.pred, dat$.s.wt)

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







































