



#### OK  .clean_effect.scale ###################################################

#' @inheritParams env-block
#' @details Modify object \code{effect.scale} in \code{env} to one of three values: "MD", "MR" and "OR".
#' @noRd

.clean_effect.scale <- function(env) {

    effect.scale <- env$effect.scale

    if (effect.scale %in% c("additive",
                            "mean difference",
                            "risk difference", "RD")) {

        env$effect.scale <- "MD"

    } else if (effect.scale %in% c("ratio",
                                   "mean ratio",
                                   "risk ratio", "rate ratio", "RR")) {

        env$effect.scale <- "MR"

    } else if (effect.scale=="odds ratio") {

        env$effect.scale <- "OR"

    } else if (!effect.scale %in% c("MD", "MR", "OR"))

        stop(paste("Effect scale",
                   effect.scale,
                   "is not understood or not supported."))
}




#### OK  .clean_boot ###########################################################

#' @inheritParams env-block
#' @details If \code{boot.seed} in \code{env} is NULL, draw a seed. Modify \code{boot.method} in \code{env} to one of two values: "cont-wt" and "resample".
#' @noRd

.clean_boot <- function(env) {

    if (env$boot.num==0)
        return()


    if (is.null(env$boot.seed))  env$boot.seed <- .get_seed()


    if (env$boot.method %in% c("continuous weights",
                               "dirichlet")) {
        env$boot.method <- "cont-wt"

    } else if (env$boot.method %in% c("Efron",
                                      "classic",
                                      "multinomial")) {
        env$boot.method <- "resample"

    } else if (!env$boot.method %in% c("cont-wt", "resample"))
        stop(paste("Boot method",
                   env$boot.method,
                   "not understood. Please specify either \"cont-wt\" (for continuous bootstrap weights method) or \"resample\" (for the classic resampling method)."))

}




#### OK  .get_means.and.effects ################################################

#' Compute potential outcome means from pseudo samples
#'
#' Compute potential outcome means from pseudo samples
#' @param w.dat Pseudo sample data
#' @param effect.scale The contrast type
#' @return A named vector of potential outcome means and effects
#' @noRd

.get_means.and.effects <- function(w.dat,
                                   effect.scale) {

    ps.names <- unique(w.dat$.samp)

    po.means <- lapply(ps.names, function(z) {
        dat <- w.dat[w.dat$.samp==z, ]
        .wtd_mean(dat$.y, dat$.f.wt)
    })

    names(po.means) <- ps.names



    effects <- c(TE = .get_contrast(po.means$p11, po.means$p00,
                                    type = effect.scale))

    if ("p10" %in% names(po.means)) {
        effects <- c(effects,
                     NDE0 = .get_contrast(po.means$p10, po.means$p00,
                                          type = effect.scale))
        effects <- c(effects,
                     NIE1 = .get_contrast(po.means$p11, po.means$p10,
                                          type = effect.scale))
    }

    if ("p01" %in% names(po.means)) {
        effects <- c(effects,
                     NIE0 = .get_contrast(po.means$p01, po.means$p00,
                                          type = effect.scale))
        effects <- c(effects,
                     NDE1 = .get_contrast(po.means$p11, po.means$p01,
                                          type = effect.scale))
    }

    c(unlist(po.means), effects)
}



#### OK  .get_contrast #########################################################

#' Contrast two potential outcome means
#'
#' @param a The mean of the first potential outcome.
#' @param b The mean of the second potential outcome.
#' @param type The type of contrast. Options allowed: "MD", "MR", "OR".
#' @return The contrast value.
#' @keywords internal

.get_contrast <- function(a, b, type) {

    if (type=="MD")
        return(a - b)


    if (type=="MR") {
        if (a<0 || b<0) {
            stop("An outcome mean is negative. Ratio effect scale is not allowed.")
        } else if (b<.000001) {
            stop("A denominator outcome mean is zero. Ratio effect scale is not appropriate.")
        } else
            return(a / b)
    }


    if (type=="OR") {
        if (a<=0 || b<=0 || a>=1 || b>=1) {
            stop("An outcome mean is outside of the (0,1) interval. Odds ratio effect scale is not appropriate")
        } else
            return((a / (1-a)) / (b / (1-b)))
    }

}










