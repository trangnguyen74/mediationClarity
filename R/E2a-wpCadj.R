



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

    y.cm.form,
    y.link,

    plot    = TRUE,
    c.order = NULL,
    c.std   = NULL
) {

    c.vars <- y.family <- NULL
}
