

#### estimate_effects ##########################################################

# TODO: figure out how to cite the paper and refer to the paper
# TODO: fill out the @return documentation

#' Estimate natural (in)direct effects using an estimator of choice
#'
#' A wrapper function that calls on estimator-specific functions.
#' @param data A data frame.
#' @param s.wt.var Optional, name of variable containing sampling weights.
#' @param estimator The estimator of choice. See the paper. -> TOADD citation. Defaults to "Y2pred".
#' @param cross.world The cross-world condition involved in the effect decomposition of choice. Should be "10" if want the (NDE0, NIE1) pair, "01" if want the (NIE0, NDE1) pair, or "both" if want both decompositions.
#' @param effect.scale The scale of effect of choice. Defaults to "MD" (i.e., mean/risk difference or additive). If outcome is non-negative, also allows "mean ratio" (which could also be specified as "ratio", "MR", "risk ratio", "rate ratio", "RR"). If outcome is binary or bounded within the (0,1) interval, also allows "odds ratio" (which could also be specified as "OR").
#' @param boot.num Number of bootstrap samples used for interval estimation, defaults to 999. If just want point estimate, set to 0.
#' @param boot.seed Optional, specify bootstrap seed for reproducibility.
#' @param boot.stratify Whether bootstrap samples are drawn stratified by treatment variable. Defaults to TRUE.
#' @param boot.method Method for drawing bootstrap samples. Options: "cont-wt" for continuous weights bootstrap, and "resample" for bootstrap by simple resampling (i.e., integer weights bootstrap).
#' @param ... Inputs specific to the estimator.
#' @return A list of objects, including\itemize{
#' \item{}{}
#' }
#' @export

estimate_effects <- function(
    data,
    s.wt.var = NULL,

    estimator = "Y2predR",
    cross.world = "10",
    effect.scale = "MD",

    boot.num = 999,
    boot.seed = NULL,
    boot.method = "cont-wt",
    boot.stratify = TRUE,
    ...

) {

    estimator <- .clean_estimator(estimator = estimator)

    args.specific <- .grab_args(estimator = estimator,
                                arg.list = ls(...))

    args.generic <- list(data = data,
                         s.wt.var = s.wt.var,
                         cross.world = cross.world,
                         effect.scale = effect.scale,
                         boot.num = boot.num,
                         boot.seed = boot.seed,
                         boot.method = boot.method,
                         boot.stratify = TRUE)

    do.call(paste0("estimate_", estimator),
            c(args.generic, args.specific))

}





################################################################################

# TODO: handle the remnant below

# remnant_estimate_effects <- function() {
#     if (is.null(s.wt)) { data$s.wt <- 1
#     } else             { data$s.wt <- data[, s.wt];  s.wt <- "s.wt"
#     }
#
#     arg.list <- clean_args(estimator = estimator,
#                            cross.world = cross.world,
#                            ...)
#
#     arg.list[["effect.scale"]] <- clean_scale(effect.scale)
#
#
#     # point estimation
#     point.args <- arg.list
#     point.args[["data"]] <- data
#
#     point <- do.call(paste0("point_", estimator), point.args)
#
#
#     # bootstrap for CI and SE
#     if (boot.stratify) { boot.stratify <- arg.list$a.var
#     } else             { boot.stratify <- NULL
#     }
#
#     if (is.null(boot.seed)) boot.seed <- get_seed()
#
#
#     ci.se <- boot_ci.se(data = data,
#                         boot.stratify = boot.stratify,
#                         boot.num = boot.num,
#                         boot.seed = boot.seed,
#                         boot.method = boot.method,
#                         FUN = paste0("point_", estimator),
#                         arg.list = arg.list)
#
#
#
#     list(effects = cbind(estimate = point, ci.se),
#          boot.seed = boot.seed)
# }









##### .clean_estimator #########################################################

#' Clean typos in estimator name
#'
#' @param estimator The estimator name as typed by user.
#' @return The right estimator name (correcting anticipated typos).

.clean_estimator <- function(estimator) {

    if (tolower(estimator) %in%
        c("psypred"))
        return("psYpred")

    if (tolower(estimator) %in%
        c("psypredr", "psypredmr", "psypred.r", "psypred.mr"))
        return("psYpredMR")

    if (tolower(estimator) %in%
        c("ypred"))
        return("Ypred")

    if (tolower(estimator) %in%
        c("ypredr", "ypredmr", "ypred.r", "ypred.mr"))
        return("YpredMR")

    if (tolower(estimator) %in%
        c("msimypred"))
        return("MsimYpred")

    if (tolower(estimator) %in%
        c("msimypredr", "msimypredmr", "msimypred.r", "msimypred.mr"))
        return("MsimYpredMR")

    if (tolower(estimator) %in%
        c("y2pred"))
        return("Y2pred")

    if (tolower(estimator) %in%
        c("y2predr", "y2predmr", "y2pred.r", "y2pred.mr"))
        return("Y2predR")

    if (tolower(estimator) %in%
        c("ndepred"))
        return("NDEpred")

    if (tolower(estimator) %in%
        c("ndepredr", "ndepredmr", "ndepred.r", "ndepred.mr"))
        return("NDEpredMR")

    if (tolower(estimator) %in%
        c("wt-cadj1", "wt.cadj1", "wtcadj1"))
        return("wtCadj1")

    if (tolower(estimator) %in%
        c("wt-cadj2", "wt.cadj2", "wtcadj2"))
        return("wtCadj2")

    if (tolower(estimator) %in%
        c("wp-cadj1", "wp.cadj1", "wpcadj1"))
        return("wpCadj1")

    if (tolower(estimator) %in%
        c("wp-cadj2", "wp.cadj2", "wpcadj2"))
        return("wpCadj2")

    if (tolower(estimator) %in%
        c("wp-mr-cadj1", "wp.mr-cadj1",
          "wp-mr.cadj1", "wp.mr.cadj1", "wpmrcadj1",
          "wp-r-cadj1", "wp.r-cadj1",
          "wp-r.cadj1", "wp.r.cadj1", "wprcadj1"))
        return("wpMRCadj1")

    if (tolower(estimator) %in%
        c("wp-mr-cadj2", "wp.mr-cadj2",
          "wp-mr.cadj2", "wp.mr.cadj2", "wpmrcadj2",
          "wp-r-cadj2", "wp.r-cadj2",
          "wp-r.cadj2", "wp.r.cadj2", "wprcadj2"))
        return("wpMRCadj2")

    stop("Estimator not recognized. Please refer to documentation.")


}



#### .grab_args ################################################################

# TODO: should this be called grab inputs (or grab specific inputs)? Does args refer to names of objects or to the objects themselves?

#' Filter a list of objects for those that are required arguments for an estimator
#'
#' @param estimator Name of estimator.
#' @param arg.list A list of objects containing (but not necessarily limited to) the inputs needed for the estimator.
#' @noRd

.grab_args <- function(estimator,
                       arg.list) {

    required <- .get_fun_arg.names(paste0(".point_", estimator))

    arg.list[names(arg.list) %in% required]
}



#### .get_fun_arg.names ########################################################

#' Get required argument names from a function
#'
#' @param fun.name Name of a function.
#' @return A vector of names of the required arguments of the function
#' @noRd

.get_fun_arg.names <- function(fun.name) {
    required <- names(as.list(args(fun.name)))
    required[! required %in% c("...", "")]
}











