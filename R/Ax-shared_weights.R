



#### OK  .setup_data ###########################################################

#' Initial data preparation
#'
#' Internal function called by \code{.prep_get_weights.} and \code{.prep_estimate_} functions
#' @inheritParams env-block
#' @details Modify \code{data} in \code{env}, front-appending variables: \code{.row}, \code{.samp}, \code{.s.wt}, \code{.w.wt}, \code{.f.wt}, \code{.a}, \code{.y}. Set \code{.row} to row number, \code{.s.wt} to values of variable named in object \code{s.wt.var} in \code{env}, and others to NA (placeholders for sample indicator, distribution-morphing weights, and final weights).
#' @noRd

.setup_data <- function(env) {

    env$data <- cbind(.row  = 1:nrow(env$data),
                      .samp = NA,
                      .s.wt = NA,
                      .w.wt = NA,
                      .f.wt = NA,
                      .a = NA,
                      .y = NA,
                      env$data)

    s.wt.var <- env$s.wt.var

    if (is.null(s.wt.var)) {
        env$data$.s.wt <- 1
    } else if (!is.character(s.wt.var) || length(s.wt.var)>1) {
        stop("s.wt.var must be a variable name, ie a single character string.")
    } else if (s.wt.var=="") {
        env$data$.s.wt <- 1
    } else if (!s.wt.var %in% names(env$data)) {
        stop(paste("Variable", s.wt.var, "not found in dataset"))
    } else
        env$data$.s.wt <- env$data[, s.wt.var]
}




#### OK  .clean_cross.world ####################################################

#' Clean cross.world input
#'
#' Internal function called by \code{.prep_get_weights.} and \code{.prep_estimate_} functions
#' @inheritParams env-block
#' @return If \code{cross.world} value is "both", change to c("10", "01").
#' @noRd

.clean_cross.world <- function(env) {

    cross.world <- env$cross.world

    if ("both" %in% cross.world) {

        env$cross.world <- c("10", "01")

    } else if (!all(cross.world %in% c("10", "01")))

        stop("Please use allowed cross.world values: \"10\", \"01\", \"both\".")

}




#### OK  stabilize_weight ######################################################

#' Standardize a weight variable by group
#' @param weight A weight vector to be standardized
#' @param s.weight Optional. A vector of sampling weight.
#' @param group Optional. A character vector indicating group membership.


stabilize_weight <- function(weight,
                             group = NULL,
                             s.weight = NULL) {

    if (is.null(s.weight)) s.weight <- rep(1, length(weight))

    if (is.null(group))
        return(weight / .wtd_mean(weight, s.weight))

    tmp <- data.frame(id = 1:length(weight),
                      wt = weight,
                      s.wt = s.weight)

    tmp <- lapply(unique(group), function(z) {
        g.dat <- tmp[group==z, ]

        g.dat$wt <- g.dat$wt / .wtd_mean(g.dat$wt, g.dat$s.wt)

        g.dat[, c("id", "wt")]
    })

    tmp <- do.call("rbind", tmp)
    tmp <- tmp[order(tmp$id), ]

    tmp$wt
}




#### OK  .plot_wt_dist #########################################################

#' .plot_wt_dist
#'
#' Internal function called by \code{weights_\*()} and \code{estimate_\*()} to make weight distribution plot(s).
#' @param w.dat The weighted dataset.
#' @param point.alpha,jitter.width Graphical parameters for \code{geom_point()}.
#' @importFrom ggplot2 ggplot aes geom_jitter geom_violin labs theme_bw theme
#' @importFrom rlang .data
#' @return If constant sampling weights, plot the densities of distribution-morphing weights in the pseudo samples
#' @return If non-constant sampling weights, output two plots, one of the distribution-morphing weights, one of the final weights (product of sampling and distribution-morphing weights)
#' @name dot-plot_wt_dist
NULL

#' @rdname dot-plot_wt_dist
#' @order 1

.plot_wt_dist <- function(
    w.dat,
    point.alpha = .1,
    jitter.width = .3
) {

    w.dat$.w.wt <- stabilize_weight(weight   = w.dat$.w.wt,
                                    group    = w.dat$.samp,
                                    s.weight = w.dat$.s.wt)

    if (is_constant(w.dat$.s.wt)) {

        dat <- w.dat[, c(".samp", ".w.wt")]

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





#### .plot_balance ### PERHAPS DEGENERATE ######################################


#' Plot balance
#'
#' Internal function that makes balance plot for mediation weights and TE weights.
#' @param w.dat abc
#' @param c.vars Names of covariates.
#' @param m.vars Names of mediators.
#' @param c.vars.std Names of covariates whose mean differences are to be standardized.
#' @param m.vars.std Names of mediators whose mean differences are to be standardized.
#' @param full.label Whether to add "(anchor)" and "(for <effect>)" notes to plot labels to help interpretation. Defaults to TRUE
#' @importFrom ggplot2 ggplot aes geom_vline geom_point scale_color_manual scale_shape_manual labs theme_bw facet_wrap xlim
#' @importFrom rlang .data
#' @noRd

.plot_balance <- function(w.dat,
                          c.vars,
                          m.vars = NULL,
                          c.vars.std,
                          m.vars.std,
                          full.label = TRUE) {


    # dummy code categorical variables and modify plot variable lists
    tmp <- .make_dummies(w.dat,
                         columns = c.vars,
                         output.names = TRUE,
                         warning = FALSE)
    w.dat <- tmp$data
    plot.c <- tmp$columns;  rm(tmp)

    if (!is.null(m.vars)) {
        tmp <- .make_dummies(w.dat,
                             columns = m.vars,
                             output.names = TRUE,
                             warning = FALSE)
        w.dat <- tmp$data
        plot.m <- tmp$columns;  rm(tmp)
    } else
        plot.m <- NULL



    # compute SMDs
    smd.dat <- .get_c.smd(w.dat      = w.dat,
                          vars        = plot.c,
                          standardize = intersect(plot.c, c.vars.std))

    if (!is.null(plot.m))
        smd.dat <- rbind(smd.dat,
                         .get_m.smd(w.dat      = w.dat,
                                    vars        = plot.m,
                                    standardize = intersect(plot.m,
                                                            m.vars.std)))


    plot.vars <- c(plot.c, plot.m)
    plot.vars <- plot.vars[length(plot.vars):1]

    smd.dat$variable <- factor(smd.dat$variable, levels = plot.vars)

    p <- ggplot(data = smd.dat,
                aes(x = .data$mean.diff,
                    y = .data$variable)) +
        geom_vline(xintercept = 0,
                   color = "gray60")

    if ("mediator" %in% unique(smd.dat$var.type)) {
        p <- p +
            geom_point(aes(color = .data$var.type,
                           shape = .data$contrast.type),
                       fill = "white",
                       size = 1.5,
                       stroke = .5) +
            scale_color_manual(name = "", values = c(1, "magenta"))

    } else {
        p <- p +
            geom_point(aes(shape = .data$contrast.type),
                       fill = "white",
                       size = 1.5,
                       stroke = .5)
    }

    p <- p +
        labs(x = "differences in means", y = "") +
        scale_shape_manual(name = "", values = c(21, 19)) +
        theme_bw() +
        xlim(min(-.3, min(smd.dat$mean.diff)),
             max(.3, max(smd.dat$mean.diff)))

    if (full.label) {
        p <- p + facet_wrap(~ contrast.labs, ncol = 3)
    } else {
        p <- p + facet_wrap(~ contrast, ncol = 3)
    }

    p
}




