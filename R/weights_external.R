

#' @name get_med_weights
#' @description This function estimates the typical weights for estimation of natural (in)direct effects including inverse probability weights that form two pseudo samples mimicking the covariate distribution of the full sample and cross-world weights that form (one or two) pseudo samples that also mimic the mediator given covariate distribution of control units.
#' @param
#' data The sample you do analysis on
#' s.wt Sampling weights. Defaults to NULL. This argument is also used in bootstrapping (in the boot_an_estimator() function) for listing a variable that hold the bootstrap weights.
#' @return
#' @export



get_med_weights <- function(
    data,
    s.wt = NULL,
    cross.world = "10", # options: "none", "10", "01", "both"
    a.c.form,
    a.cm.form,

    max.stabilized.wt = 30,
    compact = FALSE

) {

    if (cross.world=="none")
        return(get_ipw_weights(data = data,
                               s.wt = s.wt,
                               a.c.form = a.c.form,
                               max.stabilized.wt = max.stabilized.wt,
                               compact = compact))



    # argument checks:
    if (!cross.world %in% c("10", "01", "both"))
        stop("Please use allowed \"cross.world\" values: \"10\", \"01\", \"both\".")

    if (!coherentA_a.c.form_a.cm.form(a.c.form, a.cm.form))
        stop("Treatment variable is not the same in formulas a.c.form and a.cm.form.")

    if (!coherentC_a.c.form_a.cm.form(a.c.form, a.cm.form))
        stop("Some C variables in a.c.form do not appear in a.cm.form.")



    # data check: treatment coded as 0/1
    a.var <- as.character(formula(a.c.form)[[2]])

    if (!is_binary01(data[, a.var]))
        stop("Variable named in \"a.var\" must be numeric and in binary 0/1 form.")



    # working (and template) dataset
    dat <- data.frame(rownum = 1:nrow(data),
                      s.wt = NA,       # sampling weight
                      sample = NA,     # which pseudo sample
                      w.wt = NA,       # dist. morphing weights (omega in paper)
                      f.wt = NA,         # s.wt * w.wt
                      data)

    if (is.null(s.wt)) { dat$s.wt <- rep(1, nrow(data))
    } else             { dat$s.wt <- data[, s.wt]
    }

    ## need to instruct folks to avoid naming anything f.wt or sample and avoid naming anything that is not sampling weight s.wt --> TO DO




    # models
    a.c.fu <- glm(formula = a.c.form,
                  data    = dat,
                  weights = s.wt,
                  family  = quasibinomial)
    a.cm.fu <- glm(formula = a.cm.form,
                   data    = dat,
                   weights = s.wt,
                   family  = quasibinomial)



    # prep for w.wt truncation: max.wt values and trunc_wt() function
    if (is.null(max.stabilized.wt)) {
        max.wt <- list(control = Inf, treat   = Inf)
    } else {
        max.wt <- lapply(list(control=0, treat=1), function(z) {
            max.stabilized.wt *
                (sum(dat$s.wt) / sum(dat$s.wt * (dat[, a.var]==z)))
        })
    }

    trunc_wt <- function(vec, max.val) {
        (vec <= max.val) * vec + (vec > max.val) * max.val
    }



    # p00 and p11 pseudo samples
    p00 <- dat[dat[, a.var]==0, ];  p00$sample <- "p00"
    p11 <- dat[dat[, a.var]==1, ];  p11$sample <- "p11"

    p00$w.wt <- 1 + exp( predict(a.c.fu, newdata = p00, type = "link"))
    p11$w.wt <- 1 + exp(-predict(a.c.fu, newdata = p11, type = "link"))

    p00$w.wt <- trunc_wt(p00$w.wt, max.wt$control)
    p11$w.wt <- trunc_wt(p11$w.wt, max.wt$treat)

    p00$f.wt <- p00$s.wt * p00$w.wt
    p11$f.wt <- p11$s.wt * p11$w.wt

    out <- rbind(p00, p11)



    # cross-world weights

    if (cross.world=="10" | cross.world=="both") {

        p10 <- dat[dat[, a.var]==1, ];  p10$sample <- "p10"

        p10$w.wt <-
            exp(-predict(a.cm.fu, newdata = p10, type = "link")) *   # odds term
            (1 + exp(predict(a.c.fu, newdata = p10, type = "link"))) # inv.prob term

        p10$w.wt <- trunc_wt(p10$w.wt, max.wt$treat)

        p10$f.wt <- p10$s.wt * p10$w.wt

        out <- rbind(out, p10)
    }

    if (cross.world=="01" | cross.world=="both") {

        p01 <- dat[dat[, a.var]==0, ];  p01$sample <- "p01"

        p01$w.wt <-
            exp(predict(a.cm.fu, newdata = p01, type = "link")) *
            (1 + exp(-predict(a.c.fu, newdata = p01, type = "link")))

        p01$w.wt <- trunc_wt(p01$w.wt, max.wt$control)

        p01$f.wt <- p01$s.wt * p01$w.wt

        out <- rbind(out, p01)
    }

    # manage output
    if (!compact) {
        attr(out, "ps.type") <- "stack"
    } else {
        out <- stack_to_compact(out)
        attr(out, "ps.type") <- "compact"
    }

    out

}



#' @name get_ipw_weights
#' @description Estimates inverse probability weights that form two pseudo samples mimicking the covariate distribution of the full sample.
#' @param
#' data The sample you do analysis on
#' s.wt Sampling weights. Defaults to NULL. This argument is also used in bootstrapping (in the boot_an_estimator() function) for listing a variable that hold the bootstrap weights.
#' @return
#' @export

get_ipw_weights <- function(
    data,
    s.wt = NULL,
    a.c.form,

    max.stabilized.wt = 30

) {

    # data check: treatment coded as 0/1
    a.var <- as.character(formula(a.c.form)[[2]])

    if (!is_binary01(data[, a.var]))
        stop("Variable named in \"a.var\" must be numeric and in binary 0/1 form.")



    # working (and template) dataset
    dat <- data.frame(rownum = 1:nrow(data),
                      s.wt = NA,       # sampling weight
                      sample = NA,     # which pseudo sample
                      w.wt = NA,       # dist. morphing weights (omega in paper)
                      f.wt = NA,         # s.wt * w.wt
                      data)

    if (is.null(s.wt)) { dat$s.wt <- rep(1, nrow(data))
    } else             { dat$s.wt <- data[, s.wt]
    }



    # model
    a.c.fu <- glm(formula = a.c.form,
                  data    = dat,
                  weights = s.wt,
                  family  = quasibinomial)



    # prep for w.wt truncation: max.wt values and trunc_wt() function
    if (is.null(max.stabilized.wt)) {
        max.wt <- list(control = Inf, treat   = Inf)
    } else {
        max.wt <- lapply(list(control=0, treat=1), function(z) {
            max.stabilized.wt *
                (sum(dat$s.wt) / sum(dat$s.wt * (dat[, a.var]==z)))
        })
    }

    trunc_wt <- function(vec, max.val) {
        (vec <= max.val) * vec + (vec > max.val) * max.val
    }



    # p00 and p11 pseudo samples
    p00 <- dat[dat[, a.var]==0, ];  p00$sample <- "p00"
    p11 <- dat[dat[, a.var]==1, ];  p11$sample <- "p11"

    p00$w.wt <- 1 + exp( predict(a.c.fu, newdata = p00, type = "link"))
    p11$w.wt <- 1 + exp(-predict(a.c.fu, newdata = p11, type = "link"))

    p00$w.wt <- trunc_wt(p00$w.wt, max.wt$control)
    p11$w.wt <- trunc_wt(p11$w.wt, max.wt$treat)

    p00$f.wt <- p00$s.wt * p00$w.wt
    p11$f.wt <- p11$s.wt * p11$w.wt



    out <- rbind(p00, p11)

    attr(out, "ps.type") <- "regs"

    out

}




#' @name get_Ypred_weights
#' @description Estimates odds weights for the Ypred estimator
#' @param
#' data The sample you do analysis on
#' s.wt Sampling weights. Defaults to NULL. This argument is also used in bootstrapping (in the boot_an_estimator() function) for listing a variable that hold the bootstrap weights.
#' @return
#' @export

get_Ypred_weights <- function(
    data,
    s.wt = NULL,
    cross.world = "10",
    a.cm.form,

    max.stabilized.wt = 30
) {

    print("under development")
}





#' @name plot_balance
#' @importFrom ggplot2 ggplot aes geom_vline geom_point scale_color_manual scale_shape_manual labs theme_bw facet_wrap xlim

plot_balance <- function(ps.dat,
                         sample = "sample",
                         s.wt = "s.wt",
                         f.wt = "f.wt",
                         c.vars,
                         m.vars = NULL,
                         standardize,
                         full.label = FALSE) {

    suppressWarnings(
        tmp <- make_dummies(ps.dat,
                            columns = c.vars,
                            output.new.varnames = TRUE)
    )
    ps.dat <- tmp$data
    c.vars <- tmp$columns
    rm(tmp)

    if (!is.null(m.vars)) {
        suppressWarnings(
            tmp <- make_dummies(ps.dat,
                                columns = m.vars,
                                output.new.varnames = TRUE)
        )
        ps.dat <- tmp$data
        m.vars <- tmp$columns
        rm(tmp)
    }





    SMDs <- compute_SMDs(ps.dat = ps.dat,
                         sample = sample,
                         s.wt = s.wt,
                         f.wt = f.wt,
                         c.vars = c.vars,
                         m.vars = m.vars,
                         standardize = standardize,
                         full.label = full.label)

    plot.vars <- c(c.vars, m.vars)
    plot.vars <- plot.vars[length(plot.vars):1]

    SMDs$variable <- factor(SMDs$variable, levels = plot.vars)

    p <- ggplot(data = SMDs,
                aes(x = smd,
                    y = variable)) +
        geom_vline(xintercept = 0,
                   color = "gray60")

    if ("mediator" %in% unique(SMDs$var.type)) {
        p <- p +
            geom_point(aes(color = var.type,
                           shape = type),
                       fill = "white",
                       size = 1.5,
                       stroke = .5) +
            scale_color_manual(name = "", values = c(1, "magenta"))

    } else {
        p <- p +
            geom_point(aes(shape = type),
                       fill = "white",
                       size = 1.5,
                       stroke = .5)
    }

    p <- p +
        labs(x = "differences in means", y = "") +
        scale_shape_manual(name = "", values = c(21, 19)) +
        theme_bw() +
        facet_wrap(~ contrast, ncol = 3) +
        xlim(min(-.3, min(SMDs$smd)),
             max(.3, max(SMDs$smd)))

    p
}



