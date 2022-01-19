#' Estimate natural (in)direct effects using an estimator of choice
#'
#' blah blah
#' @param
#' data The sample you do analysis on
#' s.wt Name of the sampling weights. Defaults to NULL. This argument is also used in bootstrapping (in the boot_an_estimator() function) for listing a variable that hold the bootstrap weights.
#' @return
#' @export

estimate_effects <- function(
    data,
    s.wt = NULL,

    estimator = "Y2pred",
    cross.world = "10",
    effect.scale = "additive",

    boot.stratify = TRUE,
    boot.num = 999,
    boot.seed = NULL,
    boot.method = "dirichlet",

    ...

) {

    if (is.null(s.wt)) { data$s.wt <- 1
    } else             { data$s.wt <- data[, s.wt];  s.wt <- "s.wt"
    }

    estFUN <- get_estFUN(estimator = estimator)

    arg.list <- clean_args(estimator = estimator,
                           cross.world = cross.world,
                           ...)

    arg.list[["effect.scale"]] <- clean_scale(effect.scale)


    # point estimation
    point.args <- arg.list
    point.args[["data"]] <- data

    point <- do.call(pasteo("point_", estimator), point.args)


    # bootstrap for CI and SE
    if (boot.stratify) { boot.stratify <- arg.list$a.var
    } else             { boot.stratify <- NULL
    }

    if (is.null(boot.seed)) boot.seed <- get_seed()


    ci.se <- boot_ci.se(data = data,
                        boot.stratify = boot.stratify,
                        boot.num = boot.num,
                        boot.seed = boot.seed,
                        boot.method = boot.method,
                        FUN = pasteo("point_", estimator),
                        arg.list = arg.list)



    list(effects = cbind(estimate = point, ci.se),
         boot.seed = boot.seed)
}

