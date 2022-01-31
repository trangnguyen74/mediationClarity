



#### synth #####################################################################

#' This is data to be included in my package
#'
#' @name synth
#' @docType data
#' @author My Name \email{blahblah@@roxygen.org}
#' @references \url{data_blah.com}
#' @keywords data
NULL



#### reshape_gather ############################################################

#' A cross between \code{gather} of package \code{dplyr} and \code{reshape}
#'
#' This function is like \code{gather} and \code{reshape(direction="long")}. Like \code{gather}, it outputs names of variables that are gathered (ie turned from wide to long format). (\code{reshape} outputs a numeric variable for assumed discrete time.) In addition to variable name and value, it outputs (like \code{reshape}) a third variable retaining the wide format row information. \code{reshape} names this variable \code{id}; the current function names the variable \code{wide.row}.
#' @param data A data frame
#' @param columns Columns to be gathered
#' @param key New variable to store names of columns being gathered
#' @param value New variable to contain values of columns being gathered
#' @param wide.row Whether to include a variable named \code{wide.row} containing the row numbers from the original wide data. Defaults to TRUE.
#' @return A data frame with the two new variables named in arguments \code{key} and \code{value}, plus a third new variable named \code{wide.row}. The remaining variables of the input dataset are carried over.
#' @export
#'
reshape_gather <- function(data,
                           columns,
                           key,
                           value,
                           wide.row = TRUE) {

    if ("id" %in% names(data))
        names(data)[names(data)=="id"] <- "id.original"

    data <- stats::reshape(data,
                           varying = columns,
                           timevar = key,
                           v.names = value,
                           direction = "long")

    if (wide.row) {
        names(data)[names(data)=="id"] <- "wide.row"
    } else {
        data <- data[, names(data)!="id"]
    }

    data[, key] <- factor(data[, key], labels = columns)
    data[, key] <- as.character(data[, key])

    if ("id.original" %in% names(data))
        names(data)[names(data)=="id.original"] <- "id"

    rownames(data) <- NULL

    data
}



#### is_constant ###############################################################

#' Check if a vector is constant
#'
#' @param x Vector or matrix to be tested.
#' @return TRUE if constant, FALSE if not.
#' @noRd

is_constant <- function(x) {

    x <- c(x)

    all(x==x[1])
}




#### is_binary01 ###############################################################

#' Check if a vector is binary coded as 0/1
#'
#' @param x Vector (or matrix, array) to be tested.
#' @noRd

is_binary01 <- function(x) {
    is.numeric(x) && all(x==1 | x==0)
}




#### maybe_continuous ##########################################################

#' Check if a variable may be continuous
#'
#' @noRd

maybe_continuous <- function(x) {
    is.numeric(x) && !is_binary01(x)
}




#### .trunc_right ##############################################################

#' Right truncation
#' Right truncation
#' @param vec A numeric vector.
#' @param max A single value of a numeric vector of the same size as \code{vec} indicating the max value(s) used to right-truncate \code{vec}.
#' @return Vector of truncated values.

.trunc_right <- function(vec, max) {
    (vec <= max) * vec + (vec > max) * max
}




#### .dummy_one_variable #######################################################

# Turn a character/factor variables into dummies

.dummy_one_variable <- function(data, variable) {


    if (is.factor(data[, variable]))    vals <- levels(data[, variable])
    if (is.character(data[, variable])) vals <- sort(unique(data[, variable]))


    if (length(vals)<=2) {

        data[, variable] <- as.integer(as.character(data[, variable]))

        out <- list(data = data,
                    dummies = variable)

    } else {
        out.variables <- paste(variable, vals, sep = "_")

        dummies <- sapply(vals, function(z) {

            1 * (data[, variable]==z)
        })
        colnames(dummies) <- out.variables

        out <- list(data = cbind(data, dummies),
                    dummies = out.variables)
    }

    out
}



#### .make_dummies #############################################################


.make_dummies <- function(data,
                          columns,
                          output.names = FALSE,
                          warning = TRUE) {

    # TO DO: need better name for argument 3

    out.columns <- NULL

    for (v in 1:length(columns)) {

        if (!is.character(data[, columns[v]]) && !is.factor(data[, columns[v]])) {
            if (warning)
                warning(paste("Dummies not created for variable",
                              columns[v],
                              "(not a factor/character variable)."))
            out.columns <- c(out.columns, columns[v])
        } else {
            tmp <- .dummy_one_variable(data, columns[v])
            out.columns <- c(out.columns, tmp$dummies)
            data <- tmp$data;  rm(tmp)
        }
    }


    if (output.names==FALSE) out <- data
    if (output.names==TRUE)  out <- list(data = data,
                                         columns = out.columns)

    out
}



.dummies_2sets <- function(data, columns1, columns2) {
    tmp1 <- .make_dummies(data = data,
                          columns = columns1,
                          output.names = TRUE,
                          warning = FALSE)

    tmp2 <- .make_dummies(data = tmp1$data,
                          columns = columns2,
                          output.names = TRUE,
                          warning = FALSE)

    list(data = tmp2$data,
         columns1 = tmp1$columns,
         columns2 = tmp2$columns)
}








#### .wtd_mean & .wtd_var ######################################################


#' Weighted statistics
#'
#' @param x A numeric vector
#' @param w A vector of weights
#' @return \code{.wtd_mean()} returns weighted mean. \code{.wtd_var()} returns weighted variance.
#' @name dot-wtd_
NULL


#' @rdname dot-wtd_
.wtd_mean <- function(x, w) {
    mean(x*w, na.rm = TRUE) / mean(w, na.rm = TRUE)
}

#' @rdname dot-wtd_
.wtd_var <- function(x, w) {
    m <- .wtd_mean(x, w)
    v <- .wtd_mean((x-m)^2, w)
    v * length(x) / (length(x) - 1)
}





#### .get_seed #################################################################

#' Sample a seed or seed vector
#'
#' This function obtains a vector of integers (or a single integer) to be used as seeds. It preserves RNG state while doing so.
#' @param size Number of seeds wanted.
#' @param master Master seed for reproducibility
#' @return A vector of sampled integers (or a single integer if  \code{size=1}).

.get_seed <- function(size = 1,
                      master = NULL) {

    preRNG <- .current_RNG()

    {
        if (!is.null(master)) set.seed(master)
        seed <- sample.int(n = .Machine$integer.max,
                           size = size,
                           replace = FALSE)
    }

    .restore_RNG(preRNG);  rm(preRNG)

    seed
}






#### .get_cur_RNG & .restore_RNG ###############################################


# TODO: figure out a function (or some device) to wrap a chunk of code that may change RNG state so as to preserve RNG code.


#' A pair of functions to protect RNG state
#'
#' Includes \code{.current_RNG()} which obtains current RNG state, and \code{.restore_RNG()} which restores to a previous RNG state
#' @param RNGstate A RNG state, e.g., the output of a previous \code{get_cur_RNG()} call.

.current_RNG <- function() {

    if (exists(".Random.seed", .GlobalEnv)) {
        RNGstate <- .GlobalEnv$.Random.seed
    } else                                  {
        RNGstate <- NULL
    }
    RNGstate
}

#' @describeIn  dot-current_RNG
.restore_RNG <- function(RNGstate) {

    if (!is.null(RNGstate)) { .GlobalEnv$.Random.seed <- RNGstate
    } else                  { rm(".Random.seed", envir = .GlobalEnv)
    }
}




#### .get_sd.pooled ############################################################

#' Compute pooled SD
#'
#' Compute pooled SD for a variable.
#' @param dat1,dat0 Treated and control subsamples.
#' @param variable Name of the variable (covariate or mediator) for which to compute pooled SD.
#' @return Pooled SD.

.get_sd.pooled <- function(variable, dat1, dat0) {

    size1 <- sum(dat1$.s.wt)
    size0 <- sum(dat0$.s.wt)

    s2.1 <- .wtd_var(dat1[, variable], dat1$.s.wt)
    s2.0 <- .wtd_var(dat0[, variable], dat0$.s.wt)

    s2.pooled <- (s2.1 * size1 + s2.0 * size0) / (size1 + size0)

    sqrt(s2.pooled)
}














