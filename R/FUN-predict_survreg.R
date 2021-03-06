#' Predictions from Parametric Survival Models
#'
#' Predicts survival or quantiles from Parametric Survival Models that were fitted
#'  using \code{\link[survival]{survreg}}. Unlike with \code{\link[survival]{predict.survreg.penal}},
#'  predictions are population averages (i.e. the frailty effect is fixed at its
#'  mean). The confidence bounds of the predictions are determined using nonparametric
#'  bootstrap resampling.
#'
#' @param model An object of class \code{\link[survival]{survreg}}.
#' @param data The data used to fit \code{model}.
#' @param strata The name of the variable in \code{data} that should be used as
#'  \code{strata} by \code{\link[boot]{boot}}. (It seems to be a good idea to use
#'  the frailty variable as \code{strata}.)
#' @param type A character string: one of \code{c("quantile", "survival")}.
#' @param quantiles A numeric vector of survival quantiles for which the times
#'  are predicted.
#' @param times A numeric vector of survival times for which the survival percentage
#'  should be predicted.
#' @param conf.int \code{NULL} or numeric vector of length one giving the desired
#'  confidence level.
#' @param R The number of bootstrap replicates. Passed to \code{\link[boot]{boot}}.
#' @param parallel Passed to \code{\link[boot]{boot}}.
#' @param ncpus Passed to \code{\link[boot]{boot}}.
#' @param ... Further arguments passed to \code{\link[boot]{boot}}.
#'
#' @note Predictions are returned for all levels of all predictors (quantitative
#'  predictors not tested yet).
#' @note Unlike \code{\link[boot]{boot}}, if \code{ncpus = NULL} and \code{parallel != "no"}
#'  all but one or two CPU are used on UNIX-Systems and on other platforms, respectively.
#'
#' @return A list with one element per predictor level combination that contains a
#'  data frame with one row and columns for the predictors, and a one or three column
#'  matrix (point prediction, lower and upper confidence limit). For each desired
#'  time / survival quantile, one row in the matrix is returned.
#'
#' @examples
#' \donttest{
#' intS2 <- with(MIC, create_int2Surv(concentration, inhibition))
#' psm <- survival::survreg(as.formula("intS2 ~ region +
#'  frailty(herd, sparse = FALSE)"), data = cbind(intS2, MIC))
#' predict_survreg(psm, data = cbind(intS2, MIC), strata = "herd",
#'  conf.int = 0.95, parallel = "snow")
#' predict_survreg(psm, data = cbind(intS2, MIC), strata = "herd",
#'  type = "survival", times = c(0.5, 1), conf.int = 0.95, parallel = "snow")
#' rm(psm, intS2)
#' }
#'
#' @import survival
#' @importFrom magrittr '%>%'
#'
#' @export
predict_survreg <- function(model,
                            data,
                            strata = NULL,
                           type = "quantile",
                           quantiles = c(0.5, 0.1),
                           times = NULL,
                           conf.int = NULL,
                           R = 500,
                           parallel = "no",
                           ncpus = NULL,
                           ...) {

  assertive::assert_is_all_of(model, "survreg")
  if (!is.null(strata)) {

    assertive::assert_is_character(strata)
    assertive::assert_is_of_length(strata, 1)
    assertive::assert_all_are_true(strata %in% names(data))

  } else {

    strata <- rep(1, nrow(data))

  }
  assertive::assert_all_are_true(type %in% c("quantile", "survival"))
  assertive::assert_is_numeric(quantiles)
  assertive::assert_is_vector(quantiles)
  assertive::assert_all_are_in_range(quantiles, 0, 1,
                                     lower_is_strict = TRUE, upper_is_strict = TRUE)
  if (!is.null(times)) {

    assertive::assert_is_numeric(times)
    assertive::assert_is_vector(times)
    assertive::assert_all_are_positive(times)

  }
  if (!is.null(conf.int)) {

    assertive::assert_is_numeric(conf.int)
    assertive::assert_is_of_length(conf.int, 1)
    assertive::assert_all_are_in_range(conf.int, 0, 1,
                                       lower_is_strict = TRUE, upper_is_strict = TRUE)

  }

  if (is.null(ncpus) & parallel != "no") {

    if (assertive::is_unix()) {

      ncpus <- system2(command = "lscpu", args = "-p", stdout = TRUE) %>%
      {length(.) - 4} %>%
        -1

    } else {

      message("'cpus' was not specified. 'cpus = 2' is used\n")
      ncpus <- 2

    }

  }



  extract_model <- function(mod) {

    predictors <- names(mod$xlevels) %>%
    {.[which(!grepl("frailty.", .))]} %>%
    {mod$xlevels[.]}

    if (length(predictors) == 0) {

      ndat <- data.frame(linPred = 0, V1 = 1)
      mm <- model.matrix(as.formula("linPred ~ 1"), ndat)

    } else {

      ndat <- cbind(linPred = 0, expand.grid(predictors))
      mm <- model.matrix(as.formula(paste("linPred",
                                          paste0(names(predictors), collapse = " + "),
                                          sep = " ~ ")),
                         ndat)

    }

    linPred <- mm %*% coef(mod)[dimnames(mm)[[2]]]

    b <- mod$scale

    return(list(ndat = ndat, linPred = linPred, b = b))

  }

  extr <- extract_model(mod = model)



  vals <- switch(type,
                 quantile = quantiles,
                 survival = times)

  selection <- expand.grid(type = c("quantile", "survival"),
                           distr = c("weibull", "exponential", "lognormal", "loglogistic")) %>%
                           {which(.[, "type"] == type & .[, "distr"] == model$dist)}
  do_predict <- switch(selection,
                       function(lp, B, pv) {
                         exp(B * log(-log(pv)) + lp)
                       },  # quantile     weibull
                       function(lp, B, pv) {
                         exp(-exp((log(pv) - lp) / B))
                       },  # survival     weibull
                       function(lp, B, pv ) {
                         exp(log(-log(pv)) + lp)
                       },  # quantile exponential
                       function(lp, B, pv) {
                         exp(-exp((log(pv) - lp)))
                       },  # survival exponential
                       function(lp, B, pv) {
                         exp(B * qnorm(1 - pv, mean = 0, sd = 1) + lp)
                       },  # quantile   lognormal
                       function(lp, B, pv) {
                         1 - pnorm((log(pv) - lp) / B, mean = 0, sd = 1)
                       },  # survival   lognormal
                       function(lp, B, pv) {
                         exp(B * log(1 / pv - 1) + lp)
                       },  # quantile loglogistic
                       function(lp, B, pv) {
                         1 / (1 + exp((log(pv) - lp) / B))
                       })  # survival loglogistic



  # function to calculate bootstrap 'statistics'
  bfun <- function(dd, id, iLP, fm = formula(model), ds = model$dist, bv = vals) {

    # packages have to be loaded for all clusters
    require(survival)
    require(magrittr)

    bm <- try(survreg(fm, data = dd[id, ], dist = ds), silent = TRUE)

    if ("try-error" %in% class(bm)) {

      return(NA)

    } else {

      bex <- extract_model(bm)
      do_predict(bex[["linPred"]][iLP], B = bex[["b"]], pv = bv) %>%
        return

    }

  }



  ppred <- lapply(extr[["linPred"]], do_predict, B = extr[["b"]], pv = vals)  # point predictions

  if (is.null(conf.int)) {

    outnames <- switch(type,
                       quantile = list(quantile = vals, ptime = "point"),
                       survival = list(time = vals, psurvival = "point"))

    out <- lapply(seq_along(ppred), function(x) {

      X <- data.frame(extr[["ndat"]][x, -1], stringsAsFactors = FALSE)
      names(X) <- names(extr[["ndat"]])[-1]
      list(X = X,
           predictions = matrix(ppred[[x]], ncol = 1,
                                dimnames = outnames)) %>%
        return

    })

  } else {

    outnames <- switch(type,
                       quantile = list(quantile = vals,
                                       ptime = c("point",
                                                 paste(format(c(0.5 - conf.int / 2,
                                                                0.5 + conf.int / 2) * 100,
                                                              trim = TRUE,
                                                              nsmall = 1),
                                                       "%"))),
                       survival = list(time = vals,
                                       psurvival = c("point",
                                                     paste(format(c(0.5 - conf.int / 2,
                                                                    0.5 + conf.int / 2) * 100,
                                                                  trim = TRUE,
                                                                  nsmall = 1),
                                                           "%"))))



    out <- lapply(seq_along(ppred), function(x) {

      bint <- boot::boot(data = data,
                         statistic = bfun,
                         R = R,
                         strata = data[, strata],
                         iLP = x,
                         parallel = parallel, ncpus = ncpus)

      bint <- sapply(seq_along(vals), function(xx) {

        quantile(bint$t[, xx],
                 probs = c(0.5 - conf.int / 2,
                           0.5 + conf.int / 2),
                 na.rm = TRUE)

      }) %>%
      {matrix(., ncol = 2, byrow = TRUE)}



      X <- data.frame(extr[["ndat"]][x, -1], stringsAsFactors = FALSE)
      names(X) <- names(extr[["ndat"]])[-1]



      bout <- list(X = X,
                   predictions = cbind(ppred[[x]], bint))
      dimnames(bout[["predictions"]]) <- outnames

      return(bout)

    })

  }



  return(out)

}



