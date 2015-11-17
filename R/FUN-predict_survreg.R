#' Predictions from Parametric Survival Models
#'
#' Predicts survival or quantiles from Parametric Survival Models that were fitted
#'  using \code{\link[survival]{survreg}}. Unlike with \code{\link[survival]{predict.survreg.penal}},
#'  predictions are population averages (i.e. the frailty effect is set to its
#'  \code{mean = 0}).
#'
#' @param mod An object of class \code{\link[survival]{survreg}}.
#' @param type A character string: one of \code{c("quantile", "survival")}.
#' @param quantiles A numeric vector of survival quantiles for which the times
#'  are predicted.
#' @param times A numeric vector of survival times for which the survival percentage
#'  should be predicted.
#' @param conf.int Logical. Shall the confidence intervals be predicted as well?
#'  (currently not supported)
#'
#' @note Predictions are returned for all levels of all predictors (quantitative
#'  predictors not tested yet).
#'
#' @return For \code{type = survival} a data frame with one row per predictor level
#'  and one column per requested time.
#' @return For \code{type = quantile} a data frame with one row per predictor level
#'  and one column per requested quantile.
#'
#' @examples
#' intS2 <- with(MIC, create_int2Surv(concentration, inhibition))
#' psm <- survival::survreg(as.formula("intS2 ~ region +
#'  frailty(herd, sparse = FALSE)"), data = cbind(intS2, MIC))
#' predict_survreg(psm)
#' predict_survreg(psm, type = "survival", times = c(0.5, 1))
#' rm(psm, intS2)
#'
#' @import survival
#' @importFrom magrittr '%>%'
#'
#' @export
predict_survreg <- function(mod,
                           type = "quantile",
                           quantiles = c(0.5, 0.9),
                           times = NULL,
                           conf.int = FALSE) {

  # TODO: include conf.int = TRUE (using bootstrap ?)

  assertive::assert_is_all_of(mod, "survreg")
  assertive::assert_all_are_true(type %in% c("quantile", "survival"))
  assertive::assert_is_numeric(quantiles)
  assertive::assert_is_vector(quantiles)
  assertive::assert_all_are_in_range(quantiles, 0, 1,
                                     lower_is_strict = TRUE, upper_is_strict = TRUE)
  if (!is.null(times)) {

    assertive::assert_is_numeric(times)
    assertive::assert_is_vector(times)

  }
  assertive::assert_is_logical(conf.int)
  assertive::assert_is_scalar(conf.int)



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

  selection <- expand.grid(type = c("quantile", "survival"),
                           distr = c("weibull", "exponential", "lognormal", "loglogistic")) %>%
                           {which(.[, "type"] == type & .[, "distr"] == mod$dist)}

  do_predict <- switch(selection,
                       function(qs, B = b, lp = linPred ) {
                         exp(B * log(-log(qs)) + lp)
                       },  # quantile     weibull
                       function(ts, B = b, lp = linPred) {
                         exp(-exp((log(ts) - lp) / B))
                       },  # survival     weibull
                       function(qs, B = b, lp = linPred ) {
                         exp(log(-log(qs)) + lp)
                       },  # quantile exponential
                       function(ts, B = b, lp = linPred) {
                         exp(-exp((log(ts) - lp)))
                       },  # survival exponential
                       function(qs, B = b, lp = linPred ) {
                         exp(B * qnorm(1 - qs, mean = 0, sd = 1) + lp)
                       },  # quantile   lognormal
                       function(ts, B = b, lp = linPred) {
                         1 - pnorm((log(ts) - lp) / B, mean = 0, sd = 1)
                       },  # survival   lognormal
                       function(qs, B = b, lp = linPred ) {
                         exp(B * log(1 / qs - 1) + lp)
                       },  # quantile loglogistic
                       function(ts, B = b, lp = linPred) {
                         1 / (1 + exp((log(ts) - lp) / B))
                       })  # survival loglogistic

  if (type == "quantile") {

    if (length(predictors) == 0) {

      out <- sapply(quantiles, do_predict) %>%
      {matrix(., ncol = length(.))} %>%
        data.frame
      names(out) <- paste("Q", quantiles, sep = "")

    } else {

      out <- data.frame(ndat[, -1], sapply(quantiles, do_predict))
      names(out) <- c(names(ndat)[-1], paste("Q", quantiles, sep = ""))

    }

  }

  if (type == "survival") {

    if (length(predictors) == 0) {

      out <- sapply(times, do_predict) %>%
      {matrix(., ncol = length(.))} %>%
        data.frame
      names(out) <- paste("T", times, sep = "")

    } else {

      out <- data.frame(ndat[, -1], sapply(times, do_predict))
      names(out) <- c(names(ndat)[-1], paste("T", times, sep = ""))

    }

  }

  return(out)

}



