#' Cox-Snell-Residual Plot
#'
#' Plots the Cox-Snell-Residuals of a model fitted using \code{\link[survival]{survreg}}
#'  against the Kaplan-Meier estimator of their cumulative hazard.
#'
#' @note If the \code{model} includes frailties, they have to be estimated with
#'  \code{sparse = FALSE} to allow the residuals to be extracted.
#'
#' @param model An object of class \code{\link[survival]{survreg}}.
#'
#' @references
#'  P. Glomb, „Statistische Modelle und Methoden in der Analyse von Lebenszeitdaten“,
#'  Diplomarbeit, Carl-von-Ossietzky-Universität, Oldenburg, 2007.
#'
#' @import survival
#' @importFrom magrittr '%>%'
#'
#' @export
CSRplot <- function(model) {

  assertive::assert_is_all_of(model, "survreg")



  # extract Cox-Snell-Residuals
  CSresiduals <- function(mod) {

    CS_transform <- function(x, distr) {

      switch(distr,
             weibull = exp(x),
             exponential = exp(x),
             lognormal = -log(1 - pnorm(x, mean = 0, sd = 1)),
             loglogistic = log(1 + exp(x)),
             stop("Unknown model distribution."))

    }

    res <- residuals(mod, "working") %>%
      CS_transform(., mod$dist) %>%
      return

  }
  res <- CSresiduals(model) %>%
  {survfit(Surv(., model$y[, ncol(model$y)] > 0) ~ 1)}



  plot(res$time, -log(res$surv), type = "b",
       lty = 1, pch = 20, col = "red",
       xlab = "Cox-Snell-Residuals", ylab = "Cumulative hazard")
  abline(0, 1, lty = 2)

}