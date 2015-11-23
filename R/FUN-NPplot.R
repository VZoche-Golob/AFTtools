#' Graphic Comparison of Parametric Survival Models to Correponding Non-Parametric Models
#'
#' Plots of the predicted survival of \code{model} and of corresponding non-parametric
#'  models.
#'
#' @details
#'  \code{\link{predict_survreg}} is used internally to calculate the predicted
#'  survival for \code{model}. As non-parametric references, the Kaplan-Meier (for
#'  right censored data) and the Turnbull (for interval censored data) models are
#'  used by calling \code{\link[survival]{survfit}}.
#'  Frailties are ignored for the non-parametric predictions.
#'
#' @note Currently only right- and interval-censored data have been tested.
#'
#' @param model An object of class \code{\link[survival]{survreg}}.
#' @param data The data used to fit \code{model}.
#'
#' @return One plot per predictor level (if more than one,
#'  \code{\link[graphics]{par}(mfrow = c(2, 2))} is used.)
#'
#' @examples
#' \donttest{
#' intS2 <- with(MIC, create_int2Surv(concentration, inhibition))
#' psm1 <- survival::survreg(as.formula("intS2 ~ region +
#'  frailty(herd, sparse = FALSE)"), data = cbind(intS2, MIC))
#' NPplot(psm1, cbind(intS2, MIC))
#'
#' psm2 <- survival::survreg(as.formula("Surv(concentration, inhibition) ~ region +
#'  frailty(herd, sparse = FALSE)"), data = MIC)
#' NPplot(psm2, MIC)
#'
#' rm(psm1, psm2, intS2)
#' }
#'
#' @references
#'  K. Goethals, B. Ampe, D. Berkvens, H. Laevens, P. Janssen, und L. Duchateau,
#'  „Modeling interval-censored, clustered cow udder quarter infection times
#'  through the shared gamma frailty model“, JABES, Bd. 14, Nr. 1, S. 1–14, März 2009.
#' @references
#'  W. R. Swindell, „ACCELERATED FAILURE TIME MODELS PROVIDE A USEFUL STATISTICAL
#'  FRAMEWORK FOR AGING RESEARCH“, Exp Gerontol, Bd. 44, Nr. 3, S. 190–200, März 2009.
#'
#' @import interval
#' @import survival
#' @importFrom magrittr '%>%'
#'
#' @export
NPplot <- function(model,
                   data) {

  assertive::assert_is_all_of(model, "survreg")



  omit_frailties <- function(x) {

    assertive::assert_is_all_of(x, "formula")

    del <- attr(x,"term.labels")

    if (is.null(del) | length(del) == 0) {

      return(update(x, x))

    } else {

      del <- del[which(grepl("frailty(.+)", del))]

      eval(parse(text = paste("return(update(x, . ~ . - ", del, "))", sep = "")))

    }
  }

  # adapted from interval::plot.icfit()
  polygon.icfit <- function(xx) {
    S <- c(1, 1 - cumsum(xx$pf))
    S <- rep(S, each = 2)[-2 * length(S)]
    time <- c(0, as.vector(xx$intmap))
    if (any(time == Inf)) {
      maxtime <- max(time[time < Inf])
      time[time == Inf] <- maxtime
      Satmax <- S[time == maxtime][1]
      polygon(c(maxtime, maxtime, 2 * maxtime, 2 *
                  maxtime, maxtime), c(Satmax, 0, 0, Satmax,
                                       Satmax), col = "lightgrey", border = NA)
    }
    tt <- rep(time, each = 2)
    tt <- c(tt[-1], tt[(length(tt) - 1):1])
    SS <- rep(S, each = 2)
    SS <- c(SS[-length(SS)], SS[length(SS):2])
    polygon(tt, SS, col = "lightgrey", border = NA)
  }



  # type of censoring
  type <- attr(with(data, eval(terms(model)[[2]])),
               "type")



  predictors <- names(model$xlevels) %>%
  {.[which(!grepl("frailty.", .))]} %>%
  {model$xlevels[.]}

  fm <- omit_frailties(formula(model))



  # non-parametric
  npL <- survfit(fm, data = data, se.fit = FALSE)

  # x-values
  xv <- seq(1e-10,
            max(npL$time[which(is.finite(npL$time))]),
            length.out = 100)

  if (type == "interval") {

    gr <- icfit(fm, data = data, conf.int = FALSE)

    # adaption of x-values
    xv <- seq(1e-10,
              max(xv, max(gr$intmap[which(is.finite(gr$intmap))])),
              length.out = 100)

  }



  # model predictions
  modpred <- predict_survreg(model = model,
                             data = data,
                             type = "survival",
                             times = xv,
                             conf.int = NULL)



  if (length(predictors) > 0) {

    dfpred <- expand.grid(predictors, KEEP.OUT.ATTRS = FALSE)

    par_set <- par(oma = c(4, 4, 3, 1),
                   mar = c(0, 0, 1.5, 0),
                   mfrow = c(2, 2))



    for (i in seq_along(dfpred[, 1])) {

      plot.new()
      plot.window(xlim = c(0, max(xv)),
                  ylim = c(0, 1),
                  new = FALSE)
      box()
      if (assertive::is_equal_to(i / 3, round(i / 3)) |
          assertive::is_equal_to(i / 4, round(i / 4)) |
          (assertive::is_even(i) & nrow(dfpred) - i == 1)) axis(1)
      if (assertive::is_odd(i)) axis(2)
      if (type == "interval") {

        polygon.icfit(gr[i])

      }
      lines(npL[i], conf.int = FALSE, mark.time = FALSE, lty = 2)
      lines(xv, modpred[[i]][[2]])

      title(main = paste0(dfpred[i, ], collapse = ", "))

      if (i == 1 | assertive::is_equal_to((i - 1) / 4, round((i - 1) / 4))) {

        legend("topright",
               legend = c("parametric",
                          "non-parametric"),
               lty = c(1, 2))

      }

    }

    title(main = "Comparison of parametric and non-parametric predictions",
          xlab = "t",
          ylab = expression(hat(S) (t)),
          outer = TRUE)

    par(par_set)

  } else {

    plot.new()
    plot.window(xlim = c(0, max(xv)),
                ylim = c(0, 1),
                new = FALSE)
    box()
    axis(1)
    axis(2)
    if (type == "interval") {

      polygon.icfit(gr)

    }
    lines(npL, conf.int = FALSE, mark.time = FALSE, lty = 2)
    lines(xv, modpred[[1]][[2]])

    title(main = "Comparison of parametric and\nnon-parametric predictions",
          xlab = "t",
          ylab = expression(hat(S) (t)),
          outer = FALSE)

      legend("topright",
             legend = c("parametric",
                        "non-parametric"),
             lty = c(1, 2))

  }

}
