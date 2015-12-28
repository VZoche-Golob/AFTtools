#' Q-Q-plot of Predicted Survival Times
#'
#' Accelerated Failure Time models are based on the assumption that different levels
#'  of a predictor (e. g. different treatments) have a multiplicative effect on the
#'  survival time. Plotting the predicted survival time quantiles of different
#'  predictor levels against each other allows to assess the validity of the assumption.
#'  If the assumption is valid, the quantile pairs will approximately lie on a straight
#'  line through the origin.
#'
#' @note If the \code{model} includes frailties, they have to be estimated with
#'  \code{sparse = FALSE} to allow predictions.
#'
#' @param model An object of class \code{\link[survival]{survreg}}.
#'
#' @references
#'  M. J. Bradburn, T. G. Clark, S. B. Love, und D. G. Altman, „Survival analysis
#'  Part III: multivariate data analysis -- choosing a model and assessing its
#'  adequacy and fit“, Br. J. Cancer, Bd. 89, Nr. 4, S. 605–611, Aug. 2003.
#' @references
#'  W. R. Swindell, „ACCELERATED FAILURE TIME MODELS PROVIDE A USEFUL STATISTICAL
#'  FRAMEWORK FOR AGING RESEARCH“, Exp Gerontol, Bd. 44, Nr. 3, S. 190–200, März 2009.
#'
#' @import survival
#' @importFrom magrittr '%>%'
#'
#' @export
AFTplot <- function(model) {

  assertive::assert_is_all_of(model, "survreg")



  par_set <- par(mar = c(4, 4, 0.5, 0.5),
                 oma = c(0, 0, 2, 0),
                 mfrow = c(2, 2))



  quantPredict <- cbind(model.frame(model),
                        Q = predict(model, type = "quantile", p = c(1:49 / 50)))



  fixEf <- names(model$xlevels) %>%
  {.[which(!grepl("frailty.", .))]}

  if (is.null(fixEf)) {

    stop("Only one fixed effect level found in model. Do not know how to produce a Q-Q-plot.")

  }

  ii <- 0
  for (i in seq_along(fixEf)) {

    vals <- lapply(model$xlevels[[i]], function(x) {

      subset(quantPredict, quantPredict[, fixEf[i]] == x) %>%
      {.[, grep("Q.", names(.))]}

    })

    cbs <- combn(length(model$xlevels[[i]]), 2)
    for (j in 1:ncol(cbs)) {

      qqplot(x = as.matrix(vals[[cbs[1, j]]]),
             y = as.matrix(vals[[cbs[2, j]]]),
             xlab = paste(fixEf[i], "==", model$xlevels[[i]][cbs[1, j]]),
             ylab = paste(fixEf[i], "==", model$xlevels[[i]][cbs[2, j]]))

       XX <- quantile(as.matrix(vals[[cbs[1, j]]]), probs = seq(0, 1, 0.1))
       YY <- quantile(as.matrix(vals[[cbs[2, j]]]), probs = seq(0, 1, 0.1))
       bb <- lm(YY ~ 0 + XX)$coefficients

       stopifnot(length(bb) == 1)

       abline(a = 0, b = bb, lty = 2)

      if (assertive::is_whole_number(ii / 4)) {

        title(main = "Q-Q-Plot of predicted survival times",
              outer = TRUE)

      }
      ii <- ii + 1

    }

  }

  par(par_set)

}