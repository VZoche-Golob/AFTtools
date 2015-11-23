#' Influential Plot
#'
#' Plot of the relative parameter change due to addition of each observation. In
#'  fact, it nicely plots \code{\link[survival]{residuals.survreg}(model, type = "dfbeta")}
#'  for each fixed effect separately.
#'
#' @note If the \code{model} includes frailties, they have to be estimated with
#'  \code{sparse = FALSE} to allow the residuals to be extracted.
#'
#' @param model An object of class \code{\link[survival]{survreg}}.
#' @param limit Relative parameter changes above \code{limit}\% are marked in red.
#'
#' @references
#'  W. R. Swindell, „ACCELERATED FAILURE TIME MODELS PROVIDE A USEFUL STATISTICAL
#'  FRAMEWORK FOR AGING RESEARCH“, Exp Gerontol, Bd. 44, Nr. 3, S. 190–200, März 2009.
#'
#' @import survival
#' @importFrom magrittr '%>%'
#'
#' @export
Iplot <- function(model, limit = 10) {

  assertive::assert_is_all_of(model, "survreg")
  assertive::assert_is_scalar(limit)
  assertive::assert_is_a_number(limit)
  assertive::assert_all_are_in_range(limit, 0, 100)



  fixedCoef <- dimnames(model.matrix(model))[[2]] %>%
  {.[which(!grepl("frailty.", .))]}

  changeRes <- residuals(model, type = "dfbeta")[, fixedCoef] %>%
  {if (is.null(dim(.))) {
    matrix(. * 100 / coef(model)[fixedCoef], nrow = 1)
  } else {
    apply(., 1, function(x) x * 100 / coef(model)[fixedCoef])
  }}

  fixedCoef <- dimnames(model.matrix(model))[[2]] %>%
  {.[which(!grepl("frailty.", .))]}


  if (length(fixedCoef) > 1) {

    par_set <- par(oma = c(4, 4, 1, 1),
                   mar = c(0, 0, 1.5, 0),
                   mfrow = c(2, 2))



    for (i in seq_along(fixedCoef)) {

      plot.new()
      plot.window(xlim = c(0, ncol(changeRes)),
                  ylim = range(changeRes),
                  new = FALSE)
      box()
      if (assertive::is_equal_to(i / 3, round(i / 3)) |
          assertive::is_equal_to(i / 4, round(i / 4)) |
          (assertive::is_even(i) & length(fixedCoef) - i == 1)) axis(1)
      if (assertive::is_odd(i)) axis(2)

      points(changeRes[i, ])
      lines(changeRes[i, ])
      if (any(abs(changeRes[i, ]) > limit)) {

        which(abs(changeRes[i, ]) > limit) %>%
          points(x = .,
                 y = changeRes[i, .],
                 col = "red")

        which(abs(changeRes[i, ]) > limit) %>%
          text(x = .,
               y = changeRes[i, .],
               labels = .,
               pos = 4,
               cex = 0.6)

      }

      title(main = fixedCoef[i])

    }

    title(xlab = "Observation",
          ylab = "Relative parameter change due to addition of an observation (%)",
          outer = TRUE)



    par(par_set)

  } else {

    plot.new()
    plot.window(xlim = c(0, ncol(changeRes)),
                ylim = range(changeRes),
                new = FALSE)
    box()
    axis(1)
    axis(2)

    points(changeRes[1, ])
    lines(changeRes[1, ])
    if (any(abs(changeRes[1, ]) > limit)) {

      which(abs(changeRes[1, ]) > limit) %>%
        points(x = .,
               y = changeRes[1, .],
               col = "red")

      which(abs(changeRes[1, ]) > limit) %>%
        text(x = .,
             y = changeRes[1, .],
             labels = .,
             pos = 4,
             cex = 0.6)

    }

    title(xlab = "Observation",
          ylab = "Relative parameter change due to addition of an observation (%)")

  }

}