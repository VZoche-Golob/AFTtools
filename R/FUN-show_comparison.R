#' Printable Comparison of Different Distributions for Parametric Survival Models
#'
#' @param comp A list returned by \code{\link{compare_PSMdist}}.
#'
#' @examples
#' intS2 <- with(MIC, create_int2Surv(concentration, inhibition))
#' dd <- compare_PSMdist(as.formula("intS2 ~ region"), Data = cbind(intS2, MIC))
#' show_comparison(dd)
#' rm(intS2, dd)
#'
#' @importFrom magrittr '%>%'
#'
#' @export
show_comparison <- function(comp) {

  assertive::assert_is_list(comp)
  assertive::assert_all_are_true(names(comp) == c("Formula",
                                                    "Distribution",
                                                    "Fit",
                                                    "Message",
                                                  "AIC",
                                                  "Effects",
                                                  "logScale",
                                                  "FrailtyVar"))
  with(comp, {

    assertive::assert_is_all_of(Formula, "formula")
    assertive::assert_is_character(Distribution)
    assertive::assert_is_logical(Fit)
    assertive::assert_is_character(Message)
    assertive::assert_is_numeric(AIC)
    assertive::assert_is_numeric(Effects)
    assertive::assert_is_numeric(logScale)
    assertive::assert_is_numeric(FrailtyVar)
    assertive::assert_are_same_length(Distribution, Fit)
    assertive::assert_are_same_length(Distribution, Message)
    assertive::assert_are_same_length(Distribution, AIC)
    assertive::assert_are_same_length(Distribution, FrailtyVar)
    assertive::assert_is_array(Effects)
    assertive::assert_is_matrix(logScale)

  })



  tab <- with(comp, {

    sapply(seq_along(Distribution), function(x) {

      ef <- Effects[,,x] %>%
      {apply(.,
             1,
             function(xx)
               paste0(format(xx[1], digits = 1, width = 5, nsmall = 2),
                      " [",
                      format(xx[2], digits = 1, width = 5, nsmall = 2),
                      ", ",
                      format(xx[3], digits = 1, width = 5, nsmall = 2),
                      "]",
                      collapse = ""))}

      c(Distribution = Distribution[x],
        Fitted = Fit[[x]],
        AIC = format(AIC[[x]], digits = 1, width = 7, nsmall = 1),
        'Var(frailty)' = format(FrailtyVar[[x]], digits = 1, width = 5, nsmall = 2),
        Coef = ef,
        'log(scale)' = paste0(format(logScale[x, 1], digits = 1, width = 5, nsmall = 2),
                              " [",
                              format(logScale[x, 2], digits = 1, width = 5, nsmall = 2),
                              ", ",
                              format(logScale[x, 3], digits = 1, width = 5, nsmall = 2),
                              "]",
                              collapse = "")) %>%
        return

    })

  })



  out <- list(
    Formula = comp$Formula,
    DistComp = tab
  )



  return(out)

}



