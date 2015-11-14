#' Create an interval censored Survival object
#'
#' Create a survival object from interval censored data with fixed observation
#'  times where only times with events or the last observation times (leading
#'  to right censoring) were reported. Observation times were defined before
#'  conducting the study and that were the same for all subjects.
#'
#' @param obsTimes A numeric vector with the observed survival times (the right
#'  interval bounds.)
#' @param ev Were events observed at \code{obsTimes}? A vector of the same length
#'  as \code{obsTimes}, either \code{logical} oder \code{integer} with 0 = censored
#'  and 1 = event.
#' @param defTimes A numeric vector with all observation times.
#'
#' @return An object of class \code{\link[survival]{Surv}} with \code{type =
#'  "interval2"}.
#'
#' @examples with(MIC, create_int2Surv(concentration, inhibition))
#'
#' @export
create_int2Surv <- function(obsTimes, ev, defTimes = unique(obsTimes)) {

  assertive::assert_is_numeric(defTimes)
  assertive::assert_is_vector(defTimes)
  assertive::assert_has_no_duplicates(defTimes)
  assertive::assert_all_are_not_na(defTimes)
  assertive::assert_is_numeric(obsTimes)
  assertive::assert_is_vector(obsTimes)
  assertive::assert_all_are_true(obsTimes %in% defTimes)
  assertive::assert_all_are_not_na(ev)
  assertive::assert_are_same_length(obsTimes, ev)

  if (!is.logical(ev)) {

    assertive::assert_all_are_true(ev %in% c(0, 1))
    ev <- as.logical(ev)

  }



  Time2 <- obsTimes
  Time2[!ev] <- NA

  Time <- sapply(obsTimes,
                 function(x, dt) {

                   lv <- sort(dt)

                   out <- which(lv == x) - 1
                   out[out == 0] <- 1
                   out <- lv[out]

                   return(out)

                 }, dt = defTimes)
  Time[Time2 == min(Time2, na.rm = TRUE)] <- NA

  return(survival::Surv(time = Time,
                        time2 = Time2,
                        type = "interval2"))

}


