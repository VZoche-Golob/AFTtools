#' Compare Different Distributions for Parametric Survival Models
#'
#' Choosing the best distribution for parametric survivial models is not trivial.
#' One possibility is to select of all plausible distributions the one with which
#' the resulting models' AIC is minimal. Therefore, several models have to be
#' fitted and compared. \code{compare_PSMdist} fits the models for all distributions
#' that come into consideration by calling \code{\link[survival]{survreg}}. The
#' fitting of the models is parallelised via \pkg{\link{snowfall}}.
#'
#' @param Formula A formula expression as for other regression models passed to
#'  \code{\link[survival]{survreg}}. The response has to be a survival object
#'  as returned by the \code{\link[survival]{Surv}} function.
#' @param Data A data frame where the variables in formula are found.
#' @param Distributions A character vector. The names of the distributions to be
#'  used. See \code{\link[survival]{survreg}} for the possible distributions.
#' @param parallel Logical determinating parallel or sequential execution, passed
#'  to \code{\link[snowfall]{sfInit}}.
#' @param cpus Passed to \code{\link[snowfall]{sfInit}}.
#' @param type Passed to \code{\link[snowfall]{sfInit}}.
#' @param level Passed to \code{\link[snowfall]{sfInit}}.
#' @param ... Further arguments passed to \code{\link[snowfall]{sfInit}},
#'  \code{\link{confint}}, and \code{\link[survival]{survreg}}.
#'
#' @details If \code{parallel = TRUE} (the default) and \code{cpus}
#'  is not specified, on Unix systems all but one CPU will be used by
#'  \code{\link[snowfall]{sfInit}}. Other systems will use \code{cpus = 2}.
#'
#' @return A list consiting of:
#' \describe{
#'  \item{\code{Formula}}{The formula used for the fitting call.}
#'  \item{\code{Distribution}}{A character vector naming the distributions used
#'    for the models. It defines the order for the following list elements.}
#'  \item{\code{Fit}}{A logical vector: if \code{FALSE} the corresponding model
#'    could not be fitted without problems and a \code{Message} is returned.}
#'  \item{\code{Message}}{A character vector: \code{NA} if the corresponding model
#'    was fitted without problems, otherwise the error message.}
#'  \item{\code{AIC}}{A numeric vector with the \code{\link[stats]{AIC}} of the
#'    models. \code{NA} if the fitting of the corresponding model caused problems.}
#'  \item{\code{Effects}}{A numeric array: the levels of the model effects define
#'    the rows, three columns (estimated coefficient, lower and upper confidence
#'    limits at 95\%-level as default), and separate layers for each distribution.
#'    \code{NA} if the fitting of the corresponding model caused problems.}
#'  \item{\code{logScale}}{A numeric matrix with one row for each distribution,
#'    and three columns (estimated coefficient, lower and upper confidence limits
#'    at 95\%-level as default). \code{NA} if the fitting of the corresponding
#'    model caused problems.}
#'  \item{\code{FrailtyVar}}{A numeric vector with the estimated variance of the
#'    frailty effect (if one is included.). \code{NA} if the fitting of the
#'    corresponding model caused problems.}
#' }
#'
#' @examples
#' \donttest{
#' intS2 <- with(MIC, create_int2Surv(concentration, inhibition))
#' compare_PSMdist(as.formula("intS2 ~ region"), Data = cbind(intS2, MIC))
#' compare_PSMdist(as.formula("intS2 ~ region + frailty(herd, sparse = FALSE)"),
#'  Data = cbind(intS2, MIC), cpus = 2, control = survreg.control(maxiter = 100))
#' compare_PSMdist(as.formula("intS2 ~ 1 + frailty(herd, sparse = FALSE)"),
#'  Data = cbind(intS2, MIC), cpus = 2, control = survreg.control(maxiter = 100))
#' }
#'
#' @seealso \code{\link{show_comparison}} produces a print- and readable table of
#'  the comparison.
#'
#' @import survival
#' @import snow
#' @import snowfall
#' @importFrom magrittr '%>%'
#'
#' @export
compare_PSMdist <- function(Formula, Data,
                    Distributions = c("weibull", "exponential", "lognormal", "loglogistic"),
                    parallel = TRUE, cpus = NULL, type = "SOCK",  level = 0.95,
                    ...) {

  assertive::assert_is_all_of(Formula, "formula")
  assertive::assert_is_data.frame(Data)
  assertive::assert_is_character(Distributions)
  assertive::assert_is_a_bool(parallel)



  if (is.null(cpus) & parallel == TRUE) {

    if (assertive::is_unix()) {

      cpus <- system2(command = "lscpu", args = "-p", stdout = TRUE) %>%
      {length(.) - 4} %>%
        -1

    } else {

      message("'cpus' was not specified. 'cpus = 2' is used\n")
      cpus <- 2

    }

  }



  dots <- list(...)
  sfInit_dots <- which(names(dots) %in% c("socketHosts",
                                               "restore",
                                               "slaveOutfile",
                                               "nostart",
                                               "useRscript")) %>%
                                               {if (length(.) > 0) {
      dots[.]
                                               } else{
                                                   list()
                                               }}



  ParWrapper <- function(distr, form = Formula, data = Data, ...) {

    optOld <- getOption("warn")
    options(warn = 2)
    on.exit(options(warn = optOld))



    lvn <- c(0.5 - level / 2, 0.5 + level / 2)

    psm <- dimnames(model.matrix(form, data = data))[[2]] %>%
    {.[which(!grepl("frailty.", .))]} %>%
    {list(
      dist = distr,
      fit = as.logical(NA),
      mes = as.character(NA),
      aic = as.numeric(NA),
      effects = matrix(as.numeric(NA), ncol = 3, nrow = length(.),
                       dimnames = list(.,
                                       c("coef",
                                         paste(format(c(0.5 - 0.95 / 2, 0.5 + 0.95 / 2) * 100,
                                                      trim = TRUE,
                                                      nsmall = 1),
                                               "%")))),
      logscale = matrix(as.numeric(NA), ncol = 3,
                        dimnames = list(distr,
                                        c("logScale",
                                          paste(format(c(0.5 - 0.95 / 2, 0.5 + 0.95 / 2) * 100,
                                                       trim = TRUE,
                                                       nsmall = 1),
                                                "%")))),
      frailtyvar = as.numeric(NA)
    )}



    mod <- try(survreg(form, data = data, dist = distr, ...), silent = TRUE)

    if ("try-error" %in% class(mod)) {

      psm[["fit"]] <- FALSE
      mod <- as.character(mod)
      attributes(mod) <- NULL
      psm[["mes"]] <- paste("On fitting:", mod)

    } else {

      modsum <- try(summary(mod), silent = TRUE)
      if ("try-error" %in% class(modsum)) {

        psm[["fit"]] <- FALSE
        modsum <- as.character(modsum)
        attributes(modsum) <- NULL
        psm[["mes"]] <- paste("On summary:", modsum)

      } else {

        if (is.nan(AIC(mod)) | is.infinite(AIC(mod)) | is.na(AIC(mod)) |
            any(is.nan(coef(mod))) | any(is.infinite(coef(mod))) |
            any(is.na(coef(mod)))) {

          psm[["fit"]] <- FALSE
          psm[["mes"]] <- paste("Strange values of AIC or coefficients.")

        } else {

          psm[["fit"]] <- TRUE

          psm[["aic"]] <- AIC(mod)

          t1 <- coef(mod) %>%
          {.[which(names(.) %in% dimnames(psm[["effects"]])[[1]])]} %>%
          {matrix(., ncol = 1)}
          t2 <- confint(mod, level = level) %>%
          {.[which(dimnames(.)[[1]] %in% dimnames(psm[["effects"]])[[1]]), ]} %>%
          {matrix(., ncol = 2)}
          psm[["effects"]] <- matrix(cbind(t1, t2), ncol = 3,
                                     dimnames = dimnames(psm[["effects"]]))

          if ("Log(scale)" %in% dimnames(summary(mod)$table)[[1]]) {

            psm[["logscale"]][1, 1:3] <- summary(mod)$table['Log(scale)', 1] +
              c(0, qnorm(lvn) * summary(mod)$table['Log(scale)', 2])

          }

          if (any(grepl("frailty(.+)", as.character(form)))) {

            psm[["frailtyvar"]] <- mod$history[[1]][[3]] %>%
              {.[[nrow(.), 1]]} %>%
              {. ^2}

          }

        }

      }

    }

    return(psm)

  }



  sfInit(parallel = parallel, cpus = cpus, type = type, sfInit_dots)
  if (parallel) {

    suppressMessages(sfLibrary(package = survival, verbose = FALSE))
    suppressMessages(sfLibrary(package = magrittr, verbose = FALSE))
    sfExport(list = list("Formula", "Data"))

  }
  ParResult <- sfLapply(Distributions, ParWrapper,...)
  sfStop()

  Result <- list(
    Formula = Formula,
    Distribution = sapply(ParResult, function(x) x[["dist"]]),
    Fit = sapply(ParResult, function(x) x[["fit"]]),
    Message = sapply(ParResult, function(x) x[["mes"]]),
    AIC = sapply(ParResult, function(x) x[["aic"]]),
    Effects = array(unlist(lapply(ParResult, function(x) x[["effects"]])),
                    dim = c(nrow(ParResult[[1]][["effects"]]),
                            3,
                            length(Distributions))),
    'logScale' = matrix(sapply(ParResult, function(x) x[["logscale"]]),
                          ncol = 3,
                          byrow = TRUE),
    FrailtyVar = sapply(ParResult, function(x) x[["frailtyvar"]])
  )
  Result <- within(Result, {
    names(Fit) <- Distribution
    names(Message) <- Distribution
    names(AIC) <- Distribution
    dimnames(Effects) <- list(dimnames(ParResult[[1]][["effects"]])[[1]],
                              dimnames(ParResult[[1]][["effects"]])[[2]],
                              Distribution)
    names(FrailtyVar) <- Distribution
  })
  dimnames(Result[["logScale"]]) <- list(Result[["Distribution"]],
                                   dimnames(ParResult[[1]][["logscale"]])[[2]])



  return(Result)

}


