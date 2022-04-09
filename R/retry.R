#' Retry function
#' @description Internal function for handling errors when accessing APIs
#' @param expr This is the function you want to catch and handle errors from
#' @param isError Function for evaluating if provided expression is throwing an error
#' @param maxErrors The maximum number of errrors it should handle from the function
#' @param sleep The amount of sleep between a caught error and the next attempt
#' @return The expression that has been either successfully ran or retried maximum number of times
#'
retry <- function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=5, sleep=0){
  ##https://stackoverflow.com/questions/20770497/how-to-retry-a-statement-on-error - code edited from asieira
  attempts <- 0
  retval <- try(eval(expr))
  while (isError(retval)) {
    attempts <- attempts + 1
    if (attempts >= maxErrors) {
      msg <- sprintf("retry: KEGGREST API might have recieved too many GET requests - try again with a lower core count or larger sleep time [[%s]]", utils::capture.output(utils::str(retval)))
      futile.logger::flog.fatal(msg)
      stop(msg)
    } else {
      msg <- sprintf("retry: error in attempt %i/%i [[%s]]", attempts, maxErrors,
                    utils::capture.output(utils::str(retval)))
      futile.logger::flog.error(msg)
      warning(msg)
    }
    if (sleep > 0) Sys.sleep(sleep)
    retval <- try(eval(expr))
  }
  return(retval)
}

