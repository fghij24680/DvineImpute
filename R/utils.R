#' @title Data transformation functions
#' @export
logit_100 <- function(data) return(log(data / 100 / (1 - data / 100)))

#' @export
inv_logit_100 <- function(data) return(100/(exp(-data)+1))

#' @title Error metrics
#' @export
rmse <- function(a,b) sqrt(mean((a-b)^2))

#' @export
rrse <- function(a,b) sqrt(sum((a-b)^2)/sum((mean(b)-b)^2))

#' @title Distribution helpers
u_fit <- function(data,fit)
  return(do.call(get(paste("p",fit$family[1],sep="")),c(list(data),fit_params(fit))))

q_fit <- function(data,fit)
  return(do.call(get(paste("q",fit$family[1],sep="")), c(list(data),fit_params(fit))))

