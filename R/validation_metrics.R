#' Validation metrics
#' Calculate RMSE
#' @param y Dataframe without missing values
#' @param y_005 5% percentile imputation
#' @param y_050 Median imputation
#' @param y_095 95% percentile imputation
#' @export
rmse <- function(a,b) return(sqrt(mean((a-b)^2)))

#' Calculate RRSE
#' @export
rrse <- function(a,b){
  return(sqrt(sum((a-b)^2)/sum((mean(b)-b)^2)))
}

#' Calculate relative coverage probability
#' @export
rcp <- function(y,y_005,y_095,alpha=0.9) {
  covered <- (y >= y_005) & (y <= y_095)
  return(abs(alpha-mean(covered)))
}

#' Calculate coverage probability
#' @export
coverage_prob <- function(y,y_005,y_095) {
  covered <- (y >= y_005) & (y <= y_095)
  return(mean(covered))
}

#' Calculate sharpness
#' @export
sharpness <- function(y_005, y_095) {
  return(mean(y_095 - y_005))
}

#' Calculate interval score
#' @export
interval_score <- function(y, y_005, y_095, alpha = 0.05) {
  interval_width <- y_095 - y_005
  penalty_below <- (y_005 - y) * (y < y_005)
  penalty_above <- (y - y_095) * (y > y_095)
  return(mean(interval_width + (2 / alpha) * (penalty_below + penalty_above)))
}

#' Calculate CRPS
#' @export
crps <- function(y,y_005,y_050,y_095){
  crps_values <- sapply(1:length(y), function(i) {
    scoringRules::crps_sample(
      y[i],c(y_005[i], y_050[i], y_095[i]) 
    )
  })
  mean_crps <- mean(crps_values)
  return(mean_crps)
}

#' Main scoring metrics interface
#' @export
scoring_metrics <- function(y,y_005,y_050,y_095,alpha = 0.05){
  a=coverage_prob(y,y_005,y_095)
  b=sharpness(y_005, y_095)
  c=interval_score(y, y_005, y_095, alpha)
  d=crps(y,y_005,y_050,y_095)
  return(c(a,b,c,d))
}
