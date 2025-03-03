#' Univariate imputation
#' @param x Dataframe with missing values
#' @param uscale Dataframe is u-scale or x-scale
#' @param print_pair If need to print pairs plot

#' Internal function of getting fitted parameter
#' @param fit A gamlss object
fit_params <- function(fit){
  fit_param=vector()
  for(i in seq(1:length(fit$parameters))) fit_param = append(fit_param,fit[fit$parameters[i]])
  return(fit_param)
}

#' D-vine median imputation
#' @export
uni.impute.dvinemed <- function(x,uscale=FALSE,print_pair=FALSE,alpha=0.5){
  missing_vec = seq(length(x[,1]))[is.na(x[,1])]
  x_train = x[-missing_vec,]
  
  u_train = x_train
  
  x_pred = x[missing_vec,-1]
  u_pred = x_pred
  
  if(!uscale){
    # Fit marginal distribution
    fit_list = list()
    for(i in seq(1,ncol(x_train),1)){
      fit_dist = gamlss::fitDist(unlist(x_train[,i]), type="realline")
      fit_list = append(fit_list,list(fit_dist))
    }
    
    # Transform to u-scale
    for(i in seq(1,ncol(x_train),1)){
      u_train[,i] = u_fit(x_train[,i],fit_list[[i]])
    }
    
    if(print_pair) rvinecoplib::pairs_copula_data(x_train,col='red')
    if(is.null(ncol(x_pred))){u_pred = u_fit(x_pred,fit_list[[1+1]])}
    else{
      for(i in seq(1,ncol(x_pred),1)){
        u_pred[,i] = u_fit(x_pred[,i],fit_list[[i+1]])
      }
    }
  }
  
  # Fit model and impute
  fitstr = paste(colnames(x)[1], "~")
  for(i in seq(2,ncol(x),1)){
    fitstr = paste(fitstr,"+",colnames(x)[i])
  }
  
  fit = vinereg::vinereg(stats::as.formula(fitstr),data=u_train, uscale = TRUE)
  if(is.null(ncol(x_pred))){
    u_pred = data.frame(u_pred)
    names(u_pred) = colnames(x)[2]
  }
  
  u_impute = stats::predict(fit, newdata = u_pred, alpha = alpha)
  if(!uscale)
    x_impute = q_fit(unlist(u_impute),fit_list[[1]])
  else 
    x_impute = unlist(u_impute)
  
  x[missing_vec,1] = x_impute
  return(x)
}

#' D-vine mean imputation
#' @export
uni.impute.dvinemean <- function(x,uscale=FALSE,print_pair=FALSE){
  missing_vec = seq(length(x[,1]))[is.na(x[,1])]
  x_train = x[-missing_vec,]
  
  u_train=x_train
  
  x_pred = x[missing_vec,-1]
  u_pred = x_pred
  
  if(!uscale){
    # Fit marginal distribution
    fit_list = list()
    for(i in seq(1,ncol(x_train),1)){
      fit_dist = gamlss::fitDist(unlist(x_train[,i]), type="realline")
      fit_list = append(fit_list,list(fit_dist))
    }
    
    # Transform to u-scale
    
    for(i in seq(1,ncol(x_train),1)){
      u_train[,i] = u_fit(x_train[,i],fit_list[[i]])
    }
    
    if(print_pair) rvinecoplib::pairs_copula_data(x_train,col='red')
    if(is.null(ncol(x_pred))){u_pred = u_fit(x_pred,fit_list[[1+1]])}
    else{
      for(i in seq(1,ncol(x_pred),1)){
        u_pred[,i] = u_fit(x_pred[,i],fit_list[[i+1]])
      }
    }
  }
  
  
  # Fit model and impute
  fitstr = paste(colnames(x)[1], "~")
  for(i in seq(2,ncol(x),1)){
    fitstr = paste(fitstr,"+",colnames(x)[i])
  }
  fit = vinereg::vinereg(stats::as.formula(fitstr),data=u_train, uscale = TRUE)
  if(is.null(ncol(x_pred))){
    u_pred = data.frame(u_pred)
    names(u_pred) = colnames(x)[2]
  }
  u_impute = stats::predict(fit, newdata = u_pred, alpha = NA)
  if(!uscale)
    x_impute = q_fit(unlist(u_impute),fit_list[[1]])
  else 
    x_impute = unlist(u_impute)
  x[missing_vec,1] = x_impute
  return(x)
}

#' Linear regression imputation
#' @export
uni.impute.lm <- function(x,print_pair=FALSE){
  missing_vec = seq(length(x[,1]))[is.na(x[,1])]
  if(print_pair) GGally::ggpairs(x[-missing_vec])
  fit = stats::lm(stats::as.formula(paste(colnames(x)[1], "~.")),data=x[-missing_vec,])
  x_impute = stats::predict(fit,newdata = x[missing_vec,])
  x[missing_vec,1] = x_impute
  return(x)
}

#' KNN imputation
#' @export
uni.impute.knn <- function(x,print_pair=FALSE){
  missing_vec = seq(length(x[,1]))[is.na(x[,1])]
  if(print_pair) GGally::ggpairs(x[-missing_vec])
  return(data.frame(impute::impute.knn(as.matrix(x))$data))
}

#' Linear quantile regression imputation
#' @export
uni.impute.rq <- function(x,print_pair=FALSE,tau=0.5){
  missing_vec = seq(length(x[,1]))[is.na(x[,1])]
  if(print_pair) GGally::ggpairs(x[-missing_vec])
  fit = quantreg::rq(stats::as.formula(paste(colnames(x)[1], "~.")),data=x[-missing_vec,],tau=tau)
  x_impute = stats::predict(fit,newdata = x[missing_vec,])
  x[missing_vec,1] = x_impute
  return(x)
}

#' Robust regression imputation
#' @export
uni.impute.rlm <- function(x,print_pair=FALSE){
  missing_vec = seq(length(x[,1]))[is.na(x[,1])]
  if(print_pair) GGally::ggpairs(x[-missing_vec])
  fit = MASS::rlm(stats::as.formula(paste(colnames(x)[1], "~.")),data=x[-missing_vec,])
  x_impute = stats::predict(fit,newdata = x[missing_vec,])
  x[missing_vec,1] = x_impute
  return(x)
}

#' Robustr quantile regression imputation
#' @export
uni.impute.lqr <- function(x,print_pair=FALSE,p=0.5){
  missing_vec = seq(length(x[,1]))[is.na(x[,1])]
  if(print_pair) GGally::ggpairs(x[-missing_vec])
  fit_lqr = lqr::lqr(stats::as.formula(paste(colnames(x)[1], "~.")),data=x[-missing_vec,],p=p)
  x_impute = fit_lqr$beta[1]+rowSums(x[missing_vec,-1, drop = FALSE]*
                                       fit_lqr$beta[-1][col(x[missing_vec,-1, drop = FALSE])])
  x[missing_vec,1] = x_impute
  return(x)
}

#' Main imputation interface
#' @param method Imputation method (dvinemed/dvinemean/lm/rq/rlm/lqr)
#' @import gamlss
#' @export
uni.impute <- function(x,method='dvinemed',uscale=FALSE,print_pair=FALSE){
  if(method == 'dvinemed'){
    return(uni.impute.dvinemed(x,uscale,print_pair))
  } else if(method == 'dvinemean'){
    return(uni.impute.dvinemean(x,uscale,print_pair))
  } else if(method == 'lm'){
    return(uni.impute.lm(x,print_pair))
  } else if(method == 'knn'){
    return(uni.impute.knn(x,print_pair))
  } else if(method == 'rq'){
    return(uni.impute.rq(x,print_pair))
  } else if(method == 'rlm'){
    return(uni.impute.rlm(x,print_pair))
  } else {
    print("Please use the correct imputation method!")
    return(FALSE)
  }
}

#' Main imputation interface for interval
#' @export
uni.impute.q <- function(x,method='dvinemed',uscale=FALSE,print_pair=FALSE,q=0.5){
  if(method == 'dvinemed'){
    return(uni.impute.dvinemed(x,uscale,print_pair,alpha=q))
  } else if(method == 'rq'){
    return(uni.impute.rq(x,print_pair,tau=q))
  } else if(method == 'lqr'){
    return(uni.impute.lqr(x,print_pair,p=q))
  } else {
    print("Please use the correct imputation method!")
    return(FALSE)
  }
}