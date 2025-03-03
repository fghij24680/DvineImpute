#' Multivariate imputation
#' Check monotonic missing pattern
#' @export
is.mono <- function(x,plot=FALSE){
  
  pt = mice::md.pattern(x,plot = plot)
  
  pt = pt[1:nrow(pt)-1,1:ncol(pt)-1]
  
  mono_flag = TRUE
  
  for(i in seq(1,nrow(pt)-1,1)){
    if(!(
      dplyr::count(pt[i,],1)==dplyr::count(pt[i+1,],1) | 
      dplyr::count(pt[i,],1)==dplyr::count(pt[i+1,],1)+1))
    {mono_flag = FALSE}
  }
  
  for(i in seq(1,ncol(pt)-1,1)){
    if(!(
      dplyr::count(pt[,i],1)==dplyr::count(pt[,i+1],1) | 
      dplyr::count(pt[,i],1)==dplyr::count(pt[,i+1],1)+1))
    {mono_flag = FALSE}
  }
  
  return(mono_flag)
}

#' Multi-variable imputation
#' @param order Imputation order (seq/hire)
#' @param method Imputation method (dvinemed/dvinemean/lm/rq/rlm/lqr)
#' @export
multi.impute <- function(x,order='seq',method='dvinemed'){
  if(!is.mono(x)){
    print("Data is not monotone missing in the correct order!")
    return(FALSE)
  }else if(method=='dvinemed' | method=='dvinemean'){
    return(multi.impute.dvine(x,order=order,method=method))
  }else if(order=='seq'){
    # Sequential imputation
    x_impute = x
    for(i in seq(ncol(x)-1,1,-1)){
      if(sum(is.na(x_impute[,i]))!=0)
        x_impute[,i:ncol(x)] = 
          uni.impute(x_impute[,i:ncol(x)],method=method)
    }
    return(x_impute)
  }else if(order=='hierarchy'){
    # Hierarchical imputation
    x_impute = x
    for(j in seq(1,ncol(x)-1,1)){
      for(i in seq(1,ncol(x)-1,1)){
        if(sum(is.na(x_impute[,i]))!=0){
          miss_all = which(is.na(x_impute[,i]))
          to_impute = 
            miss_all[miss_all %in% which(!is.na(x_impute[,i+1]))]
          model_impute = c(which(!is.na(x_impute[,i])),to_impute)
          x_impute[model_impute,i:ncol(x)] = 
            uni.impute(x_impute[model_impute,i:ncol(x)],method=method)
        }
      }
    }
    return(x_impute)
  }else {
    print("Please provide correct parameters!")
    return(FALSE)
  }
}

#' Multi-variable D-vine imputation
#' @export
multi.impute.dvine <- function(x,order='seq',method='dvinemed'){
  fit_list = list()
  for(i in seq(1,ncol(x),1)){
    fit_dist = gamlss::fitDist(unlist(x[,i]), type="realline")
    fit_list = append(fit_list,list(fit_dist))
  }
  
  x_impute = x
  u_impute = x
  
  for(i in seq(1,ncol(x),1)){
    u_impute[which(!is.na(x[,i])),i] = u_fit(unlist(x[which(!is.na(x[,i])),i]),fit_list[[i]])
  }
  
  if(order=='seq'){
    # Sequential imputation
    for(i in seq(ncol(x)-1,1,-1)){
      if(sum(is.na(u_impute[,i]))!=0)
        u_impute[,i:ncol(x)] = 
          uni.impute(u_impute[,i:ncol(x)],method=method,uscale=TRUE)
    }
  }else if(order=='hierarchy'){
    # Hierarchical imputation
    for(j in seq(1,ncol(x)-1,1)){
      for(i in seq(1,ncol(x)-1,1)){
        if(sum(is.na(u_impute[,i]))!=0){
          miss_all = which(is.na(u_impute[,i]))
          to_impute = 
            miss_all[miss_all %in% which(!is.na(u_impute[,i+1]))]
          model_impute = c(which(!is.na(u_impute[,i])),to_impute)
          u_impute[model_impute,i:ncol(x)] = 
            uni.impute(u_impute[model_impute,i:ncol(x)],method=method,uscale=TRUE)
        }
      }
    }
  }else {
    print("Please provide correct parameters!")
    return(FALSE)
  }
  
  for(i in seq(1,ncol(x),1)){
    x_impute[,i] = q_fit(unlist(u_impute[,i]),fit_list[[i]])
  }
  
  return(x_impute)
}

#' Multi-variable imputation interval
#' @export
multi.impute.q <- function(x,order='seq',method='dvinemed',q=0.5){
  if(!is.mono(x)){
    print("Data is not monotone missing in the correct order!")
    return(FALSE)
  }else if(method=='dvinemed' | method=='dvinemean'){
    return(multi.impute.q.dvine(x,order=order,method=method,q=q))
  }else if(order=='seq'){
    # Sequential imputation
    x_impute = x
    for(i in seq(ncol(x)-1,1,-1)){
      if(sum(is.na(x_impute[,i]))!=0)
        x_impute[,i:ncol(x)] = 
          uni.impute.q(x_impute[,i:ncol(x)],method=method)
    }
    if(q==0.5){
      return(x_impute)
    } else {
      x_quantile = x
      for(i in seq(ncol(x)-1,1,-1)){
        if(sum(is.na(x_quantile[,i]))!=0){
          missvec = seq(length(x[,i]))[is.na(x[,i])]
          x_impute_tmp = x_impute
          x_impute_tmp[missvec,i] = NA
          x_quantile[,i] = 
            uni.impute.q(x_impute_tmp[,i:ncol(x)],method=method,q=q)[,1]
        }
      }
      return(x_quantile)
    }
  }else if(order=='hierarchy'){
    # Hierarchical imputation
    x_impute = x
    for(j in seq(1,ncol(x)-1,1)){
      for(i in seq(1,ncol(x)-1,1)){
        if(sum(is.na(x_impute[,i]))!=0){
          miss_all = which(is.na(x_impute[,i]))
          to_impute = 
            miss_all[miss_all %in% which(!is.na(x_impute[,i+1]))]
          model_impute = c(which(!is.na(x_impute[,i])),to_impute)
          x_impute[model_impute,i:ncol(x)] =
            uni.impute.q(x_impute[model_impute,i:ncol(x)],method=method)
        }
      }
    }
    if(q==0.5){
      return(x_impute)
    } else {
      x_quantile = x
      for(j in seq(1,ncol(x)-1,1)){
        for(i in seq(1,ncol(x)-1,1)){
          if(sum(is.na(x_quantile[,i]))!=0){
            miss_all = which(is.na(x_quantile[,i]))
            to_impute = 
              miss_all[miss_all %in% which(!is.na(x_quantile[,i+1]))]
            x_impute_tmp = x_impute
            x_impute_tmp[miss_all,i] = NA
            # model_impute = c(which(!is.na(x_quantile[,i])),to_impute)
            # x_quantile[to_impute,i] = 
            #   uni.impute.q(x_impute_tmp[model_impute,i:ncol(x)],method=method,q=q)[to_impute,1]
            x_quantile[to_impute,i] = 
              uni.impute.q(x_impute_tmp[,i:ncol(x)],method=method,q=q)[to_impute,1]
          }
        }
      }
      return(x_quantile)
    }
  }else {
    print("Please provide correct parameters!")
    return(FALSE)
  }
}

#' Multi-variable D-vine imputation interval
#' @export
multi.impute.q.dvine <- function(x,order='seq',method='dvinemed',q=0.5){
  fit_list = list()
  for(i in seq(1,ncol(x),1)){
    fit_dist = gamlss::fitDist(unlist(x[,i]), type="realline")
    fit_list = append(fit_list,list(fit_dist))
  }
  
  x_impute = x
  u_impute = x
  
  for(i in seq(1,ncol(x),1)){
    u_impute[which(!is.na(x[,i])),i] = u_fit(unlist(x[which(!is.na(x[,i])),i]),fit_list[[i]])
  }
  
  u=u_impute
  
  if(order=='seq'){
    # Sequential imputation
    for(i in seq(ncol(x)-1,1,-1)){
      if(sum(is.na(u_impute[,i]))!=0)
        u_impute[,i:ncol(x)] = 
          uni.impute.q(u_impute[,i:ncol(x)],method=method,uscale=TRUE)
    }
    if(q==0.5){
      u_quantile = u_impute
    } else {
      u_quantile = u
      for(i in seq(ncol(x)-1,1,-1)){
        if(sum(is.na(u_quantile[,i]))!=0){
          missvec = seq(length(x[,i]))[is.na(x[,i])]
          u_impute_tmp = u_impute
          u_impute_tmp[missvec,i] = NA
          u_quantile[,i] = 
            uni.impute.q(u_impute_tmp[,i:ncol(x)],method=method,uscale=TRUE,q=q)[,1]
        }
      }
    }
  }else if(order=='hierarchy'){
    # Hierarchical imputation
    for(j in seq(1,ncol(x)-1,1)){
      for(i in seq(1,ncol(x)-1,1)){
        if(sum(is.na(u_impute[,i]))!=0){
          miss_all = which(is.na(u_impute[,i]))
          to_impute = 
            miss_all[miss_all %in% which(!is.na(u_impute[,i+1]))]
          model_impute = c(which(!is.na(u_impute[,i])),to_impute)
          u_impute[model_impute,i:ncol(x)] =
            uni.impute.q(u_impute[model_impute,i:ncol(x)],method=method,uscale=TRUE)
        }
      }
    }
    if(q==0.5){
      u_quantile = u_impute
    } else {
      u_quantile = u
      for(j in seq(1,ncol(x)-1,1)){
        for(i in seq(1,ncol(x)-1,1)){
          if(sum(is.na(u_quantile[,i]))!=0){
            miss_all = which(is.na(u_quantile[,i]))
            to_impute = 
              miss_all[miss_all %in% which(!is.na(u_quantile[,i+1]))]
            u_impute_tmp = u_impute
            u_impute_tmp[miss_all,i] = NA
            # model_impute = c(which(!is.na(x_quantile[,i])),to_impute)
            # x_quantile[to_impute,i] = 
            #   uni.impute.q(x_impute_tmp[model_impute,i:ncol(x)],method=method,q=q)[to_impute,1]
            u_quantile[to_impute,i] = 
              uni.impute.q(u_impute_tmp[,i:ncol(x)],method=method,uscale=TRUE,q=q)[to_impute,1]
          }
        }
      }
    }
  }else {
    print("Please provide correct parameters!")
    return(FALSE)
  }
  
  x_quantile = x
  
  for(i in seq(1,ncol(x),1)){
    x_quantile[,i] = q_fit(unlist(u_quantile[,i]),fit_list[[i]])
  }
  
  return(x_quantile)
}
