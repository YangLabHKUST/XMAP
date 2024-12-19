get_coef_raw <- function(ldscore_mat,sumstat1,sumstat2,max_int){
  #########
  # obtain raw coefficients for h2 / genetic covariance estimates
  #########
  p <- nrow(ldscore_mat)
  ncoef <- ncol(ldscore_mat)
  
  nprod <- sqrt(sumstat1$N)*sqrt(sumstat2$N)
  zprod <- sumstat1$Z*sumstat2$Z
  score_sum <- rowSums(matrix(ldscore_mat[,-1],nrow = p))
  
  # get averages
  mean_ldscore <- mean(score_sum)
  mean_zprod <- mean(zprod)
  mean_nprod <- mean(nprod)
  
  # estimate intercept
  zprod_sorted <- sort(zprod,decreasing = F)
  idx <- 0.95*p
  intercept <- mean(zprod_sorted[1:idx])
  intercept <- ifelse(intercept>max_int,max_int,intercept)
  
  # get raw coefs
  coef <- (mean_zprod-intercept) / (mean_nprod*mean_ldscore)
  
  return(list(coef=coef,intercept=intercept))
}

get_pred  <- function(coef,x,n1,n2,intercept){
  #########
  # obtain predicted chi-squared for constructing regression weight
  #########
  p <- nrow(x)
  
  nprod <- sqrt(n1)*sqrt(n2)
  score_sum <- rowSums(matrix(x[,-1],nrow = p))
  pred <- coef*score_sum*nprod + intercept
  
  return(pred)
}


update_weight <- function(w,pred){
  #########
  # update weight using predicted chi-squared for one trait
  #########
  var <- 2*pred^2
  
  w <- w/var
  
  return(w)
}

update_weightx <- function(w1,w2,wx,pred1,pred2,predx){
  #########
  # update weight using predicted chi-squared for two pops
  #########
  var1 <- 2*pred1^2
  var2 <- 2*pred2^2
  varx <- 2*pred1*pred2 + predx^2
  
  w1 <- w1/var1
  w2 <- w2/var2
  wx <- wx/varx
  
  return(list(w1=w1,w2=w2,wx=wx))
}

regression <- function(x,y,constrain_intercept,subtract,nblocks=200,jknife=T){
  #########
  # perform least squared regression to get coefs
  #########
  
  ncoef <- ncol(x)
  
  xtx <- t(x)%*%x
  xty <- t(x)%*%y
  
  # obtain coefs
  if(constrain_intercept){
    coefs <- solve(xtx[-1,-1],xty[-1])
    coefs <- c(subtract,coefs)
  } else {
    coefs <- solve(xtx,xty)
  }
  
  if(jknife){
    seperator <- floor(seq(from=1,to=length(y),length.out=nblocks+1))
    from <- seperator[-length(seperator)]
    to <- c(seperator[2:length(y)]-1,length(y))
    
    coefs_jk <- matrix(0,nblocks,ncoef)
    for(i in 1:nblocks){
      xtx_blk <- t(x[from[i]:to[i],])%*%x[from[i]:to[i],]
      xty_blk <- t(x[from[i]:to[i],])%*%y[from[i]:to[i]]
      
      xtx_loo <- xtx-xtx_blk
      xty_loo <- xty-xty_blk
      
      if(constrain_intercept){
        coefs_jk[i,-1] <- solve(xtx_loo[-1,-1],xty_loo[-1])
      } else{
        coefs_jk[i,] <- solve(xtx_loo,xty_loo)
      }
    }
    jk_cov <- cov(coefs_jk)*(nblocks-1)#cov(t(nblocks*c(coefs)-t((nblocks-1)*coefs_jk))) / nblocks
    jk_se <- sqrt(diag(jk_cov))
  } else {
    jk_se <- NA
  }
  
  return(list(coefs=coefs,
              coefs_se=jk_se))
}

get_coef <- function(score,sumstat1,sumstat2,w,constrain_intercept,subtract,nblocks=200,jknife=T){
  #########
  # wrapper function to get coefs
  #########
  p <- nrow(score)
  nprod <- sqrt(sumstat1$N)*sqrt(sumstat2$N)
  zprod <- sumstat1$Z*sumstat2$Z
  if(constrain_intercept){
    zprod <- zprod-subtract
  }
  score[,-1] <- score[,-1]*nprod
  
  # scale the matrix to improve matrix condition
  nbar <- mean(nprod)
  score[,-1] <- score[,-1]/nbar
  
  # apply weight to data
  score_w <- score*sqrt(w)
  zprod_w <- zprod*sqrt(w)
  
  # perform regression
  reg <- regression(score_w,zprod_w,constrain_intercept,subtract,nblocks = nblocks,jknife = jknife)
  
  # re-scale coefs
  reg$coefs [-1] <- reg$coefs[-1]/nbar
  
  if(jknife){
    reg$coefs_se[-1] <- reg$coefs_se[-1]/nbar
  }
  
  return(reg)
}

estimate_intercept <- function(chiSQ,score,ld_w){
  pred <- (mean(chiSQ)-1)/mean(score)*score+1
  w <- ld_w/(2*pred^2)
  x <- cbind(1,score)*sqrt(w)
  y <- chiSQ*sqrt(w)
  
  # don't need n here, slope has no interpretation
  intercept <- solve(t(x)%*%x,t(x)%*%y)[1]
  return(intercept)
}

estimate_gc <- function(sumstat1,sumstat2,ldscore1,ldscore2,ldscorex,
                        reg_w1=1,reg_w2=1,reg_wx=1,constrain_intercept=F,int1=1,int2=1,intx=0,
                        nblocks=200,jknife=T){
  
  # add an additional column for intercept
  ldscore1 <- cbind(1,ldscore1)
  ldscore2 <- cbind(1,ldscore2)
  ldscorex <- cbind(1,ldscorex)
  
  # get initial estimate of coefs for constructing weights
  raw_coef1 <- get_coef_raw(ldscore1,sumstat1,sumstat1,1)
  raw_coef2 <- get_coef_raw(ldscore2,sumstat2,sumstat2,1)
  raw_coefx <- get_coef_raw(ldscorex,sumstat1,sumstat2,1)
  
  # get predicted chi-squared for constructing weights
  n1 <- sumstat1$N
  n2 <- sumstat2$N
  pred1 <- get_pred(raw_coef1$coef,ldscore1,n1,n1,raw_coef1$intercept)
  pred2 <- get_pred(raw_coef2$coef,ldscore2,n2,n2,raw_coef2$intercept)
  predx <- get_pred(raw_coefx$coef,ldscorex,n1,n2,raw_coefx$intercept)
  
  # update weight
  reg_w <- update_weightx(reg_w1,reg_w2,reg_wx,pred1,pred2,predx)
  
  # conduct step 1 regression to obtain coefs
  coef1 <- get_coef(ldscore1,sumstat1,sumstat1,reg_w$w1,constrain_intercept,int1,nblocks = nblocks,jknife = jknife)
  coef2 <- get_coef(ldscore2,sumstat2,sumstat2,reg_w$w2,constrain_intercept,int2,nblocks = nblocks,jknife = jknife)
  coefx <- get_coef(ldscorex,sumstat1,sumstat2,reg_w$wx,constrain_intercept,intx,nblocks = nblocks,jknife = jknife)
  
  return(list(tau1=coef1,tau2=coef2,theta=coefx))
}



estimate_h2 <- function(sumstat,ldscore,reg_w=1,constrain_intercept=F,int=1){
  
  # add an additional column for intercept
  ldscore <- cbind(1,ldscore)
  
  # get initial estimate of coefs for constructing weights
  raw_coef <- get_coef_raw(ldscore,sumstat,sumstat,1)
  
  # get predicted chi-squared for constructing weights
  n <- sumstat$N
  pred <- get_pred(raw_coef$coef,ldscore,n,n,raw_coef$intercept)
  
  # update weight
  reg_w <- update_weight(reg_w,pred)
  
  # conduct step 1 regression to obtain coefs
  coef <- get_coef(ldscore,sumstat,sumstat,reg_w,constrain_intercept,int)
  
  return(coef)
}
