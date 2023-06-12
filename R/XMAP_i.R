XMAP_i <- function(X1, X2, y1, y2, K = 10, K1 = 0, K2 = 0, Sigma1, Sigma2, rho, Eta1 = NULL, Eta2 = NULL, Omega, Sig_E1, Sig_E2, prior_weights, prior_weights1 = NULL, prior_weights2 = NULL, maxIter = 1000, tol = 5e-6, estimate_prior_variance = T, estimate_residual_variance = T, estimate_background_variance = F, initialize = F) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  n <- n1 + n2
  
  p <- ncol(X1)
  if (p != ncol(X2)) {
    stop("X1 and X2 musth have the same number of predictors!")
  }
  
  X1y1 <- t(X1) %*% y1
  X2y2 <- t(X2) %*% y2
  
  X1X1 <- t(X1) %*% X1
  X2X2 <- t(X2) %*% X2
  
  
  if (length(Sigma1) == 1) {
    Sigma1 <- rep(Sigma1, K)
    Sigma1 <- Sigma1 / (1:K)
  }
  if (length(Sigma2) == 1) {
    Sigma2 <- rep(Sigma2, K)
    Sigma2 <- Sigma2 / (1:K)
  }
  if (length(rho) == 1) rho <- rep(rho, K)
  
  
  if (length(Eta1) == 1) {
    Eta1 <- rep(Eta1, K1)
    Eta1 <- Eta1 / (1:K1)
  }
  if (length(Eta2) == 1) {
    Eta2 <- rep(Eta2, K2)
    Eta2 <- Eta2 / (1:K2)
  }
  
  
  Sigma <- array(rbind(Sigma1, rho * sqrt(Sigma1 * Sigma2), rho * sqrt(Sigma1 * Sigma2), Sigma2), dim = c(2, 2, K))
  invSigma <- array(apply(Sigma, 3, solve), c(2, 2, K))
  
  
  gamma <- matrix(1 / p, K, p)              # posterior probability of shared effects K by p
  mu1 <- mu2 <- matrix(0, K, p)     # posterior mean of shared effects K by p
  
  post_Sig <- aperm(replicate(p, Sigma), c(3, 4, 2, 1))    # posterior covariance of shared effects K by p by 2 by 2
  
  sum_gammamu1 <- colSums(mu1 * gamma)
  sum_gammamu2 <- colSums(mu2 * gamma)
  
  post_mumu <- list()                 # posterior second monent of beta each component will be 2*2*p matrix
  
  Lambdaj <- matrix(0, 2, 2)
  
  if (K1 > 0) {
    tgamma1 <- matrix(1 / p, K1, p)           # posterior probability of pop1-specific effects K1 by p
    tmu1 <- matrix(0, K1, p)          # posterior mean of pop1-specific effects K1 by p
    post_Eta1 <- replicate(p, Eta1)      # posterior variance of pop1-specific effects K1 by p
    sum_gammamu1 <- sum_gammamu1 + colSums(tmu1 * tgamma1)
    post_mumu1 <- matrix(0, K1, p)  # posterior second monent of alpha1 each component will be p vector
  } else {
    tgamma1 <- tmu1 <- NULL
  }
  
  if (K2 > 0) {
    tgamma2 <- matrix(1 / p, K2, p)           # posterior probability of pop2-specific effects K2 by p
    tmu2 <- matrix(0, K2, p)          # posterior mean of pop2-specific effects K2 by p
    post_Eta2 <- replicate(p, Eta2)      # posterior variance of pop2-specific effects K2 by p
    sum_gammamu2 <- sum_gammamu2 + colSums(tmu2 * tgamma2)
    post_mumu2 <- matrix(0, K2, p)  # posterior second monent of alpha2 each component will be p vector
  } else {
    tgamma2 <- tmu2 <- NULL
  }
  
  nuj <- matrix(0, p, 2)
  
  Xnu1 <- X1 %*% nuj[,1]
  Xnu2 <- X2 %*% nuj[,2]
  
  ELBO <- rep(-Inf, maxIter)
  # ELBO[1] <- get_ELBO(n1,n2,X1,X2,y1,y2,K,K1,K2,invSigma,Eta1,Eta2,Sig_E1,Sig_E2,prior_weights,prior_weights1,prior_weights2,X1X1,X2X2,gamma,tgamma1,tgamma2,mu1,mu2,tmu1,tmu2,post_Sig,post_Eta1,post_Eta2,sum_gammamu1,sum_gammamu2)
  # cat("Initial ELBO value:",ELBO[1],".\n")
  
  
  ########################## Initialization ############################
  if (initialize) {
    update_times <- rep(0, 3)
    # E-step for shared effects
    for (k in 1:(K * 3)) {
      if (update_times[1] < 5) {
        k12 <- update_times[1] + 1
        sum_gammamu1k <- sum_gammamu1 - (mu1[k12,] * gamma[k12,])
        sum_gammamu2k <- sum_gammamu2 - (mu2[k12,] * gamma[k12,])
        
        gamma_k <- rep(0, p)
        post_mumu[[k12]] <- array(0, c(2, 2, p)) # posterior second monent of beta each component is 2*2*p matrix
        for (j in 1:p) {
          tmp1 <- X1y1[j] - sum(X1X1[, j] * sum_gammamu1k)
          tmp2 <- X2y2[j] - sum(X2X2[, j] * sum_gammamu2k)
          invSig_kj <- invSigma[, , k12] + diag(c(X1X1[j, j] / Sig_E1, X2X2[j, j] / Sig_E2))
          post_Sig[k12, j, ,] <- solve(invSig_kj)
          mu_kj <- post_Sig[k12, j, ,] %*% c(tmp1 / Sig_E1, tmp2 / Sig_E2)
          mu1[k12, j] <- mu_kj[1]
          mu2[k12, j] <- mu_kj[2]
          
          gamma_k[j] <- log(prior_weights[j]) - log_det(invSig_kj) / 2 + t(mu_kj) %*% invSig_kj %*% mu_kj / 2
          
          post_mumu[[k12]][, , j] <- mu_kj %*% t(mu_kj) + post_Sig[k12, j, ,]
          
        }
        gamma_k <- gamma_k - max(gamma_k)
        gamma[k12,] <- exp(gamma_k) / sum(exp(gamma_k))
        # gamma[gamma > 0.999999999] <- 0.999999999
        # gamma[gamma < 0.000000001] <- 0.000000001
        
        # sum_gammamu1 <- sum_gammamu1k+(mu1[k,]*gamma[k,])
        # sum_gammamu2 <- sum_gammamu2k+(mu2[k,]*gamma[k,])
        ELBO12 <- get_ELBO(n1, n2, X1, X2, y1, y2, K, K1, K2, invSigma, Eta1, Eta2, Sig_E1, Sig_E2, prior_weights, prior_weights1, prior_weights2, X1X1, X2X2, gamma, tgamma1, tgamma2, mu1, mu2, tmu1, tmu2, post_Sig, post_Eta1, post_Eta2, sum_gammamu1k + (mu1[k12,] * gamma[k12,]), sum_gammamu2k + (mu2[k12,] * gamma[k12,]))
        
      }else {
        ELBO12 <- -Inf
      }
      # cat("Initial iteration ELBO=", ELBO12, "\n")
      
      
      if (update_times[2] < 5) {
        k1 <- update_times[2] + 1
        sum_gammamu1k <- sum_gammamu1 - (tmu1[k1,] * tgamma1[k1,])
        
        tgamma1_k <- rep(0, p)
        for (j in 1:p) {
          tmp1 <- X1y1[j] - sum(X1X1[, j] * sum_gammamu1k)
          post_Eta1[k1, j] <- 1 / (1 / Eta1[k1] + X1X1[j, j] / Sig_E1)
          tmu1[k1, j] <- post_Eta1[k1, j] * tmp1 / Sig_E1
          
          tgamma1_k[j] <- log(prior_weights1[j]) +
            log(post_Eta1[k1, j]) / 2 +
            tmu1[k1, j]^2 / post_Eta1[k1, j] / 2
          
          post_mumu1[k1, j] <- tmu1[k1, j]^2 + post_Eta1[k1, j]
          
        }
        tgamma1_k <- tgamma1_k - max(tgamma1_k)
        tgamma1[k1,] <- exp(tgamma1_k) / sum(exp(tgamma1_k))
        # gamma[gamma > 0.999999999] <- 0.999999999
        # gamma[gamma < 0.000000001] <- 0.000000001
        
        # sum_gammamu1 <- sum_gammamu1k+(tmu1[k,]*tgamma1[k,])
        ELBO1 <- get_ELBO(n1, n2, X1, X2, y1, y2, K, K1, K2, invSigma, Eta1, Eta2, Sig_E1, Sig_E2, prior_weights, prior_weights1, prior_weights2, X1X1, X2X2, gamma, tgamma1, tgamma2, mu1, mu2, tmu1, tmu2, post_Sig, post_Eta1, post_Eta2, sum_gammamu1k + (tmu1[k1,] * tgamma1[k1,]), sum_gammamu2)
      }else {
        ELBO1 <- -Inf
      }
      # cat("Initial iteration ELBO=", ELBO1, "\n")
      
      if (update_times[3] < 5) {
        k2 <- update_times[3] + 1
        sum_gammamu2k <- sum_gammamu2 - (tmu2[k2,] * tgamma2[k2,])
        
        tgamma2_k <- rep(0, p)
        for (j in 1:p) {
          tmp2 <- X2y2[j] - sum(X2X2[, j] * sum_gammamu2k)
          post_Eta2[k2, j] <- 1 / (1 / Eta2[k2] + X2X2[j, j] / Sig_E2)
          tmu2[k2, j] <- post_Eta2[k2, j] * tmp2 / Sig_E2
          
          tgamma2_k[j] <- log(prior_weights2[j]) +
            log(post_Eta2[k2, j]) / 2 +
            tmu2[k2, j]^2 / post_Eta2[k2, j] / 2
          
          post_mumu2[k2, j] <- tmu2[k2, j]^2 + post_Eta2[k2, j]
          
        }
        tgamma2_k <- tgamma2_k - max(tgamma2_k)
        tgamma2[k2,] <- exp(tgamma2_k) / sum(exp(tgamma2_k))
        # gamma[gamma > 0.999999999] <- 0.999999999
        # gamma[gamma < 0.000000001] <- 0.000000001
        
        # sum_gammamu2 <- sum_gammamu2k+(tmu2[k,]*tgamma2[k,])
        ELBO2 <- get_ELBO(n1, n2, X1, X2, y1, y2, K, K1, K2, invSigma, Eta1, Eta2, Sig_E1, Sig_E2, prior_weights, prior_weights1, prior_weights2, X1X1, X2X2, gamma, tgamma1, tgamma2, mu1, mu2, tmu1, tmu2, post_Sig, post_Eta1, post_Eta2, sum_gammamu1, sum_gammamu2k + (tmu2[k2,] * tgamma2[k2,]))
      }else {
        ELBO2 <- -Inf
      }
      # cat("Initial iteration ELBO=", ELBO2, "\n")
      
      ELBO_init <- c(ELBO12, ELBO1, ELBO2)
      cat("Initialization step ELBO: ", ELBO12, ",\t", ELBO1, ",\t", ELBO2, ". Update ", c("shared", "pop-1", "pop-2")[which.max(ELBO_init)], " component \n")
      l <- c(3, 1, 2)[which.max(ELBO_init)]
      if (l == 3) {
        update_times[1] <- update_times[1] + 1
        k12 <- update_times[1]
        sum_gammamu1k <- sum_gammamu1 - (mu1[k12,] * gamma[k12,])
        sum_gammamu2k <- sum_gammamu2 - (mu2[k12,] * gamma[k12,])
        
        gamma_k <- rep(0, p)
        post_mumu[[k12]] <- array(0, c(2, 2, p)) # posterior second monent of beta each component is 2*2*p matrix
        for (j in 1:p) {
          tmp1 <- X1y1[j] - sum(X1X1[, j] * sum_gammamu1k)
          tmp2 <- X2y2[j] - sum(X2X2[, j] * sum_gammamu2k)
          # invSig_kj <- invSigma[,,k]+diag(c(X1X1[j,j]/Sig_E1,X2X2[j,j]/Sig_E2))
          # mu_kj <- solve(invSig_kj,c(tmp1/Sig_E1,tmp2/Sig_E2))
          invSig_kj <- invSigma[, , k12] + diag(c(X1X1[j, j] / Sig_E1, X2X2[j, j] / Sig_E2))
          post_Sig[k12, j, ,] <- solve(invSig_kj)
          mu_kj <- post_Sig[k12, j, ,] %*% c(tmp1 / Sig_E1, tmp2 / Sig_E2)
          mu1[k12, j] <- mu_kj[1]
          mu2[k12, j] <- mu_kj[2]
          
          gamma_k[j] <- log(prior_weights[j]) - log_det(invSig_kj) / 2 + t(mu_kj) %*% invSig_kj %*% mu_kj / 2
          
          post_mumu[[k12]][, , j] <- mu_kj %*% t(mu_kj) + post_Sig[k12, j, ,]
          
        }
        gamma_k <- gamma_k - max(gamma_k)
        gamma[k12,] <- exp(gamma_k) / sum(exp(gamma_k))
        # gamma[gamma > 0.999999999] <- 0.999999999
        # gamma[gamma < 0.000000001] <- 0.000000001
        
        sum_gammamu1 <- sum_gammamu1k + (mu1[k12,] * gamma[k12,])
        sum_gammamu2 <- sum_gammamu2k + (mu2[k12,] * gamma[k12,])
      }
      
      if (l == 1) {
        update_times[2] <- update_times[2] + 1
        k1 <- update_times[2]
        sum_gammamu1k <- sum_gammamu1 - (tmu1[k1,] * tgamma1[k1,])
        
        tgamma1_k <- rep(0, p)
        for (j in 1:p) {
          tmp1 <- X1y1[j] - sum(X1X1[, j] * sum_gammamu1k)
          post_Eta1[k1, j] <- 1 / (1 / Eta1[k1] + X1X1[j, j] / Sig_E1)
          tmu1[k1, j] <- post_Eta1[k1, j] * tmp1 / Sig_E1
          
          tgamma1_k[j] <- log(prior_weights1[j]) +
            log(post_Eta1[k1, j]) / 2 +
            tmu1[k1, j]^2 / post_Eta1[k1, j] / 2
          
          post_mumu1[k1, j] <- tmu1[k1, j]^2 + post_Eta1[k1, j]
          
        }
        tgamma1_k <- tgamma1_k - max(tgamma1_k)
        tgamma1[k1,] <- exp(tgamma1_k) / sum(exp(tgamma1_k))
        # gamma[gamma > 0.999999999] <- 0.999999999
        # gamma[gamma < 0.000000001] <- 0.000000001
        
        sum_gammamu1 <- sum_gammamu1k + (tmu1[k1,] * tgamma1[k1,])
      }
      
      if (l == 2) {
        update_times[3] <- update_times[3] + 1
        k2 <- update_times[3]
        sum_gammamu2k <- sum_gammamu2 - (tmu2[k2,] * tgamma2[k2,])
        
        tgamma2_k <- rep(0, p)
        for (j in 1:p) {
          tmp2 <- X2y2[j] - sum(X2X2[, j] * sum_gammamu2k)
          post_Eta2[k2, j] <- 1 / (1 / Eta2[k2] + X2X2[j, j] / Sig_E2)
          tmu2[k2, j] <- post_Eta2[k2, j] * tmp2 / Sig_E2
          
          tgamma2_k[j] <- log(prior_weights2[j]) +
            log(post_Eta2[k2, j]) / 2 +
            tmu2[k2, j]^2 / post_Eta2[k2, j] / 2
          
          post_mumu2[k2, j] <- tmu2[k2, j]^2 + post_Eta2[k2, j]
          
        }
        tgamma2_k <- tgamma2_k - max(tgamma2_k)
        tgamma2[k2,] <- exp(tgamma2_k) / sum(exp(tgamma2_k))
        # gamma[gamma > 0.999999999] <- 0.999999999
        # gamma[gamma < 0.000000001] <- 0.000000001
        
        sum_gammamu2 <- sum_gammamu2k + (tmu2[k2,] * tgamma2[k2,])
      }
      
    }
    
    
    # Update variance (M-step)
    if (estimate_prior_variance) {
      # shared effects
      for (k in 1:K) {
        Sigma[, , k] <- apply(post_mumu[[k]] * array(rep(gamma[k,], each = 4), c(2, 2, p)), c(1, 2), sum)
        invSigma <- array(apply(Sigma, 3, solve), c(2, 2, K))
      }
      
      # pop1-specific effects
      Eta1 <- rowSums(tgamma1 * post_mumu1)
      # pop2-specific effects
      Eta2 <- rowSums(tgamma2 * post_mumu2)
    }
    
    if (estimate_residual_variance) {
      post_bb1 <- t(sapply(post_mumu, function(x) x[1, 1,])) * gamma # K by p matrix of second moment of b1
      post_bb2 <- t(sapply(post_mumu, function(x) x[2, 2,])) * gamma # K by p matrix of second moment of b2
      post_aa1 <- post_mumu1 * tgamma1            # K by p matrix of second moment of a1
      post_aa2 <- post_mumu2 * tgamma2            # K by p matrix of second moment of a2
      
      Sig_E1 <- 1 / n1 * (sum((y1 - X1 %*% sum_gammamu1)^2) - sum((X1 %*% t(mu1 * gamma))^2) + sum(diag(X1X1) * colSums(post_bb1)) -
                            sum((X1 %*% t(tmu1 * tgamma1))^2) + sum(diag(X1X1) * colSums(post_aa1)))
      Sig_E2 <- 1 / n2 * (sum((y2 - X2 %*% sum_gammamu2)^2) - sum((X2 %*% t(mu2 * gamma))^2) + sum(diag(X2X2) * colSums(post_bb2)) -
                            sum((X2 %*% t(tmu2 * tgamma2))^2) + sum(diag(X2X2) * colSums(post_aa2)))
    }
    ELBO[1] <- get_ELBO(n1, n2, X1, X2, y1, y2, K, K1, K2, invSigma, Eta1, Eta2, Sig_E1, Sig_E2, prior_weights, prior_weights1, prior_weights2, X1X1, X2X2, gamma, tgamma1, tgamma2, mu1, mu2, tmu1, tmu2, post_Sig, post_Eta1, post_Eta2, sum_gammamu1, sum_gammamu2)
    cat("Initial ELBO value:", ELBO[1], ".\n")
  }
  
  ######################################################################
  
  
  # E-step for shared effects
  for (iter in 2:maxIter) {
    # update q (E-step)
    # XX <- matrix(0,2*p,2*p)
    # XX[1:p, 1:p] <- X1X1 / Sig_E1
    # XX[(p + 1):(2 * p), (p + 1):(2 * p)] <- X2X2 / Sig_E2
    #
    # invLambda <- XX + (solve(Omega) %x% diag(p))
    # cholinvLambda <- chol(invLambda)
    # Lambda <- chol2inv(cholinvLambda)
    # Lambdaj[1, 1] <- sum(diag(Lambda)[1:p])
    # Lambdaj[2, 2] <- sum(diag(Lambda)[(p + 1):(2 * p)])
    # Lambdaj[1, 2] <- Lambdaj[2, 1] <- sum(diag(Lambda[1:p, (p + 1):(2 * p)]))
    #
    # nu <- Lambda %*% c(1 / Sig_E1 * (X1y1 - X1X1 %*% sum_gammamu1), 1 / Sig_E2 * (X2y2 - X2X2 %*% sum_gammamu2))
    # nuj <- matrix(nu, p, 2)
    #
    # Xnu1 <- X1 %*% nu[1:p]
    # Xnu2 <- X2 %*% nu[(p + 1):(2 * p)]
    
    for (k in 1:K) {
      sum_gammamu1k <- sum_gammamu1 - (mu1[k,] * gamma[k,])
      sum_gammamu2k <- sum_gammamu2 - (mu2[k,] * gamma[k,])
      
      gamma_k <- rep(0, p)
      post_mumu[[k]] <- array(0, c(2, 2, p)) # posterior second monent of beta each component is 2*2*p matrix
      for (j in 1:p) {
        tmp1 <- X1y1[j] -
          sum(X1X1[, j] * sum_gammamu1k) -
          sum(X1[, j] * Xnu1)
        tmp2 <- X2y2[j] -
          sum(X2X2[, j] * sum_gammamu2k) -
          sum(X2[, j] * Xnu2)
        # invSig_kj <- invSigma[,,k]+diag(c(X1X1[j,j]/Sig_E1,X2X2[j,j]/Sig_E2))
        # mu_kj <- solve(invSig_kj,c(tmp1/Sig_E1,tmp2/Sig_E2))
        invSig_kj <- invSigma[, , k] + diag(c(X1X1[j, j] / Sig_E1, X2X2[j, j] / Sig_E2))
        post_Sig[k, j, ,] <- solve(invSig_kj)
        mu_kj <- post_Sig[k, j, ,] %*% c(tmp1 / Sig_E1, tmp2 / Sig_E2)
        mu1[k, j] <- mu_kj[1]
        mu2[k, j] <- mu_kj[2]
        
        gamma_k[j] <- log(prior_weights[j]) - log_det(invSig_kj) / 2 + t(mu_kj) %*% invSig_kj %*% mu_kj / 2
        
        post_mumu[[k]][, , j] <- mu_kj %*% t(mu_kj) + post_Sig[k, j, ,]
        
      }
      gamma_k <- gamma_k - max(gamma_k)
      gamma[k,] <- exp(gamma_k) / sum(exp(gamma_k))
      # gamma[gamma > 0.999999999] <- 0.999999999
      # gamma[gamma < 0.000000001] <- 0.000000001
      
      sum_gammamu1 <- sum_gammamu1k + (mu1[k,] * gamma[k,])
      sum_gammamu2 <- sum_gammamu2k + (mu2[k,] * gamma[k,])
      
      
      if (K1 > 0) {
        sum_gammamu1k <- sum_gammamu1 - (tmu1[k,] * tgamma1[k,])
        
        tgamma1_k <- rep(0, p)
        for (j in 1:p) {
          tmp1 <- X1y1[j] -
            sum(X1X1[, j] * sum_gammamu1k) -
            sum(X1[, j] * Xnu1)
          post_Eta1[k, j] <- 1 / (1 / Eta1[k] + X1X1[j, j] / Sig_E1)
          tmu1[k, j] <- post_Eta1[k, j] * tmp1 / Sig_E1
          
          tgamma1_k[j] <- log(prior_weights1[j]) +
            log(post_Eta1[k, j]) / 2 +
            tmu1[k, j]^2 / post_Eta1[k, j] / 2
          
          post_mumu1[k, j] <- tmu1[k, j]^2 + post_Eta1[k, j]
          
        }
        tgamma1_k <- tgamma1_k - max(tgamma1_k)
        tgamma1[k,] <- exp(tgamma1_k) / sum(exp(tgamma1_k))
        # gamma[gamma > 0.999999999] <- 0.999999999
        # gamma[gamma < 0.000000001] <- 0.000000001
        
        sum_gammamu1 <- sum_gammamu1k + (tmu1[k,] * tgamma1[k,])
      }
      
      
      if (K2 > 0) {
        sum_gammamu2k <- sum_gammamu2 - (tmu2[k,] * tgamma2[k,])
        
        tgamma2_k <- rep(0, p)
        for (j in 1:p) {
          tmp2 <- X2y2[j] -
            sum(X2X2[, j] * sum_gammamu2k) -
            sum(X2[, j] * Xnu2)
          post_Eta2[k, j] <- 1 / (1 / Eta2[k] + X2X2[j, j] / Sig_E2)
          tmu2[k, j] <- post_Eta2[k, j] * tmp2 / Sig_E2
          
          tgamma2_k[j] <- log(prior_weights2[j]) +
            log(post_Eta2[k, j]) / 2 +
            tmu2[k, j]^2 / post_Eta2[k, j] / 2
          
          post_mumu2[k, j] <- tmu2[k, j]^2 + post_Eta2[k, j]
          
        }
        tgamma2_k <- tgamma2_k - max(tgamma2_k)
        tgamma2[k,] <- exp(tgamma2_k) / sum(exp(tgamma2_k))
        # gamma[gamma > 0.999999999] <- 0.999999999
        # gamma[gamma < 0.000000001] <- 0.000000001
        
        sum_gammamu2 <- sum_gammamu2k + (tmu2[k,] * tgamma2[k,])
      }
      
    }
    XX <- matrix(0, 2 * p, 2 * p)
    XX[1:p, 1:p] <- X1X1 / Sig_E1
    XX[(p + 1):(2 * p), (p + 1):(2 * p)] <- X2X2 / Sig_E2
    
    invLambda <- XX + (solve(Omega) %x% diag(p))
    cholinvLambda <- chol(invLambda)
    Lambda <- chol2inv(cholinvLambda)
    Lambdaj[1, 1] <- sum(diag(Lambda)[1:p])
    Lambdaj[2, 2] <- sum(diag(Lambda)[(p + 1):(2 * p)])
    Lambdaj[1, 2] <- Lambdaj[2, 1] <- sum(diag(Lambda[1:p, (p + 1):(2 * p)]))
    
    nu <- Lambda %*% c(1 / Sig_E1 * (X1y1 - X1X1 %*% sum_gammamu1), 1 / Sig_E2 * (X2y2 - X2X2 %*% sum_gammamu2))
    nuj <- matrix(nu, p, 2)
    
    Xnu1 <- X1 %*% nu[1:p]
    Xnu2 <- X2 %*% nu[(p + 1):(2 * p)]
    
    ELBO[iter] <- get_ELBO_omega12(n1, n2, X1, X2, y1, y2, K, K1, K2, invSigma, Eta1, Eta2, Omega, Sig_E1, Sig_E2, prior_weights, prior_weights1, prior_weights2, X1X1, X2X2, gamma, tgamma1, tgamma2, mu1, mu2, tmu1, tmu2, post_Sig, post_Eta1, post_Eta2, sum_gammamu1, sum_gammamu2, nuj, Lambda, cholinvLambda)
    
    # Update variance (M-step)
    if (estimate_background_variance) {
      Omega <- (t(nuj) %*% nuj + Lambdaj) / p
    }
    
    if (estimate_prior_variance) {
      # shared effects
      for (k in 1:K) {
        Sigma[, , k] <- apply(post_mumu[[k]] * array(rep(gamma[k,], each = 4), c(2, 2, p)), c(1, 2), sum)
        invSigma <- array(apply(Sigma, 3, solve), c(2, 2, K))
      }
      
      if (K1 > 0) {
        # pop1-specific effects
        Eta1 <- rowSums(tgamma1 * post_mumu1)
      }
      
      if (K2 > 0) {
        # pop2-specific effects
        Eta2 <- rowSums(tgamma2 * post_mumu2)
      }
    }
    
    if (estimate_residual_variance) {
      post_bb1 <- t(sapply(post_mumu, function(x) x[1, 1,])) * gamma # K by p matrix of second moment of b1
      post_bb2 <- t(sapply(post_mumu, function(x) x[2, 2,])) * gamma # K by p matrix of second moment of b2
      
      Sig_E1 <- 1 / n1 * (sum((y1 - X1 %*% sum_gammamu1 - Xnu1)^2) - sum((X1 %*% t(mu1 * gamma))^2) +
                            sum(diag(X1X1) * colSums(post_bb1)) +
                            sum(X1X1 * Lambda[1:p, 1:p]))
      Sig_E2 <- 1 / n2 * (sum((y2 - X2 %*% sum_gammamu2 - Xnu2)^2) - sum((X2 %*% t(mu2 * gamma))^2) +
                            sum(diag(X2X2) * colSums(post_bb2)) +
                            sum(X2X2 * Lambda[(p + 1):(2 * p), (p + 1):(2 * p)]))
      
      
      if (K1 > 0) {
        post_aa1 <- post_mumu1 * tgamma1            # K by p matrix of second moment of a1
        Sig_E1 <- Sig_E1 - 1 / n1 * sum((X1 %*% t(tmu1 * tgamma1))^2) + 1 / n1 * sum(diag(X1X1) * colSums(post_aa1))
      }
      
      if (K2 > 0) {
        post_aa2 <- post_mumu2 * tgamma2            # K by p matrix of second moment of a2
        Sig_E2 <- Sig_E2 - 1 / n2 * sum((X2 %*% t(tmu2 * tgamma2))^2) + 1 / n2 * sum(diag(X2X2) * colSums(post_aa2))
      }
    }
    # Evaluate ELBO and check convergence
    # ELBO[iter] <- get_ELBO_omega12(n1, n2, X1, X2, y1, y2, K, K1, K2, invSigma, Eta1, Eta2, Omega, Sig_E1, Sig_E2, prior_weights, prior_weights1, prior_weights2, X1X1, X2X2, gamma, tgamma1, tgamma2, mu1, mu2, tmu1, tmu2, post_Sig, post_Eta1, post_Eta2, sum_gammamu1, sum_gammamu2, nuj, Lambda, cholinvLambda)
    cat(iter, "-th iteration finished.\t ELBO=", ELBO[iter], ",\t Diff=", ELBO[iter] - ELBO[iter - 1], "\n")
    if (ELBO[iter] - ELBO[iter - 1] < (tol * abs(ELBO[iter - 1]))) {
      ELBO <- ELBO[1:iter]
      break()
    }
  }
  
  ret <- list(gamma = gamma, tgamma1 = tgamma1, tgamma2 = tgamma2,
              mu1 = mu1, mu2 = mu2, tmu1 = tmu1, tmu2 = tmu2, nu = nu,
              Sigma = Sigma, Eta1 = Eta1, Eta2 = Eta2,
              Omega = Omega,
              Sig_E1 = Sig_E1, Sig_E2 = Sig_E2,
              ELBO = ELBO,
              K = K, K1 = K1, K2 = K2,
              prior_weights = prior_weights, prior_weights1 = prior_weights1, prior_weights2 = prior_weights2)
  return(ret)
}


get_ELBO_omega12 <- function(n1, n2, X1, X2, y1, y2, K, K1, K2, invSigma, Eta1, Eta2, Omega, Sig_E1, Sig_E2, prior_weights, prior_weights1, prior_weights2, X1X1, X2X2, gamma, tgamma1 = NULL, tgamma2 = NULL, mu1, mu2, tmu1 = NULL, tmu2 = NULL, post_Sig, post_Eta1 = NULL, post_Eta2 = NULL, sum_gammamu1, sum_gammamu2, nuj, Lambda, cholinvLambda) {
  p <- ncol(X1)
  llh_data <- -n1 * 0.5 * log(2 * pi * Sig_E1) -
    n2 * 0.5 * log(2 * pi * Sig_E2) -
    0.5 * sum((y1 - X1 %*% sum_gammamu1 - X1 %*% nuj[, 1])^2) / Sig_E1 -
    0.5 * sum((y2 - X2 %*% sum_gammamu2 - X2 %*% nuj[, 2])^2) / Sig_E2
  var_term <- 0
  beta_term <- 0
  q_beta <- 0
  for (k in 1:K) {
    for (j in 1:p) {
      # Sig_kj <- solve(invSigma[,,k]+diag(c(X1X1[j,j]/Sig_E1,X2X2[j,j]/Sig_E2)))
      mu_kj <- c(mu1[k, j], mu2[k, j])
      mumu <- mu_kj %*% t(mu_kj)
      var_term <- var_term -
        0.5 *
        gamma[k, j] *
        X1X1[j, j] *
        (mumu[1, 1] + post_Sig[k, j, 1, 1]) / Sig_E1 -
        0.5 *
        gamma[k, j] *
        X2X2[j, j] *
        (mumu[2, 2] + post_Sig[k, j, 2, 2]) / Sig_E2
      
      beta_term <- beta_term - 0.5 *
        gamma[k, j] *
        sum(diag(invSigma[, , k] %*% (mumu + post_Sig[k, j, ,])))
      
      q_beta <- q_beta + 0.5 * (gamma[k, j] * log_det(post_Sig[k, j, ,]) + gamma[k, j] * log_det(invSigma[, , k]))
    }
    var_term <- var_term +
      0.5 * sum((X1 %*% (mu1[k,] * gamma[k,]))^2) / Sig_E1 +
      0.5 * sum((X2 %*% (mu2[k,] * gamma[k,]))^2) / Sig_E2
  }
  
  if (K1 > 0) {
    for (k in 1:K1) {
      for (j in 1:p) {
        
        # Eta1_kj <- 1/(1/Eta1[k]+X1X1[j,j]/Sig_E1)
        
        var_term <- var_term -
          0.5 *
          tgamma1[k, j] *
          X1X1[j, j] *
          (tmu1[k, j]^2 + post_Eta1[k, j]) / Sig_E1
        
        beta_term <- beta_term -
          0.5 *
          tgamma1[k, j] *
          (tmu1[k, j]^2 + post_Eta1[k, j]) / Eta1[k]
        
        q_beta <- q_beta +
          0.5 * (tgamma1[k, j] * log(post_Eta1[k, j]) - tgamma1[k, j] * log(Eta1[k]))
      }
      var_term <- var_term + 0.5 * sum((X1 %*% (tmu1[k,] * tgamma1[k,]))^2) / Sig_E1
    }
  }
  
  
  if (K2 > 0) {
    for (k in 1:K2) {
      for (j in 1:p) {
        
        # Eta2_kj <- 1/(1/Eta2[k]+X2X2[j,j]/Sig_E2)
        
        var_term <- var_term -
          0.5 *
          tgamma2[k, j] *
          X2X2[j, j] *
          (tmu2[k, j]^2 + post_Eta2[k, j]) / Sig_E2
        
        beta_term <- beta_term -
          0.5 *
          tgamma2[k, j] *
          (tmu2[k, j]^2 + post_Eta2[k, j]) / Eta2[k]
        
        q_beta <- q_beta +
          0.5 * (tgamma2[k, j] * log(post_Eta2[k, j]) - tgamma2[k, j] * log(Eta2[k]))
      }
      var_term <- var_term + 0.5 * sum((X2 %*% (tmu2[k,] * tgamma2[k,]))^2) / Sig_E2
    }
  }
  
  
  gamma_term <- sum(t(gamma) * log(prior_weights))
  q_term <- -sum(gamma * log(gamma)) + q_beta
  
  if (K1 > 0) {
    gamma_term <- gamma_term + sum(t(tgamma1) * log(prior_weights1))
    q_term <- q_term - sum(tgamma1 * log(tgamma1))
  }
  
  if (K1 > 0) {
    gamma_term <- gamma_term + sum(t(tgamma2) * log(prior_weights2))
    q_term <- q_term - sum(tgamma2 * log(tgamma2))
  }
  
  XX <- matrix(0, 2 * p, 2 * p)
  XX[1:p, 1:p] <- X1X1 / Sig_E1
  XX[(p + 1):(2 * p), (p + 1):(2 * p)] <- X2X2 / Sig_E2
  
  bg_term <- -0.5 * sum(diag(t(nuj) %*% nuj %*% solve(Omega))) -
    0.5 * sum((XX + (solve(Omega) %x% diag(p))) * Lambda) -
    sum(log(diag(cholinvLambda))) -
    p * sum(log(diag(chol(Omega))))
  
  elbo <- llh_data +
    var_term +
    beta_term +
    gamma_term +
    q_term +
    bg_term
  
  return(elbo)
}

#
# get_ELBO <- function(n1,n2,X1,X2,y1,y2,K,invSigma,Sig_E1,Sig_E2,prior_weights,X1y1,X2y2,X1X1,X2X2,gamma,mu1,mu2,sum_gammamu1,sum_gammamu2){
#   # llh1 <- -n1*0.5*log(2*pi*Sig_E1) - n2*0.5*log(2*pi*Sig_E2) - 0.5*(sum((y1-X1%*%sum_gammamu1)^2)/Sig_E1+sum((y2-X2%*%sum_gammamu2)^2)/Sig_E2)
#
#   gammamu1 <- mu1*gamma
#   gammamu2 <- mu2*gamma
#   WRW1 <- gammamu1%*%X1X1%*%t(gammamu1)/Sig_E1
#   WRW2 <- gammamu2%*%X2X2%*%t(gammamu2)/Sig_E2
#   llh1 <- -n1*0.5*log(2*pi*Sig_E1) - n2*0.5*log(2*pi*Sig_E2) +
#     sum(X1y1*sum_gammamu1)/Sig_E1+sum(X2y2*sum_gammamu2)/Sig_E2-
#     0.5*sum(colSums(gamma*mu1^2)*diag(X1X1)/Sig_E1) - 0.5*sum(colSums(gamma*mu2^2)*diag(X2X2)/Sig_E2) -
#     0.5*(sum(WRW1)+sum(WRW2)-sum(diag(WRW1))-sum(diag(WRW2)))
#   var_term <- 0
#   trace_term <- matrix(0,nrow(gamma),ncol(gamma))
#   tmp1 <- tmp2 <- 0
#   for(k in 1:K){
#     for(j in 1:p){
#       Sig_kj <- solve(invSigma[,,k]+diag(c(X1X1[j,j]/Sig_E1,X2X2[j,j]/Sig_E2)))
#       mu <- c(mu1[k,j],mu2[k,j])
#       mumu <- mu%*%t(mu)
#       # var1 <- 0.5*X1X1[j,j]/Sig_E1*(gamma[k,j]*(Sig_kj[1,1]+mumu[1,1])-gamma[k,j]^2*mumu[1,1])
#       # var2 <- 0.5*X2X2[j,j]/Sig_E2*(gamma[k,j]*(Sig_kj[2,2]+mumu[2,2])-gamma[k,j]^2*mumu[2,2])
#       # var_term <- var_term + var1 + var2
#       trace_term[k,j] <- sum(X1X1[j,j]/Sig_E1*Sig_kj[1,1]) + sum(X2X2[j,j]/Sig_E2*Sig_kj[2,2])
#
#       tmp1 <- tmp1 + 0.5*gamma[k,j]*sum(diag(invSigma[,,k]%*%(Sig_kj+mumu-diag(2))))
#       tmp2 <- tmp2 + 0.5*gamma[k,j]*(log_det(Sig_kj)+log_det(invSigma[,,k])) - 0.5*log_det(invSigma[,,k])
#     }
#   }
#
#   llh2 <- -0.5*p*log_det(apply(invSigma,c(1,2),sum)) - tmp1
#   llh3 <- sum(t(gamma)*log(prior_weights))-sum(gamma*log(gamma)) + tmp2
#
#   # elbo <- llh1 + var_term + llh2 + llh3
#
#   elbo <- llh1 - 0.5*sum(gamma*trace_term) + llh2 + llh3
#   return(elbo)
# }
