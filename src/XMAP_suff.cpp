// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "XMAP_suff.hpp"
#include <cstdlib>
#include <sys/time.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::List XMAP_suff(const arma::cube & XX, const arma::mat & Xy, const arma::colvec yty, const arma::colvec & n,  
                     int K, arma::cube Sigma, arma::mat Omega, arma::colvec & Sig_E,
                     const arma::colvec & prior_weights, int maxIter, double tol,
                     bool estimate_prior_variance, bool estimate_residual_variance, bool estimate_background_variance) {
  Timer timer;
  int p = XX.n_cols;
  int t = Sigma.n_cols;
  
  arma::cube invSigma(t,t,K);
  for(int k = 0; k < K; k++) {
    invSigma.slice(k) = Sigma.slice(k).i();
  }
  
  // posterior foreground parameters
  arma::mat gamma(K,p,arma::fill::value(1.0/p));
  arma::cube mu(t,K,p);
  arma::field<arma::cube> post_Sig(K);
  
  arma::mat sum_gammamu(t,p);
  for(int j = 0; j < t; j++) {
    sum_gammamu.row(j) = arma::sum(gamma%arma::mat(mu.row(j)),0);
  }
  
  arma::field<arma::cube> post_mumu(K);
  arma::mat XXsum_gammamu(p,t);
  for(int tt = 0; tt < t; tt++){
    XXsum_gammamu.col(tt) = XX.slice(tt) * sum_gammamu.row(tt).t();
  }
  // posterior background
  // arma::mat XXtp(t*p,t*p,arma::fill::value(0));
  // for(int tt = 0; tt < t; tt++){
  //   XXtp.submat(p*tt,p*tt,p*(tt+1)-1,p*(tt+1)-1) = XX.slice(tt) / Sig_E(tt);
  // }
  // arma::mat invLambda = XXtp + arma::kron(Omega.i(),arma::eye(p,p));
  // arma::mat cholinvLambda = arma::chol(invLambda);
  // arma::mat invU = cholinvLambda.i();
  // arma::mat Lambda = invU * invU.t();
  // 
  // arma::colvec nutp = Lambda * arma::vectorise(Xy.each_row() / Sig_E.t());
  // arma::mat nu = arma::reshape(nutp,p,t).t();
  // arma::mat XXnu(p,t);
  // for(int tt = 0; tt < t; tt++){
  //   XXnu.col(tt) = XX.slice(tt) * nu.row(tt).t();
  // }
  
  
  
  
  
  
  arma::mat nu(t,p);
  arma::mat XXnu(p,t);
  // for(int j = 0; j < t; j++){
  //   XXnu.col(j) = XX.slice(j) * nu.row(j).t();
  // }
  
  arma::colvec ELBO(maxIter);
  ELBO(0) = - arma::datum::inf;
  // std::cout << "Time for initiallization is " << timer.update_time() << " sec" << std::endl;
  
  for(int iter = 1; iter < maxIter; iter++){
    
    // E-step
    // posterior of background
    arma::mat XXtp(t*p,t*p,arma::fill::value(0));
    for(int tt = 0; tt < t; tt++){
      XXtp.submat(p*tt,p*tt,p*(tt+1)-1,p*(tt+1)-1) = XX.slice(tt) / Sig_E(tt);
    }
    arma::mat invLambda = XXtp + arma::kron(Omega.i(),arma::eye(p,p));
    arma::mat cholinvLambda = arma::chol(invLambda);
    arma::mat invU = cholinvLambda.i();
    arma::mat Lambda = invU * invU.t();

    arma::mat Lambdaj(t,t);
    for(int tt = 0; tt < t; tt++){
      Lambdaj(tt,tt) = sum(Lambda.submat(p*tt,p*tt,p*(tt+1)-1,p*(tt+1)-1).diag());
      for(int ttt = tt+1; ttt < t; ttt++){
        Lambdaj(tt,ttt) = sum(Lambda.submat(p*tt,p*ttt,p*(tt+1)-1,p*(ttt+1)-1).diag());
      }
    }
    Lambdaj = arma::symmatu(Lambdaj);
    // std::cout << "Time for getting Lambdaj is " << timer.update_time() << " sec" << std::endl;
    
    arma::mat tmp = Xy - XXsum_gammamu;
    arma::colvec nutp = Lambda * arma::vectorise(tmp.each_row() / Sig_E.t());
    nu = arma::reshape(nutp,p,t).t();
    for(int tt = 0; tt < t; tt++){
      XXnu.col(tt) = XX.slice(tt) * nu.row(tt).t();
    }
    
    
    
    // posterior of foreground
    for(int k = 0; k < K; k++){
      arma::mat sum_gammamu_k = sum_gammamu - (arma::mat(mu.col(k)).each_row() % gamma.row(k));
      
      // arma::colvec gamma_k(p);
      post_mumu(k) = arma::cube(t,t,p);
      post_Sig(k) = arma::cube(t,t,p);
      arma::colvec tmp;
      
      for(int j = 0; j < p; j++){
        arma::mat invSig_kj = invSigma.slice(k);
        invSig_kj.diag() += arma::colvec(XX.tube(j,j))/Sig_E;
        // for(int tt = 0; tt < t; tt++){
        //   invSig_kj(tt,tt) += XX(j,j,tt)/Sig_E(tt);
        // }
        post_Sig(k).slice(j) = invSig_kj.i();
        
        tmp = (Xy.row(j) - sum(arma::mat(XX.col(j)) % sum_gammamu_k.t(), 0) - XXnu.row(j)).t() / Sig_E;
        mu.slice(j).col(k) = post_Sig(k).slice(j) * tmp;
        
        post_mumu(k).slice(j) = mu.slice(j).col(k) * mu.slice(j).col(k).t() + post_Sig(k).slice(j);
        gamma(k,j) = log(prior_weights(j)) - arma::log_det_sympd(invSig_kj) / 2.0 + arma::as_scalar(mu.slice(j).col(k).t() * invSig_kj * mu.slice(j).col(k)) / 2.0;
      }
      
      gamma.row(k) -= gamma.row(k).max();
      gamma.row(k) = exp(gamma.row(k)) / sum(exp(gamma.row(k)));
      
      sum_gammamu = sum_gammamu_k + (arma::mat(mu.col(k)).each_row() % gamma.row(k));
    }
    
    // arma::mat XXsum_gammamu(p,t);
    for(int tt = 0; tt < t; tt++){
      XXsum_gammamu.col(tt) = XX.slice(tt) * sum_gammamu.row(tt).t();
    }

    ELBO(iter) = get_ELBO_XMAP_suff(XX, Xy, yty, XXtp, n, K, invSigma, Omega, Sig_E, prior_weights, gamma, mu, post_Sig, sum_gammamu, 
         XXsum_gammamu, nu, XXnu, Lambda, cholinvLambda);
    // std::cout << "Time for evaluating ELBO is " << timer.update_time() << " sec" << std::endl;
    
    // M-step
    if(estimate_background_variance) {
      Omega = (nu * nu.t() + Lambdaj) / p;
      // std::cout << "Time for estimating Omega is " << timer.update_time() << " sec" << std::endl;
    }
    
    if(estimate_prior_variance) {
      for(int k = 0; k < K; k++) {
        arma::mat tmp(t,t,arma::fill::value(0));
        for(int j = 0; j < p; j++){
          tmp += post_mumu(k).slice(j) * gamma(k,j);
        }
        Sigma.slice(k) = tmp;
        invSigma.slice(k) = Sigma.slice(k).i();
      }
      // std::cout << "Time for estimating Sigma is " << timer.update_time() << " sec" << std::endl;
    }
    
    if(estimate_residual_variance){
      for(int tt = 0; tt < t; tt++){
        arma::mat post_bb(K,p);
        for(int k = 0; k < K; k++){
          post_bb.row(k) = arma::rowvec(post_mumu(k).tube(tt,tt)) % gamma.row(k);
        }
        
        arma::mat tmp = arma::mat(mu.row(tt)) % gamma;
        Sig_E(tt) = 1/n(tt) * (yty(tt) + 
          sum(sum_gammamu.row(tt) % XXsum_gammamu.col(tt).t()) + 
          sum(nu.row(tt) % XXnu.col(tt).t()) -
          2 * sum(Xy.col(tt).t() % sum_gammamu.row(tt)) -
          2 * sum(Xy.col(tt).t() % nu.row(tt)) +
          2 * sum(nu.row(tt) % XXsum_gammamu.col(tt).t()) -
          accu(tmp.t() % (XX.slice(tt) * tmp.t())) + 
          sum(XX.slice(tt).diag() % sum(post_bb,0).t()) + 
          accu(XX.slice(tt) % Lambda.submat(p*tt,p*tt,p*(tt+1)-1,p*(tt+1)-1)));
        
        
        // std::cout << nu.submat(0,0,1,5) << std::endl;
        // std::cout << yty(tt) << std::endl;
        // std::cout << sum(sum_gammamu.row(tt) % XXsum_gammamu.col(tt).t()) << std::endl;
        // std::cout << sum(nu.row(tt) % XXnu.col(tt).t()) << std::endl;
        // std::cout << 2 * sum(Xy.col(tt).t() % sum_gammamu.row(tt))<< std::endl;
        // std::cout << 2 * sum(Xy.col(tt).t() % nu.row(tt)) << std::endl;
        // std::cout << 2 * sum(nu.row(tt) % XXsum_gammamu.col(tt).t()) << std::endl;
        // std::cout << accu(tmp.t() % (XX.slice(tt) * tmp.t())) << std::endl;
        // std::cout << sum(XX.slice(tt).diag() % sum(post_bb,0).t()) << std::endl;
        // std::cout << accu(XX.slice(tt) % Lambda.submat(p*tt,p*tt,p*(tt+1)-1,p*(tt+1)-1)) << std::endl;
        // std::cout << Sig_E(tt) << std::endl;
      }
      // std::cout << "Time for estimating Sigma_E is " << timer.update_time() << " sec" << std::endl;
    }
    
    // ELBO(iter) = get_ELBO_XMAP_suff(XX, Xy, yty, XXtp, n, K, invSigma, Omega, Sig_E, prior_weights, gamma, mu, post_Sig, sum_gammamu, 
    //                                  XXsum_gammamu, nu, XXnu, Lambda, cholinvLambda);
    printf("%d-th iteration finished.\t ELBO = %f,\t Diff = %f \n", iter+1, ELBO(iter), ELBO(iter)-ELBO(iter-1));
    if(ELBO(iter)-ELBO(iter-1) < tol * abs(ELBO(iter-1))){
      ELBO = ELBO.subvec(0,iter);
      break;
    }
    // printf("%d-th iteration \n", iter+1);
    
  }
  
  return Rcpp::List::create(Rcpp::Named("gamma")=gamma,
                            Rcpp::Named("mu")=mu,
                            Rcpp::Named("nu")=nu,
                            Rcpp::Named("Sigma")=Sigma,
                            Rcpp::Named("Omega")=Omega,
                            Rcpp::Named("Sig_E")=Sig_E,
                            Rcpp::Named("ELBO")=ELBO,
                            Rcpp::Named("K")=K,
                            Rcpp::Named("prior_weights")=prior_weights);
}

