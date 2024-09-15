#ifndef XMAP_suff_hpp
#define XMAP_suff_hpp

double get_ELBO_XMAP_suff(const arma::cube & XX, const arma::mat & Xy, const arma::colvec yty, arma::mat & XXtp,
                          const arma::colvec & n, int K, arma::cube & invSigma, arma::mat Omega, arma::colvec & Sig_E,
                          const arma::colvec & prior_weights,
                          arma::mat gamma, arma::cube mu, arma::field<arma::cube> post_Sig, arma::mat sum_gammamu, arma::mat XXsum_gammamu, 
                          arma::mat nu, arma::mat XXnu, arma::mat Lambda, arma::mat cholinvLambda); 
double get_ELBO_XMAP_suff_old(const arma::cube & XX, const arma::mat & Xy, const arma::colvec yty, arma::mat & XXtp,
                          const arma::colvec & n, int K, arma::cube & invSigma, arma::mat Omega, arma::colvec & Sig_E,
                          const arma::colvec & prior_weights,
                          arma::mat gamma, arma::cube mu, arma::field<arma::cube> post_Sig, arma::mat sum_gammamu, arma::mat XXsum_gammamu, 
                          arma::mat nu, arma::mat XXnu, arma::mat Lambda, arma::mat cholinvLambda); 


class Timer {
private:
  double prevtime, curtime;
  
public:
  
  Timer(void);
  double update_time(void);
  double get_time(void);
};

#endif /* XMAP_suff_hpp */
