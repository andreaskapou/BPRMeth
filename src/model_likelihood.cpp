// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' @rdname model_log_likelihood
//'
//' @export
// [[Rcpp::export]]
double bpr_log_likelihood(const arma::vec& w, const arma::mat& X,
                          const arma::mat& H, const double lambda,
                          const bool is_nll){
    int nrows = X.n_rows;
    int ncols = X.n_cols;
    double ll = 0.0;

    // Predictions of the target variables
    Rcpp::NumericVector g = Rcpp::wrap(H * w);
    // Compute the cdf of N(0,1) distribution (i.e. probit function)
    Rcpp::NumericVector Phi = Rcpp::pnorm(g);

    for (int i = 0; i < nrows; i++){
        // In extreme cases where probit is 0 or 1, subtract a tiny number
        // so we can evaluate the log(0) when computing the likelihood
        if (Phi[i] > (1 - 1e-15)){ Phi[i] = 1 - 1e-15;
        }else if (Phi[i] < 1e-15){ Phi[i] = 1e-15; }
        // Compute the log likelihood
        if (ncols == 3){ // If Binomial distributed data
            ll += R::dbinom(X(i, 2), X(i, 1), Phi[i], true);
        }else{           // If Bernoulli distributed data
            ll += R::dbinom(X(i, 1), 1, Phi[i], true);
        }
    }

    // Compute ridge regression likelihood
    ll = ll - lambda * arma::as_scalar(w.t() * w);
    // If we require the Negative Log Likelihood
    if (is_nll == true){ ll = (-1) * ll; }
    return ll;
}


//' @rdname model_log_likelihood
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector bpr_gradient(const arma::vec& w, const arma::mat& X,
                                 const arma::mat& H, const double lambda,
                                 const bool is_nll){

    int nrows = X.n_rows;
    int ncols = X.n_cols;
    int M = w.size();

    // Predictions of the target variables
    Rcpp::NumericVector g = Rcpp::wrap(H * w);
    // Compute the cdf of N(0,1) distribution (i.e. probit function)
    Rcpp::NumericVector Phi = Rcpp::pnorm(g);
    // Compute the density of a N(0,1) distribution
    Rcpp::NumericVector N = Rcpp::dnorm(g);
    Rcpp::NumericVector gr(M);

    for (int i = 0; i < nrows; i++){
        // In extreme cases where probit is 0 or 1, subtract a tiny number
        // so we can evaluate the log(0) when computing the likelihood
        if (Phi[i] > (1 - 1e-15)){ Phi[i] = 1 - 1e-15;
        }else if (Phi[i] < 1e-15){ Phi[i] = 1e-15; }
        if (N[i] < 1e-15){ N[i] = 1e-15; }
        // Compute the gradient vector w.r.t the coefficients w
        if (ncols == 3){  // If Binomial distributed data
            for (int m = 0; m < M; m++){
                gr[m] += N[i] * H(i, m) * (X(i, 2) / Phi[i] - (X(i, 1) -
                    X(i, 2)) / (1 - Phi[i]) );
            }
        }else{            // If Bernoulli distributed data
            for (int m = 0; m < M; m++){
                gr[m] += N[i] * H(i, m) * (X(i, 1) / Phi[i] -
                    (1 - X(i, 1)) / (1 - Phi[i]) );
            }
        }
    }
    for (int m = 0; m < M; m++){
        // Compute ridge regression likelihood
        gr[m] -= 2 * lambda * w[m];
        // If we require the Negative Log Likelihood
        if (is_nll == true){ gr[m] *= -1; }
    }
    return gr;
}


//' @rdname model_log_likelihood
//'
//' @export
// [[Rcpp::export]]
double betareg_log_likelihood(const arma::vec& w, arma::mat& X,
                              const arma::mat& H, const double lambda,
                              const bool is_nll){
    int nrows = X.n_rows;
    double ll = 0.0;

    // Predictions of the target variables
    Rcpp::NumericVector tmp = Rcpp::wrap(H * w);
    // Compute the cdf of N(0,1) distribution (i.e. probit function)
    Rcpp::NumericVector Phi = Rcpp::pnorm(tmp);

    for (int i = 0; i < nrows; i++){
        // In extreme cases where probit is 0 or 1, subtract a tiny number
        // so we can evaluate the log(0) when computing the likelihood
        if (Phi[i] > (1 - 1e-15)){ Phi[i] = 1 - 1e-15;
        }else if (Phi[i] < 1e-15){ Phi[i] = 1e-15; }
        // Do the same for the actual observations
        if (X(i,1) > (1 - 1e-15)){ X(i,1) = 1 - 1e-15;
        }else if (X(i,1) < 1e-15){ X(i,1) = 1e-15; }

        ll += R::lgammafn(X(i,2)) - R::lgammafn(Phi[i]*X(i,2)) -
            R::lgammafn((1-Phi[i])*X(i,2)) + (Phi[i]*X(i,2)-1)*log(X(i,1)) +
            ((1-Phi[i])*X(i,2)-1)*log(1-X(i,1));
    }

    // Compute ridge regression likelihood
    ll = ll - lambda * arma::as_scalar(w.t() * w);
    // If we require the Negative Log Likelihood
    if (is_nll == true){ ll = (-1) * ll; }
    return ll;
}



//' @rdname model_log_likelihood
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector betareg_gradient(const arma::vec& w, arma::mat& X,
                                     const arma::mat& H, const double lambda,
                                     const bool is_nll){
    int nrows = X.n_rows;
    int M = w.size();

    // Predictions of the target variables
    Rcpp::NumericVector g = Rcpp::wrap(H * w);
    // Compute the cdf of N(0,1) distribution (i.e. probit function)
    Rcpp::NumericVector Phi = Rcpp::pnorm(g);
    // Compute the density of a N(0,1) distribution
    Rcpp::NumericVector N = Rcpp::dnorm(g);
    Rcpp::NumericVector gr(M);

    for (int i = 0; i < nrows; i++){
        // In extreme cases where probit is 0 or 1, subtract a tiny number
        // so we can evaluate the log(0) when computing the likelihood
        if (Phi[i] > (1 - 1e-15)){ Phi[i] = 1 - 1e-15;
        }else if (Phi[i] < 1e-15){ Phi[i] = 1e-15; }
        // Do the same for the actual observations
        if (X(i,1) > (1 - 1e-15)){ X(i,1) = 1 - 1e-15;
        }else if (X(i,1) < 1e-15){ X(i,1) = 1e-15; }
        if (N[i] < 1e-15){ N[i] = 1e-15; }

        // Compute the gradient vector w.r.t the coefficients w
        for (int m = 0; m < M; m++){
            gr[m] += N[i] * X(i,2) * ( log(X(i,1)) - log(1-X(i,1)) -
                R::digamma(Phi[i]*X(i,2)) + R::digamma((1-Phi[i])*X(i,2)) ) * H(i,m);
        }
    }
    for (int m = 0; m < M; m++){
        // Compute ridge regression likelihood
        gr[m] -= 2 * lambda * w[m];
        // If we require the Negative Log Likelihood
        if (is_nll == true){ gr[m] *= -1; }
    }
    return gr;
}


//' @rdname model_log_likelihood
//'
//' @export
// [[Rcpp::export]]
double sum_weighted_bpr_lik(const arma::vec& w, const Rcpp::List& X_list,
                            const Rcpp::List& H_list,
                            const arma::vec& r_nk, const double lambda,
                            const bool is_nll){
    // Number of regions
    int N = X_list.size();
    Rcpp::NumericVector res(N);
    for (int i = 0; i < N; i++){
        // Extract observations in each region
        arma::mat X = X_list[i];
        // Extract deisgn matrix of each region
        arma::mat H = H_list[i];
        // Compute BPR likelihood
        res[i] = bpr_log_likelihood(w, X, H, lambda, is_nll);
    }
    // Inner product with the weight vector of posterior probabilities
    arma::vec ll = as<arma::vec>(res);
    return arma::as_scalar(r_nk.t() * ll);
}


//' @rdname model_log_likelihood
//'
//' @export
// [[Rcpp::export]]
arma::rowvec sum_weighted_bpr_grad(const arma::vec& w, const Rcpp::List& X_list,
                                   const Rcpp::List& H_list,
                                   const arma::vec& r_nk,
                                   const double lambda, const bool is_nll){
    // Number of regions
    int N = X_list.size();
    // Number of basis functions
    int M = w.size();
    Rcpp::NumericMatrix res(N, M);
    for (int i = 0; i < N; i++){
        // Extract observations in each region
        arma::mat X = X_list[i];
        // Extract deisgn matrix of each region
        arma::mat H = H_list[i];
        // Compute the gradient of BPR model
        res(i, _) = bpr_gradient(w, X, H, lambda, is_nll);
    }
    // Inner product with the weight vector of posterior probabilities
    arma::mat ll = as<arma::mat>(res);
    arma::rowvec w_lik = r_nk.t() * ll;
    return w_lik;
}


//' @rdname model_log_likelihood
//'
//' @export
// [[Rcpp::export]]
double sum_weighted_betareg_lik(const arma::vec& w, const Rcpp::List& X_list,
                                const Rcpp::List& H_list,
                                const arma::vec& r_nk, const double lambda,
                                const bool is_nll){
    // Number of regions
    int N = X_list.size();
    Rcpp::NumericVector res(N);
    for (int i = 0; i < N; i++){
        // Extract observations in each region
        arma::mat X = X_list[i];
        // Extract deisgn matrix of each region
        arma::mat H = H_list[i];
        // Compute Beta regression likelihood
        res[i] = betareg_log_likelihood(w, X, H, lambda, is_nll);
    }
    // Inner product with the weight vector of posterior probabilities
    arma::vec ll = as<arma::vec>(res);
    return arma::as_scalar(r_nk.t() * ll);
}


//' @rdname model_log_likelihood
//'
//' @export
// [[Rcpp::export]]
arma::rowvec sum_weighted_betareg_grad(const arma::vec& w, const Rcpp::List& X_list,
                                       const Rcpp::List& H_list,
                                       const arma::vec& r_nk,
                                       const double lambda, const bool is_nll){
    // Number of regions
    int N = X_list.size();
    // Number of basis functions
    int M = w.size();
    Rcpp::NumericMatrix res(N, M);
    for (int i = 0; i < N; i++){
        // Extract observations in each region
        arma::mat X = X_list[i];
        // Extract deisgn matrix of each region
        arma::mat H = H_list[i];
        // Compute the gradient of Beta regression model
        res(i, _) = betareg_gradient(w, X, H, lambda, is_nll);
    }
    // Inner product with the weight vector of posterior probabilities
    arma::mat ll = as<arma::mat>(res);
    arma::rowvec w_lik = r_nk.t() * ll;
    return w_lik;
}



// @rdname model_log_likelihood
//
// @export
// TODO [[Rcpp::export]]
Rcpp::NumericVector bpr_lik_region(const arma::vec& w, const Rcpp::List& X_list,
                                   const Rcpp::List& H_list,
                                   const double lambda, const bool is_nll){
    // Number of regions
    int N = X_list.size();
    Rcpp::NumericVector res(N);
    for (int i = 0; i < N; i++){
        // Extract observations in each region
        arma::mat X = X_list[i];
        // Extract deisgn matrix of each region
        arma::mat H = H_list[i];
        //         Rcpp::List des_obj = H_list[i];
        //         arma::mat H = as<arma::mat>(des_obj["H"]);
        // Compute BPR likelihood
        res[i] = bpr_log_likelihood(w, X, H, lambda, is_nll);
    }
    return res;
}


// @rdname model_log_likelihood
//
// @export
// TODO [[Rcpp::export]]
Rcpp::NumericMatrix bpr_lik_resp(const arma::mat& w, const Rcpp::List& X_list,
                                 const Rcpp::List& H_list, const arma::vec pi_k,
                                 const double lambda, const bool is_nll){
    // Number of regions
    int N = X_list.size();
    int K = w.n_cols;
    int k;
    Rcpp::NumericMatrix res(N, K);
    for (int i = 0; i < N; i++){
        // Extract observations in each region
        arma::mat X = X_list[i];
        // Extract deisgn matrix of each region
        arma::mat H = H_list[i];
        // Compute BPR likelihood
        for (k = 0; k < K; k++){
            arma::vec ww = w.col(k);
            res(i, k) = pi_k[k] + bpr_log_likelihood(ww, X, H, lambda, is_nll);
        }
    }
    return res;
}
