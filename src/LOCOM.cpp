#include <RcppArmadillo.h>
#include <cmath>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::cube CalculateXX(const arma::mat & X){

    int K = X.n_cols;
    int n_sam = X.n_rows;
    arma::cube XX(K, K, n_sam);

    for (int s = 0; s < n_sam; s++){
        for (int i = 0; i < K; i++){
            XX(i, i, s) = X(s, i) * X(s, i);
            for (int j = i+1; j < K; j++){
                XX(i, j, s) = X(s, i) * X(s, j);
                XX(j, i, s) = XX(i, j, s);
            }
        }
    }
    return XX;
}


arma::cube UpdateXX(const arma::cube & XX, const arma::mat & X, arma::vec Yr){

    int n_sam = XX.n_slices;
    int K = XX.n_rows;
    arma::cube XX_perm = XX;

    for (int s = 0; s < n_sam; s++){
        XX_perm(0, 0, s) = Yr(s) * Yr(s);
        XX_perm(0, 1, s) = Yr(s);
        XX_perm(1, 0, s) = Yr(s);
        for (int j = 2; j < K; j ++){
            XX_perm(0, j, s) = Yr(s) * X(s, j);
            XX_perm(j, 0, s) = XX_perm(0, j, s);
        }
    }
    return XX_perm;
}


//' Solve Modified Estimation Equation using Quasi-Newton Algorithm
//'
//' This is a function to implement Quasi-Newton Algorithm to solve modified estimation equations.
//'
//'
//' @param freq_table frequency table
//' @param X Design matrix of all samples
//' @param beta_init Initial value of beta parameters
//' @param weight weight in the estimation equations
//' @param tol Stopping criterion. Default value is 1e-6.
//' @param iter_max Maximum iteration number

//' @examples
//' between(1:12, 7, 9)
//'
// [[Rcpp::export]]
arma::mat Newton(arma::mat freq_table, arma::mat X, arma::cube XX,
           arma::mat beta_init, arma::mat weight,
           double tol, int iter_max, double Firth_thresh,
           arma::vec prop_presence) {

    int n_otu = freq_table.n_cols;
    int n_sam = freq_table.n_rows;
    int K = X.n_cols;

    arma::mat beta = beta_init;
    arma::vec z;
    arma::vec u;
    arma::vec S;
    arma::vec J_temp;
    arma::vec J_temp_1;
    arma::mat J(K, K, fill::zeros);
    arma::mat J_inv;
    arma::vec H_temp;
    arma::mat H(K, K, fill::zeros);
    arma::vec step;

    for (int j = 0; j < n_otu; j++){

        beta.col(j) = beta_init.col(j);

        for (int i = 0; i < iter_max; i++){

            z = exp(X * beta.col(j));
            u = z / (1 + z);
            S = X.t() * ( (freq_table.col(j) - u) % (weight.col(j)) ); //*:matrix multiplication; %: element-wise multiplication

            J_temp = - (u % (1 - u)) % (weight.col(j));

            J.fill(0);
            for (int i_sam = 0; i_sam < n_sam; i_sam++){
                J = J + J_temp(i_sam) * XX.slice(i_sam);
            }
            J_inv = inv(J);

            //Firth

            if (prop_presence(j) < Firth_thresh){

                J_temp_1 = J_temp % (1 - 2*u);

                for (int k = 0; k < K; k++) {

                    H_temp = J_temp_1 % X.col(k);

                    H.fill(0);
                    for (int i_sam = 0; i_sam < n_sam; i_sam++){
                        H = H + H_temp(i_sam) * XX.slice(i_sam);
                    }
                    S(k) = S(k) + 0.5*accu(H % J_inv); //trace(H * J_inv);
                }
            }

            // update

            step = J_inv * S;

            if (sum(abs(step)) > 5*K){
                beta.col(j).fill(0);
            } else {
                beta.col(j) = beta.col(j) - step;
            }

            if (sum(abs(step) < tol) == K){
                break;
            }
        } // iter
    } // otu
    return beta;
} // Newton()


// [[Rcpp::export]]
arma::mat perm_Newton(arma::mat freq_table, arma::vec Yr, arma::mat X, arma::cube XX,
                         arma::mat beta_init, arma::mat weight,
                         arma::umat perm,
                         double tol, int iter_max, double Firth_thresh,
                         arma::vec prop_presence) {

    int n_otu = freq_table.n_cols;
    int n_perm = perm.n_cols;

    arma::mat beta_est_temp(n_perm, n_otu);
    arma::mat var_est_temp(n_perm, n_otu);

    arma::vec Yr_perm;
    arma::mat X_perm = X;
    uvec o;

    for (int i_perm = 0; i_perm < n_perm; i_perm++){

        o = perm.col(i_perm);

        Yr_perm = Yr(o);
        X_perm.col(0) = Yr_perm;
        XX = UpdateXX(XX, X, Yr_perm);

        mat res = Newton(freq_table, X_perm, XX, beta_init, weight, tol, iter_max, Firth_thresh, prop_presence);

        beta_est_temp.row(i_perm) = res.row(0);
    }

    return beta_est_temp;

} // perm_Newton()

