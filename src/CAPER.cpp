//******************************************************************************//
// Date: 2026-3-11
// Author: Shiquan Sun; Ye Li; Jin Ning 
// Condition-associated pattern extraction and recovery in multi-condition single-cell transcriptomics with CAPER
//******************************************************************************//
#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ctime>

using namespace std;
using namespace arma;
using namespace Rcpp;

#define ARMA_DONT_PRINT_ERRORS

// diag(I) where I = diag(I_{k1}, sigma_c I_{k2}, sigma_s I_{k2})
arma::vec Itilde(const int k1, const int k2, const double sigma_c, const double sigma_s){
  arma::vec diag_I = arma::zeros<arma::vec>(k1 + 2*k2);
  diag_I.subvec(0, k1-1).fill(1.0);
  diag_I.subvec(k1, k1+k2-1).fill(sigma_c);
  diag_I.subvec(k1+k2, k1+2*k2-1).fill(sigma_s);
  return diag_I;
}// end func

// compute I * Lambda^T by row-scaling Lambda^T according to diag(I)
arma::mat UpdateItildeLambdaT(const arma::mat& Lambda,
                             const int k1, const int k2,
                             const double sigma_c, const double sigma_s){
  arma::mat LambdaT = Lambda.t();
  arma::mat I_LambdaT = LambdaT;

  I_LambdaT.rows(0, k1-1) = LambdaT.rows(0, k1-1);
  I_LambdaT.rows(k1, k1+k2-1) = sigma_c * LambdaT.rows(k1, k1+k2-1);
  I_LambdaT.rows(k1+k2, k1+2*k2-1) = sigma_s * LambdaT.rows(k1+k2, k1+2*k2-1);

  return I_LambdaT;
}// end func

// (Psi + Lambda * I * Lambda^T)^(-1) via Woodbury (Psi is diagonal stored as psi_diag)
arma::mat FastInverseMatV2(const arma::vec& psi_diag,
                           const arma::mat& Lambda,
                           const arma::mat& I_LambdaT) {
  arma::vec Psi_inv_diag = 1.0/(psi_diag + 1e-6);
  arma::mat PsiInv_Lambda = Lambda.each_col() % Psi_inv_diag;
  arma::mat I_LambdaT_PsiInv = I_LambdaT.each_row() % Psi_inv_diag.t();

  arma::mat Kmat = I_LambdaT * PsiInv_Lambda;
  Kmat.diag() += 1.0;

  arma::mat solve_term = solve(Kmat, I_LambdaT_PsiInv, solve_opts::fast);
  arma::mat inv_mat = -PsiInv_Lambda * solve_term;
  inv_mat.diag() += Psi_inv_diag;

  return inv_mat;
}// end func

// E-step: posterior moments for Z_g = [z_g; z_g^c; z_g^s]
// IMPORTANT: keep the exact multiplication structure consistent with the original implementation.
void ComputeExpVar(const arma::vec& psi_diag,
                   const arma::mat& Lambda,
                   const arma::mat& I_LambdaT,
                   const arma::vec& diag_I,
                   const arma::mat& X_all,
                   arma::mat &E_Z,
                   arma::mat &Var_Z,
                   arma::mat &E_Z_part){

    arma::mat shared_term = FastInverseMatV2(psi_diag, Lambda, I_LambdaT);

    // E_Z_part = I Lambda^T (Psi + Lambda I Lambda^T)^(-1)
    E_Z_part = I_LambdaT * shared_term;

    // E[Z | X] = E_Z_part X^T
    E_Z = E_Z_part * X_all.t();

    // Var(Z | X) = I - I Lambda^T (Psi + Lambda I Lambda^T)^(-1) Lambda I
    // Same structure as the original: -E_part * (I_LambdaT)^T + diag(I)
    Var_Z = -E_Z_part * I_LambdaT.t();
    Var_Z.diag() += diag_I;
}// end func

arma::mat ComputeSVD(const arma::mat X, const int k, const bool flag_sparse){
    arma::mat U;
    arma::vec s;
    arma::mat V;
    if(flag_sparse == TRUE){
        arma::sp_mat sp_X(X);
        svds(U, s, V, sp_X, k);
    }else{
        svd(U, s, V, X);
    }
    arma::mat D = arma::zeros<arma::mat>(k, k);
    D.diag() = sqrt(s(arma::span(0, k-1)));
    arma::mat Z = U.cols(0, k-1) * D;
    arma::mat Z_norm = (Z - mean(vectorise(Z))) / stddev(vectorise(Z));
    return Z_norm;
}// end func

// initialize Lambda = [Lambda^c, Lambda_u^c, 0; Lambda^s, 0, Lambda_u^s]
arma::mat ComputeLambda(arma::mat X_c, arma::mat X_s, arma::mat &Lambda,
                        const int k1, const int k2,
                        const int n_c, const int n_s){

    arma::mat Z_c = ComputeSVD(X_c, k1, FALSE);
    arma::mat Zc_ZcT_inv_Zc = Z_c * inv(Z_c.t() * Z_c);

    arma::mat Lambda_c = X_c.t() * Zc_ZcT_inv_Zc;
    arma::mat Lambda_s = X_s.t() * Zc_ZcT_inv_Zc;

    arma::mat R_c = X_c - Z_c * Lambda_c.t();

    arma::mat Z_u_c = ComputeSVD(R_c, k2, FALSE);
    arma::mat Lambda_u_c = R_c.t() * Z_u_c * inv(Z_u_c.t() * Z_u_c);

    arma::mat Z_s = ComputeSVD(X_s, k1, FALSE);
    arma::mat R_s = X_s - Z_s * Lambda_s.t();
    arma::mat Z_u_s = ComputeSVD(R_s, k2, FALSE);
    arma::mat Lambda_u_s = R_s.t() * Z_u_s * inv(Z_u_s.t() * Z_u_s);

    Lambda.submat(0,   0,     n_c-1,      k1-1)       = Lambda_c;
    Lambda.submat(n_c, 0,     n_c+n_s-1,  k1-1)       = Lambda_s;
    Lambda.submat(0,   k1,    n_c-1,      k1+k2-1)    = Lambda_u_c;
    Lambda.submat(n_c, k1+k2, n_c+n_s-1,  k1+2*k2-1)  = Lambda_u_s;

    return Lambda;
}// end func

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
Rcpp::List CAPER_rcpp(const arma::mat& X_c,
                      const arma::mat& X_s,
                      const int k1,
                      const int k2,
                      const double sigma_c,
                      const double sigma_s,
                      const int max_iter,
                      const double epsilon,
                      const bool verbose = false){
    try {
        if (X_c.n_rows != X_s.n_rows) {
            Rcpp::stop("X_c and X_s must have the same number of rows.");
        }
        if (k1 <= 0 || k2 <= 0) {
            Rcpp::stop("k1 and k2 must be positive.");
        }        
        
        if (verbose) {
          Rcpp::Rcout << "sigma_c: " << sigma_c << "\n";
          Rcpp::Rcout << "sigma_s: " << sigma_s << "\n";
        }

        const int n_c = X_c.n_cols;
        const int n_s = X_s.n_cols;
        const int G   = X_c.n_rows;

        arma::vec diag_I = Itilde(k1, k2, sigma_c, sigma_s);

        arma::mat Lambda = arma::zeros<arma::mat>(n_c + n_s, k1 + 2*k2);
        arma::mat Lambda_new = arma::zeros<arma::mat>(n_c + n_s, k1 + 2*k2);
        arma::mat I_LambdaT = Lambda.t();

        arma::mat Lambda_c_new   = arma::zeros<arma::mat>(n_c, k1);
        arma::mat Lambda_s_new   = arma::zeros<arma::mat>(n_s, k1);
        arma::mat Lambda_u_c_new = arma::zeros<arma::mat>(n_c, k2);
        arma::mat Lambda_u_s_new = arma::zeros<arma::mat>(n_s, k2);

        if (verbose) Rcpp::Rcout << "## Initializing Lambda ...\n";
        ComputeLambda(X_c, X_s, Lambda, k1, k2, n_c, n_s);
        I_LambdaT = UpdateItildeLambdaT(Lambda, k1, k2, sigma_c, sigma_s);

        arma::mat Lambda_u_c_old = Lambda.submat(0,   k1,    n_c-1,      k1+k2-1);
        arma::mat Lambda_u_s_old = Lambda.submat(n_c, k1+k2, n_c+n_s-1,  k1+2*k2-1);

        arma::mat X_all = join_rows(X_c, X_s);
        arma::mat XtX = X_all.t() * X_all;
        arma::vec diag_XtX = XtX.diag();

        arma::mat X_all_c = X_all.cols(0, n_c-1);
        arma::mat X_all_s = X_all.cols(n_c, n_c+n_s-1);

        arma::vec psi_diag = arma::ones<arma::vec>(n_c + n_s);

        arma::mat E_Z = arma::zeros<arma::mat>(k1 + 2*k2, G);
        arma::mat E_Z_part = arma::zeros<arma::mat>(k1 + 2*k2, n_c + n_s);
        arma::mat Var_Z = arma::zeros<arma::mat>(k1 + 2*k2, k1 + 2*k2);
        arma::mat sum_E_ZZt = arma::zeros<arma::mat>(k1 + 2*k2, k1 + 2*k2);

        arma::mat Xc_E_z   = arma::zeros<arma::mat>(n_c, k1);
        arma::mat Xs_E_z   = arma::zeros<arma::mat>(n_s, k1);
        arma::mat Xc_E_z_c = arma::zeros<arma::mat>(n_c, k2);
        arma::mat Xs_E_z_s = arma::zeros<arma::mat>(n_s, k2);

        if (verbose) Rcpp::Rcout << "## Iterating ...\n";
        for(int iter = 0; iter < max_iter; iter++) {
            if (verbose && iter == 0) {
                Rcpp::Rcout << "## Start EM: max_iter=" << max_iter
                            << ", epsilon=" << epsilon
                            << ", k1=" << k1 << ", k2=" << k2 << "\n";
            }

            if (verbose && iter % 10 == 0) {
                Rcpp::Rcout << "Iteration " << iter << " of " << max_iter << "...\n";
            }
            ComputeExpVar(psi_diag, Lambda, I_LambdaT, diag_I, X_all, E_Z, Var_Z, E_Z_part);

            sum_E_ZZt = (double)G * Var_Z + E_Z * E_Z.t();

            arma::mat E_z   = E_Z.rows(0,       k1-1);
            arma::mat E_z_c = E_Z.rows(k1,      k1+k2-1);
            arma::mat E_z_s = E_Z.rows(k1+k2,   k1+2*k2-1);

            Xc_E_z   = X_all_c.t() * E_z.t();
            Xs_E_z   = X_all_s.t() * E_z.t();
            Xc_E_z_c = X_all_c.t() * E_z_c.t();
            Xs_E_z_s = X_all_s.t() * E_z_s.t();

            Lambda_c_new =
              (Xc_E_z - Lambda_u_c_old * (sum_E_ZZt.submat(0, k1, k1-1, k1+k2-1).t()))
              * inv(sum_E_ZZt.submat(0, 0, k1-1, k1-1));

            Lambda_s_new =
              (Xs_E_z - Lambda_u_s_old * (sum_E_ZZt.submat(0, k1+k2, k1-1, k1+2*k2-1).t()))
              * inv(sum_E_ZZt.submat(0, 0, k1-1, k1-1));

            Lambda_u_c_new =
              (Xc_E_z_c - Lambda_c_new * sum_E_ZZt.submat(0, k1, k1-1, k1+k2-1))
              * inv(sum_E_ZZt.submat(k1, k1, k1+k2-1, k1+k2-1));

            Lambda_u_s_new =
              (Xs_E_z_s - Lambda_s_new * sum_E_ZZt.submat(0, k1+k2, k1-1, k1+2*k2-1))
              * inv(sum_E_ZZt.submat(k1+k2, k1+k2, k1+2*k2-1, k1+2*k2-1));

            Lambda_new.zeros();
            Lambda_new.submat(0,   0,     n_c-1,      k1-1)       = Lambda_c_new;
            Lambda_new.submat(n_c, 0,     n_c+n_s-1,  k1-1)       = Lambda_s_new;
            Lambda_new.submat(0,   k1,    n_c-1,      k1+k2-1)    = Lambda_u_c_new;
            Lambda_new.submat(n_c, k1+k2, n_c+n_s-1,  k1+2*k2-1)  = Lambda_u_s_new;

            arma::mat Apsi = Lambda_new * E_Z_part;
            arma::vec diag_A_XtX = sum(Apsi % XtX.t(), 1);
            psi_diag = (diag_XtX - diag_A_XtX) / (double)G;
            psi_diag.elem(find(psi_diag < 1e-5)).fill(1e-5);

            if ((iter > 100) && (iter % 10 == 0)){
                if (sum(square(vectorise(abs(Lambda * Lambda.t() - Lambda_new * Lambda_new.t())))) < epsilon){
                    if (verbose) Rcpp::Rcout << "## Converged at iter " << iter << "\n";
                    break;
                }
            }

            Lambda = Lambda_new;
            I_LambdaT = UpdateItildeLambdaT(Lambda, k1, k2, sigma_c, sigma_s);

            Lambda_new.zeros();
            Lambda_u_c_old = Lambda_u_c_new;
            Lambda_u_s_old = Lambda_u_s_new;
        }

        return List::create(
            Rcpp::Named("Lambda")   = Lambda,
            Rcpp::Named("E_Z")      = E_Z,
            Rcpp::Named("psi_diag") = psi_diag
        );

    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "C++ exception (unknown reason)..." );
    }

    return Rcpp::List::create(); // not reached; avoids -Wreturn-type
}// end func

