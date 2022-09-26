#ifndef __algorithms_hpp__
#define __algorithms_hpp__

#include <cmath>
#include <vector>
#include <iostream>
#include <armadillo>
#include <assert.h>


arma::mat create_tridiagonal(const arma::vec& a, const arma::vec& d, const arma::vec& e);

arma::mat create_tridiagonal(int N, double a, double d, double e);

arma::mat create_symmetric_tridiagonal(int N, double a, double d);

double max_offdiag_symmetric(const arma::mat& A, int& k, int& l);

int analytical_eigenproblem(const arma::mat A, arma::vec& eigval, arma::mat& eigvec);

void jacobi_rotate(arma::mat &A, arma::mat &R, int k, int l);

void jacobi_eigensolver(arma::mat A, double eps, arma::vec &eigenvalues, arma::mat &eigenvectors, const int maxiter, int &iterations, bool &converged);



#endif 